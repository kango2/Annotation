#!/bin/bash
#PBS -N Part2_trainingset
#PBS -l ncpus=1,walltime=4:00:00,storage=gdata/if89+gdata/xl04,mem=80GB
#PBS -j oe

###################################################################################
### Part two of Genome Annotation pipeline script

### Setting up the log file
export log=${workingdir}/annotation.log
echo -e "\n\n\n" >> ${log}
date | tee >> ${log}
echo -e "Part Two of annotation pipeline: Generating training gene set\n====================\n Merging parallel output of exonerate." >> ${log}

### Loading modules
module use /g/data/if89/apps/modulefiles
module load blast/2.11.0 Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 genometools/1.6.2 bedtools/2.28.0

###################################################################################
### Process the exonerate output
echo -e " [CMD LOG]\t"perl processexonerate.pl -inputfasta ${trinity_out} -outputdir ${workingdir}/Exonerate/cDNA_out -targetgenome ${workingdir}/Simple_header/genome/Target_genome_clean.fa -querychunktotal 720 -hintpriority 4 -hintsource E >> ${log}
echo -e " [CMD LOG]\t"perl processexonerate.pl -inputfasta ${workingdir}/Trinity_filter/results/CDS_filtered_final.fa -outputdir ${workingdir}/Exonerate/CDS_filtered_out -targetgenome ${workingdir}/Simple_header/genome/Target_genome_clean.fa -querychunktotal 720 -hintpriority 4 -hintsource E >> ${log}

cd ${workingdir}/Scripts
perl processexonerate.pl -inputfasta ${trinity_out} -outputdir ${workingdir}/Exonerate/cDNA_out -targetgenome ${workingdir}/Simple_header/genome/Target_genome_clean.fa -querychunktotal 720 -hintpriority 4 -hintsource E
perl processexonerate.pl -inputfasta ${workingdir}/Trinity_filter/results/CDS_filtered_final.fa -outputdir ${workingdir}/Exonerate/CDS_filtered_out -targetgenome ${workingdir}/Simple_header/genome/Target_genome_clean.fa -querychunktotal 720 -hintpriority 4 -hintsource E

###################################################################################
### Processing and filtering exonerate output for duplicate alignments

# Exonerate can sometime have multiple alignments for the same transcript (even with bestn=1 due to identical alignment score), therefore we need to filter them out
export TRINITY_BASE=$(basename ${trinity_out} | sed 's![^.]*$!!')
cd ${workingdir}/Exonerate/annotation/working_dir
cp ${workingdir}/Exonerate/cDNA_out/${TRINITY_BASE}exonerate.target.genes.gff3 .
cp ${workingdir}/Exonerate/CDS_filtered_out/CDS_filtered_final.fa.exonerate.target.genes.gff3 .

### Filter out duplicate alignments
echo -e "\n "Filtering out transcripts with multiple alignments due to identical exonerate alignment score. >> ${log}

# Remove from CDS gff3 file
echo -e " [1]\t"Removing from CDS gff3. >> ${log}

grep "gene" CDS_filtered_final.fa.exonerate.target.genes.gff3 | cut -f9 | cut -f4 -d ";" | sed 's/from=//g' | sort | uniq -d | sort > CDS_duplicated
grep -v -w -f CDS_duplicated CDS_filtered_final.fa.exonerate.target.genes.gff3 > CDS_filtered_final.fa.exonerate.target.genes2.gff3

# remove from cDNA gff3 file
echo -e " [2]\t"Removing from cDNA gff3. >> ${log}

grep "gene" ${TRINITY_BASE}exonerate.target.genes.gff3 | cut -f9 | cut -f4 -d ";" | sed 's/from=//g' | sort | uniq -d | sort > cDNA_duplicated
grep -v -w -f cDNA_duplicated ${TRINITY_BASE}exonerate.target.genes.gff3 > ${TRINITY_BASE}exonerate.target.genes2.gff3

# also filter them out in the hints file
echo -e " [3]\t"Removing from hints file. >> ${log}

grep -v -w -f cDNA_duplicated ${workingdir}/Exonerate/cDNA_out/${TRINITY_BASE}exonerate.target.hints.gff3 > ${workingdir}/Exonerate/cDNA_out/${TRINITY_BASE}exonerate.target.hints_unique.gff3

# move and rename hints file to a more identifiable place
mv ${workingdir}/Exonerate/cDNA_out/${TRINITY_BASE}exonerate.target.hints_unique.gff3 ${workingdir}/Augustus/hints/exonerate_hints.gff3

###################################################################################
### Now modify the cDNA and CDS gff3 output from exonerate to generate a training set for training augustus
echo -e "\n "Generating training set from cDNA and CDS exonerate gff3. >> ${log}

# replacing exon to CDS in 3rd column | #replacing .e to .cds in 9th column | #selecting only the CDS
awk '$3=="exon"{$3="CDS"} 1' OFS="\t"   CDS_filtered_final.fa.exonerate.target.genes2.gff3 | sed 's/\.e/\.cds/g' | grep CDS > CDS.gff3

# selecting genes that are in both exon and CDS gff3
grep -w "gene" ${TRINITY_BASE}exonerate.target.genes2.gff3 | awk -F '\t' '{gsub(/.*from=/,"",$9); gsub(/;/,"",$9); print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' > Trinity_all.bed
grep -w "gene" CDS_filtered_final.fa.exonerate.target.genes2.gff3 | awk -F '\t' '{gsub(/.*from=/,"",$9); gsub(/;/,"",$9); print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' > CDS.bed

# -s so it's strand specific, -wa and -wb to output both input bed files.
echo -e " [1]\t"Selecting genes in both cDNA gff3 and CDS gff3, and are overlapping. >> ${log}

bedtools intersect -a Trinity_all.bed -b CDS.bed -wa -wb -s > trinity_all_vs_CDS
awk '$4 == $10' trinity_all_vs_CDS | cut -f4 > both.txt
echo '##gff-version 3' > merge.gff3
grep -w -f both.txt ${TRINITY_BASE}exonerate.target.genes2.gff3 >> merge.gff3
grep -w -f both.txt CDS.gff3 >> merge.gff3
gt gff3 -sortlines -tidy -retainids merge.gff3 > merge2.gff3
gt gff3 -retainids merge2.gff3 > merge3.gff3

### this script adds the inferred UTRs from exon-CDS differences
echo -e " [2]\t"Adding inferred UTRs from exon-CDS difference. >> ${log}

python3 ${workingdir}/Scripts/addUTRs.py ${workingdir}/Exonerate/annotation/working_dir/merge3.gff3 ${workingdir}/Exonerate/annotation/working_dir/merge3_with_UTRs.gff3
sed -i 's/exonerate:est2genome/est2genome/g' merge3_with_UTRs.gff3
grep mRNA merge3_with_UTRs.gff3 | cut -f9 | cut -d ';' -f1 | sed 's/ID=//g' > IDs_in_gff.txt
grep -w -f IDs_in_gff.txt ${workingdir}/Simple_header/trinity_all/header.map > header_filtered.map
awk '{print "s/"$2"\\b/"$1"/g"}' header_filtered.map > header_filtered2.map
sed -i -f header_filtered2.map merge3_with_UTRs.gff3

# grep-ing the relevant gff3 information from annotation file needed for training augustus
grep -v -E 'exon|gene|mRNA' merge3_with_UTRs.gff3 > ${workingdir}/Exonerate/annotation/CDS_with_UTR/CDS_with_UTRs.gff3
grep -v -E 'exon|UTR|gene|mRNA' merge3_with_UTRs.gff3 > ${workingdir}/Exonerate/annotation/CDS_only/CDS.gff3

# we now select only genes that have both 3' and 5' utrs for training augustus
echo -e " [3]\t"Selecting genes with both 3\' and 5\' UTRs. >> ${log}

grep "three_prime_UTR" merge3_with_UTRs.gff3 | cut -f9 | cut -d "." -f1 | sort | uniq > all3
grep "five_prime_UTR" merge3_with_UTRs.gff3 | cut -f9 | cut -d "." -f1 | sort | uniq > all5
comm -12 <(sort all3) <(sort all5) > bothutr.lst
sed -i "s/ID=//g" bothutr.lst

### Turning gffs into genbank format with 2000 bases of flanking region
echo -e " [4]\t"Turning GFF3 into genbank format with 2000 bases of flanking region. >> ${log}

# gff3 to genbank format
cd ${workingdir}/Exonerate/annotation/CDS_only
gff2gbSmallDNA.pl CDS.gff3 ${workingdir}/Simple_header/genome/Target_genome_clean.fa 2000 training.gb --good=${workingdir}/Exonerate/annotation/working_dir/bothutr.lst
cd ${workingdir}/Exonerate/annotation/CDS_with_UTR
gff2gbSmallDNA.pl CDS_with_UTRs.gff3 ${workingdir}/Simple_header/genome/Target_genome_clean.fa 2000 training_with_utr.gb --good=${workingdir}/Exonerate/annotation/working_dir/bothutr.lst

### randomly choose 500 for training
echo -e " [5]\tRandomly choose 500 genes for training Augustus. Note that all genes NOT chosen will be used for an initial training step, the chosen 500 will be used for optimising Augustus species-specific parameters."  >> ${log}

# this generates training.gb.test with 500 genes
# and training.gb.train with the rest ~4500 genes
cd ${workingdir}/Exonerate/annotation/CDS_only
randomSplit.pl training.gb 500

# we then further split the training.gb.test (contains 500 genes) into 200 for evaluation and training, whilst the other 300 for training only
# evaluation is the step that takes very long whereas training is fast, but since we now have a way to use the parallel perl module we can use more
echo -e " [6]\tRandomly split the 500 genes into 200 for evaluation + training, and 300 for training only. This is because evaluation takes long but not training."  >> ${log}

randomSplit.pl training.gb.test 200

# similarily for generating UTR training
cd ${workingdir}/Exonerate/annotation/CDS_with_UTR
randomSplit.pl training_with_utr.gb 500
randomSplit.pl training_with_utr.gb.test 200

###################################################################################
### Setting up commands to run for part 3 in log file
echo -e "\n Path to folders containing the training gene set files, please use the following commands to check if the files are empty.\n [!]\tls -shl ${workingdir}/Exonerate/annotation/CDS_only\n [!]\tls -shl ${workingdir}/Exonerate/annotation/CDS_with_UTR"  >> ${log}
time=$(date)
echo -e "\n #################### ${time}\n Part Two completed, please check and confirm that your training gene set file is not empty, and then run part Three with the following commands.\n CMDs to run:\n [!]\tcd ${workingdir}/Scripts/workflow\n [!]\tqsub -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part3_TrainingAugustus.sh" >> ${log}
