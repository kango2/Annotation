#!/bin/bash
#PBS -N Part4_RunningAugustus
#PBS -l ncpus=16,walltime=24:00:00,storage=gdata/if89+gdata/xl04,mem=60GB,jobfs=60GB
#PBS -j oe

###################################################################################
### Part four of Genome Annotation pipeline script

### Setting up the log file
export log=${workingdir}/annotation.log
echo -e "\n\n\n" >> ${log}
date | tee >> ${log}
echo -e "Part Four of annotation pipeline: Running Augustus for your species\n====================\n Running Augustus gene prediction in ${PBS_NCPUS} chunks." >> ${log}

### Loading modules
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022

### Setting up augustus path
export AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config

###################################################################################
### Running augustus in parallel and processing output

# This generates the augustus commands
cd ${workingdir}/Augustus/split_genome
for i in $(ls -1v Target_genome_clean.fa_chunk* | sed 's/Target_genome_clean.fa_chunk_//g'); do echo augustus --species=${species} --softmasking=1 --UTR=on --print_utr=on --hintsfile=${workingdir}/Augustus/hints/exonerate_hints.gff3 --extrinsicCfgFile=${workingdir}/Augustus/config/extrinsic/extrinsic.MPE.cfg --exonnames=off --stopCodonExcludedFromCDS=t --AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config ${workingdir}/Augustus/split_genome/Target_genome_clean.fa_chunk_"$i" ">" ${workingdir}/Augustus/annotation/aug"$i".out; done | awk '{print " [CMD LOG]\t"$0}' >> ${log}
for i in $(ls -1v Target_genome_clean.fa_chunk* | sed 's/Target_genome_clean.fa_chunk_//g'); do echo augustus --species=${species} --softmasking=1 --UTR=on --print_utr=on --hintsfile=${workingdir}/Augustus/hints/exonerate_hints.gff3 --extrinsicCfgFile=${workingdir}/Augustus/config/extrinsic/extrinsic.MPE.cfg --exonnames=off --stopCodonExcludedFromCDS=t --AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config ${workingdir}/Augustus/split_genome/Target_genome_clean.fa_chunk_"$i" ">" ${workingdir}/Augustus/annotation/aug"$i".out; done | parallel --jobs ${PBS_NCPUS} {}

# This merges the output and process them
echo -e "\n Merging results from running Augustus in parallel." >> ${log}

for i in $(ls -1v Target_genome_clean.fa_chunk* | sed 's/Target_genome_clean.fa_chunk_//g'); do cat ${workingdir}/Augustus/annotation/aug"$i".out >> ${workingdir}/Augustus/annotation/aug.out; done
cd ${workingdir}/Augustus/annotation
cat aug.out | join_aug_pred.pl > augustus.gff
getAnnoFasta.pl augustus.gff
grep -v "#" augustus.gff > augustus.gtf

# This changes the contig name in augustus.gtf back to their original name (they were changed to a cleaner name whilst running augustus)
grep -w -f <(cut -f1 augustus.gtf | sort | uniq) ${workingdir}/Simple_header/genome/header.map > filtered_header
awk '{print "s/"$1"\\b/"$2"/g"}' filtered_header > filtered_header2
sed -i -f filtered_header2 augustus.gtf

# This outputs a gff3 file from the gtf file
echo -e "\n Annotation is now ready in ${workingdir}/Augustus/annotation/annotation.gff3" >> ${log}

gtf2gff.pl < augustus.gtf --out=augustus.gff3 --gff3

# This prepares for bedintersect which is needed for stitching
grep -w 'CDS\|three_prime_utr\|five_prime_utr' augustus.gff3 > augustus_for_bedintersect.gff3

###################################################################################
### Setting up commands to run for Part 5A and 5B in the log file, user will only have to choose one
time=$(date)
echo -e "\n #################### ${time}\n Part Four completed, Lastly run\n Either part FiveA OR FiveB, explanation below.\n\n Part5A will stitch together some augustus genes based on overlapping with cDNA alignments to the ${species} genome from a related species, then blast all genes to uniprot database for identification.\n NOTE that stitching only works with cDNA fasta file from ENSEMBL currently.\n NOTE that you may want to set a custom threshold of max intron length for the stitched gene, stitched gene containing an inferred intron with length greater than your value here will be discarded.\n leaving it as 0 will stitch together every candidates.\n CMDs to run for Part5A, change path in first command:\n [!]\texport Related_species=/path/to/related_species_cDNA.fa\n [!]\texport maxintron=0\n [!]\tcd ${workingdir}/Scripts/workflow\n [!]\tqsub -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code},Related_species=\${Related_species},maxintron=\${maxintron} Part5A_StitchThenIdentify.sh\n\n Part5B will just identify the Augustus genes by blasting to uniprot database.\n CMDs to run for Part5B:\n [!]\tcd ${workingdir}/Scripts/workflow\n [!]\tqsub -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part5B_Identify.sh" >> ${log}
