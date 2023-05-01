#!/bin/bash
#PBS -N Part1_trainingset
#PBS -l ncpus=48,walltime=16:00:00,storage=gdata/if89+gdata/xl04,mem=60GB
#PBS -j oe

###################################################################################
### Part one of Genome Annotation pipeline script, this script filters for high quality transcripts in trinity output and generates the 2 sets of exonerate commands to run

### Setting up the log file
export log=${workingdir}/annotation.log
date | tee >> ${log}
echo -e "Part One of annotation pipeline: Generating training gene set\n====================\n Simplifying headers for target genome and trinity fasta files.\n" >> ${log}

### loading modules
module use /g/data/if89/apps/modulefiles
module load blast/2.11.0 Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 genometools/1.6.2 seqtk/1.3

###################################################################################
### Simplifying fasta headers for augustus, augustus can potentially fail with complex fasta headers

# simplifying genome headers
cd ${workingdir}/Simple_header/genome
cp ${targetgenome} .
simplifyFastaHeaders.pl ${targetgenome} SimpleHeaderChr Target_genome_clean.fa header.map
sed -i 's/>//g' header.map

# simplifying trinity headers
cd ${workingdir}/Simple_header/trinity_all
cp ${trinity_out} .
simplifyFastaHeaders.pl ${trinity_out} SimpleTranscript Trinity_all_clean.fa header.map
sed -i 's/>//g' header.map
sed -i 's/ .*//' header.map

###################################################################################
### Running exonerate on cDNA
echo -e " "Running exonerate on trinity output."\n [CMD LOG]\t"for i in ${trinity_out}';' do inputfasta='$(realpath $i)'';' for c in 1 241 481';' \
do qsub -P ${project_code} -o ${workingdir}/Scripts/workflow \
-v querychunktotal=720,querychunkstart='$c',querychunkend='$((c+239))',outputdir=${workingdir}/Exonerate/cDNA_out,inputfasta='${inputfasta}',targetgenome=${workingdir}/Simple_header/genome/Target_genome_clean.fa runexonerate.sh';' done';' done >> ${log}

echo -e " Exonerate job IDs:" >> ${log}
cd ${workingdir}/Scripts
for i in ${trinity_out}; do inputfasta=$(realpath $i); for c in 1 241 481; do qsub -P ${project_code} -o ${workingdir}/Scripts/workflow -v querychunktotal=720,querychunkstart=$c,querychunkend=$((c+239)),outputdir=${workingdir}/Exonerate/cDNA_out,inputfasta=${inputfasta},targetgenome=${workingdir}/Simple_header/genome/Target_genome_clean.fa runexonerate.sh; done; done | awk '{print " [JOB ID]\t"$1}' | tee >> ${log}

###################################################################################
### Filtering the trinity output for high quality transcripts
echo -e "\n Filtering trinity fasta for high quality transcript:\n [1]\tBlast to uniprot and translating output." >> ${log}

# making database
cd ${workingdir}/Trinity_filter/first_filter/database
makeblastdb -in ${workingdir}/uniprot_sprot.fasta \
-parse_seqids \
-title "uniprot_sprot" -dbtype prot \
-out ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot

# blasting trinity.out to uniprot fasta
cd ../blast
blastx -query ${trinity_out} \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-outfmt 6 \
-max_target_seqs 1 \
-num_threads ${PBS_NCPUS} \
-out ${workingdir}/Trinity_filter/first_filter/blast/output.tsv

# translating output
perl ${workingdir}/Scripts/blastxtranslations.pl ${trinity_out} \
${workingdir}/uniprot_sprot.fasta \
${workingdir}/Trinity_filter/first_filter/blast/output.tsv \
Trinity_all

# Start filtering according to blast to uniprot
cd ${workingdir}/Trinity_filter/first_filter/blast

# This is the main filter command below
grep ">" Trinity_all.cds.all.fa | sed 's/>//g' | awk '{print $1"\t"$2"\t"$NF;}' | sed 's/len=//g' | sed 's/:/\t/g' | awk '{ $9 = $4 - $3 } 1' OFS="\t" | awk '$6=="yes"{print ;}' | awk '$7=="yes"{print ;}' | awk '$8 != NA {print ;}' | awk '$8 >= 0.95 {print ;}' | awk '$8 <= 1.05 {print ;}' | sed 's/_i/\t/g' | sort -k10,10 -k3,3nr | awk 'BEGIN{FS=OFS="\t"} {if (a[$1]<$10) {a[$1]=$10; data[$1]=$0}} END{for (i in a) print data[i]}' | awk '$3 <= 10000 {print ;}' | awk '{print $1"_i"$2}' > filter_list.lst
# Let's break this long ass command down into parts and explain each
# Part one, extract the fasta header from the cds.fa, which contains a lot of useful information, we then calculate CDS length from this
# $ grep ">" Trinity_all.cds.all.fa | sed 's/>//g' | awk '{print $1"\t"$2"\t"$NF;}' | sed 's/len=//g' | sed 's/:/\t/g' | awk '{ $9 = $4 - $3 } 1' OFS="\t" > test2
# Part two, this filters for transcripts with start and stop codon, and relative length of >=0.95 and <=1.05  Finally we separate the isoform identifier in the header
# $ awk '$6=="yes"{print ;}' test2 | awk '$7=="yes"{print ;}' | awk '$8 != NA {print ;}' | awk '$8 >= 0.95 {print ;}' | awk '$8 <= 1.05 {print ;}' | sed 's/_i/\t/g' > test3
# Part three, we filter for transcript isoform with the longest CDS (and longest transcript if same longest CDS length)
# for the first sort command, -k10,10 (10 is the column of CDS length) and -k3,3nr (3 is the column of transcript length)
# for the first awk command, $10 is the CDS length, there are only two in the command, do not change other stuff
# no need to change third awk command
# $ sort -k10,10 -k3,3nr| awk 'BEGIN{FS=OFS="\t"} {if (a[$1]<$10) {a[$1]=$10; data[$1]=$0}} END{for (i in a) print data[i]}' | awk '$3 <= 10000 {print ;}' | awk '{print $1"_i"$2}' > test4

# extract list of high quality peptide for all vs all filter
seqtk subseq Trinity_all.pep.all.fa filter_list.lst > ../../second_filter/blast/peptide_filtered.fa

###################################################################################
### ALL VS ALL clustering
echo -e " [2]\tAll vs all blast." >> ${log}

cd ${workingdir}/Trinity_filter/second_filter/database
cp ${workingdir}/Trinity_filter/second_filter/blast/peptide_filtered.fa .

# making database
makeblastdb -in ${workingdir}/Trinity_filter/second_filter/database/peptide_filtered.fa \
-parse_seqids \
-title "peptide_filtered" -dbtype prot \
-out ${workingdir}/Trinity_filter/second_filter/database/peptide_filtered

# running blastp
cd ${workingdir}/Trinity_filter/second_filter/blast
blastp -query ${workingdir}/Trinity_filter/second_filter/blast/peptide_filtered.fa \
-db ${workingdir}/Trinity_filter/second_filter/database/peptide_filtered \
-parse_deflines \
-outfmt 6 \
-num_threads ${PBS_NCPUS} \
-out ${workingdir}/Trinity_filter/second_filter/blast/all_vs_all_out.tsv

### Cluster and filter ALL VS ALL results

# filtering blastp result for e-value < 1e-6 and percentage identity of >= 80%
echo -e " [3]\tRemove all but one transcript from each cluster.\n" >> ${log}

awk '$1 != $2{print ;}' all_vs_all_out.tsv | awk '$3 >= 80{print ;}' | awk '$11 < 1e-6{print ;}' > ${workingdir}/Trinity_filter/third_filter/all_vs_all_out_filtered.tsv
cd ${workingdir}/Trinity_filter/third_filter
awk '{print "1stCOLUMN\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' all_vs_all_out_filtered.tsv > all_vs_all_out_filtered2.tsv

# now we group the blast result into clusters of highly similar proteins
cat all_vs_all_out_filtered2.tsv | perl -lne '$_=~s/\"//g; @a=split("\t",$_); $g{$a[1]}{$a[1]}=""; $g{$a[2]}{$a[2]}=""; $g{$a[1]}{$a[2]}=""; $g{$a[2]}{$a[1]}=""; END {foreach $v (sort keys %g) { print join("\t", map { $_ } sort keys %{$g{$v}})}}' | sort | uniq > cluster.txt

# removing the first entry in each line (removing one transcript in each cluster, this is the transcript we keep)
perl -pe 's/.*?\t//' cluster.txt > cluster2.txt

# turning all transcripts in cluster2 into a list
tr '\t' '\n' < cluster2.txt > cluster3.txt

# removing duplicated entries, now cluster4.txt contains the list of transcript we remove from the high-quality transcript list
sort cluster3.txt | uniq -u > cluster4.txt
grep -v -w -f cluster4.txt ${workingdir}/Trinity_filter/first_filter/blast/filter_list.lst > ${workingdir}/Trinity_filter/results/filter_list_final.lst
cd ${workingdir}/Trinity_filter/results
seqtk subseq ${workingdir}/Trinity_filter/first_filter/blast/Trinity_all.cdna.all.fa filter_list_final.lst > cDNA_filtered_final.fa
seqtk subseq ${workingdir}/Trinity_filter/first_filter/blast/Trinity_all.cds.all.fa filter_list_final.lst > CDS_filtered_final.fa
seqtk subseq ${workingdir}/Trinity_filter/first_filter/blast/Trinity_all.pep.all.fa filter_list_final.lst > peptide_filtered_final.fa

###################################################################################
### Running exonerate on filtered CDS
echo -e " "Running exonerate on filtered CDS."\n [CMD LOG]\t"for i in ${workingdir}/Trinity_filter/results/CDS_filtered_final.fa';' do inputfasta='$(realpath $i)'';' for c in 1 241 481';' \
do qsub -P ${project_code}-o ${workingdir}/Scripts/workflow \
-v querychunktotal=720,querychunkstart='$c',querychunkend='$((c+239))',outputdir=${workingdir}/Exonerate/CDS_filtered_out,inputfasta='${inputfasta}',targetgenome=${workingdir}/Simple_header/genome/Target_genome_clean.fa runexonerate.sh';' done';' done >> ${log}

echo -e " Exonerate job IDs:" >> ${log}
cd ${workingdir}/Scripts
for i in ${workingdir}/Trinity_filter/results/CDS_filtered_final.fa; do inputfasta=$(realpath $i); for c in 1 241 481; do qsub -P ${project_code} -o ${workingdir}/Scripts/workflow -v querychunktotal=720,querychunkstart=$c,querychunkend=$((c+239)),outputdir=${workingdir}/Exonerate/CDS_filtered_out,inputfasta=${inputfasta},targetgenome=${workingdir}/Simple_header/genome/Target_genome_clean.fa runexonerate.sh; done; done | awk '{print " [JOB ID]\t"$1}' | tee >> ${log}

###################################################################################
### Setting up commands to run for part 2 in log file
time=$(date)
echo -e "\n #################### ${time}\n Part One completed, Run part Two with the following commands after all 6 exonerate jobs finish successfully.\n CMDs to run after your exonerate job finishes:\n [!]\tcd ${workingdir}/Scripts/workflow\n [!]\tqsub -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part2_trainingset.sh" >> ${log}
