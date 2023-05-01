#!/bin/bash
#PBS -N Part5A_StitchThenIdentify
#PBS -l ncpus=32,walltime=12:00:00,storage=gdata/if89+gdata/xl04,mem=80GB,jobfs=80GB
#PBS -j oe

###################################################################################
### Part five A of Genome Annotation pipeline script

### Setting up the log file
export log=${workingdir}/annotation.log
echo -e "\n\n\n" >> ${log}
date | tee >> ${log}
echo -e "Part FiveA of annotation pipeline: Stitching Augustus genes and Similarity search to Uniprot database (BLAST)\n====================\n Blasting cDNA from related species to ${species} genome." >> ${log}

### Loading modules
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 blast/2.11.0 genometools/1.6.2 seqkit/2.3.1 gffread/0.12.7 bedtools/2.28.0 seqtk/1.3

###################################################################################
### This script prepares the gene candidates for stitching

# First, blast related_species cDNA to genome

# making database
cd ${workingdir}/Augustus/annotation_stitch/blast/database
makeblastdb -in ${targetgenome} \
-parse_seqids \
-title "genome" -dbtype nucl \
-out ${workingdir}/Augustus/annotation_stitch/blast/database/genome

# blast related_species cDNA to genome
blastn -query ${Related_species} \
-db ${workingdir}/Augustus/annotation_stitch/blast/database/genome \
-max_target_seqs 1 \
-evalue 1e-6 \
-parse_deflines \
-outfmt 6 \
-num_threads ${PBS_NCPUS} \
-out ${workingdir}/Augustus/annotation_stitch/blast/related_cDNA_to_genome.tsv

# Add gene id immediately following the transcript id, separated by "|" , then turns the blast tsv output into a bed file
echo -e " Turning blast output into bed format." >> ${log}

cd ${workingdir}/Augustus/annotation_stitch/blast
join <(sort related_cDNA_to_genome.tsv) <(grep ">" ${Related_species} | cut -f1,4 -d' ' | sed 's/ /\t/g' | sed 's/>//g' | sort) -t $'\t' -o1.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14 | perl -pe 's/\t/|/' | perl -lne '@a=split("\t", $_); print "$a[1]\t". (($a[8] < $a[9]) ? $a[8] : $a[9] )."\t".(($a[8] < $a[9]) ? $a[9] : $a[8])."\t$a[0]\t$a[11]\t".(($a[8] < $a[9]) ? "+" : "-" )' > related_cDNA_to_genome.bed
mv related_cDNA_to_genome.bed ${workingdir}/Augustus/annotation_stitch/intersect_and_stitch/related_cDNA_to_genome.bed

# Generating the list of genes to be stitched together
echo -e "\n Generating list of genes to be stitched together." >> ${log}

cd ${workingdir}/Augustus/annotation_stitch/intersect_and_stitch
awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' ../../annotation/augustus_for_bedintersect.gff3 > augustus.bed

# -s so it's strand specific, -wa and -wb to output both input bed files.
bedtools intersect -a related_cDNA_to_genome.bed -b augustus.bed -wa -wb -s > intersect
cut -f4,10 intersect | sed 's/;/\t/g' | cut -f1,3 | sed 's/Parent=//g' | sed 's/\.t1//g' | sort | uniq | sort | awk '{ printf "%s", (NR==1 || pre!=$1? (NR>1? ORS:"")$1: "") OFS $2; pre=$1 }
    END  { print "" }' | sed 's/ /\t/g' | awk 'NF>=3' > intersect2
cut -f2- -d "|" intersect2 | sort | uniq > intersect3.tsv

# Since the annotation file are in order, the IDs we grep will also be in order, this generates a long list of .tsv file of consecutive IDs
grep gene ../../annotation/augustus.gff3 | grep "+" | cut -f9 | sed 's/;//g' | sed 's/ID=//g' | paste -s -d "\t" > consecutive_positive.tsv
grep gene ../../annotation/augustus.gff3 | grep "-" | cut -f9 | sed 's/;//g' | sed 's/ID=//g' | paste -s -d "\t" > consecutive_negative.tsv

# We now try grep-ing each intersect cluster from the consecutive IDs, cluster contains neighbouring genes (consecutive IDs) if it can be grep-ed from the consecutive IDs

# we then filter for clusters that can be grep-ed (with yes)
cut -f2- intersect3.tsv | while read a; do if grep -q "$a" consecutive_positive.tsv; then echo -e yes'\t'"$a"; else echo -e no'\t'"$a"; fi done > intersect4_pos.tsv
grep yes intersect4_pos.tsv | cut -f2- | sort | uniq > intersect5_pos.tsv
cat intersect5_pos.tsv | while read -r a; do OCCUR=$(grep "${a}" intersect5_pos.tsv | wc -l); echo -e $OCCUR'\t'"${a}"; done | awk '$1 == 1 {print ;}' | cut -f2- > intersect6_pos.tsv

cut -f2- intersect3.tsv | while read a; do if grep -q "$a" consecutive_negative.tsv; then echo -e yes'\t'"$a"; else echo -e no'\t'"$a"; fi done > intersect4_neg.tsv
grep yes intersect4_neg.tsv | cut -f2- | sort | uniq > intersect5_neg.tsv
cat intersect5_neg.tsv | while read -r a; do OCCUR=$(grep "${a}" intersect5_neg.tsv | wc -l); echo -e $OCCUR'\t'"${a}"; done | awk '$1 == 1 {print ;}' | cut -f2- > intersect6_neg.tsv

###################################################################################
### Stitching together genes in the candidate list

# This generates a "temp_all" file with all the stitched genes

# stitching all the clusters on positive strand
echo -e " Stitching together genes on positive strand." >> ${log}
StitchedGeneCount=1
cd ${workingdir}/Augustus/annotation_stitch
cat intersect_and_stitch/intersect6_pos.tsv | while read a; do last=$(echo $a | wc -w); one=1; lastone=$[last - one]; echo "$a" > temp_input; perl -i -pe "chomp if eof" temp_input; readarray -d $'\t' -t a < temp_input; for index in "${!a[@]}"; do if [ $index == 0 ]; then    \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'gene\|mRNA\|transcription_start_site\|five_prime_utr\|start_codon\|CDS\|intron' > temp); \
elif [ $index == $lastone ]; then \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'transcription_end_site\|stop_codon\|three_prime_utr\|CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
else \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
fi done; gt gff3 -sort -tidy -retainids -addintrons <(sed 's/CDS/exon/g' temp) > temp2; grep "$(cut -f3,4,5 temp2 | grep intron | sort | uniq -u)" temp2 >> temp; gt gff3 -sort -tidy -retainids temp | grep -v "#" > temp2; endcoord=$(sort -n -r -k 5 temp2 | cut -f5 | head -n 1); CurrentGeneID="StitchedPg"${StitchedGeneCount}; \
perl ${workingdir}/Scripts/AppendStitchedGene.pl ${workingdir}/Augustus/annotation_stitch/temp ${CurrentGeneID} | awk '$2=="."{$2="Inferred"} 1' OFS="\t" | awk '$3=="gene"{$5='$endcoord'} 1' OFS="\t" | awk '$3=="mRNA"{$5='$endcoord'} 1' OFS="\t" > temp3; \
output=$(awk '$2 == "Inferred"' temp3 | awk -v maxintron="$maxintron" '$5-$4 <= maxintron || maxintron == 0 {print}'); if [ -n "$output" ]; \
then cat temp3 >> temp_all; echo -e "$(cat temp_input)" | sed 's/\t/,/g' | awk -v count="${CurrentGeneID}" '{print count"\t"$1}' >> sourcegene; echo -e "$(cat temp_input)" >> sourcegene2; ((StitchedGeneCount++)); \
else grep -w -f <(tr '\t' '\n' < temp_input) ../annotation/augustus.gff3 >> temp_all; \
fi; \
rm temp; \
rm temp2; \
rm temp3; done

# stitching all the clusters on negative strand
echo -e " Stitching together genes on negative strand." >> ${log}
StitchedGeneCount=1
cd ${workingdir}/Augustus/annotation_stitch
cat intersect_and_stitch/intersect6_neg.tsv | while read a; do last=$(echo $a | wc -w); one=1; lastone=$[last - one]; echo "$a" > temp_input; perl -i -pe "chomp if eof" temp_input; readarray -d $'\t' -t a < temp_input; for index in "${!a[@]}"; do if [ $index == 0 ]; then    \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'gene\|mRNA\|transcription_end_site\|three_prime_utr\|stop_codon\|CDS\|intron' > temp); \
elif [ $index == $lastone ]; then \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'transcription_start_site\|start_codon\|five_prime_utr\|CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
else \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
fi done; gt gff3 -sort -tidy -retainids -addintrons <(sed 's/CDS/exon/g' temp) > temp2; grep "$(cut -f3,4,5 temp2 | grep intron | sort | uniq -u)" temp2 >> temp; gt gff3 -sort -tidy -retainids temp | grep -v "#" > temp2; endcoord=$(sort -n -r -k 5 temp2 | cut -f5 | head -n 1); CurrentGeneID="StitchedNg"${StitchedGeneCount}; \
perl ${workingdir}/Scripts/AppendStitchedGene.pl ${workingdir}/Augustus/annotation_stitch/temp ${CurrentGeneID} | awk '$2=="."{$2="Inferred"} 1' OFS="\t" | awk '$3=="gene"{$5='$endcoord'} 1' OFS="\t" | awk '$3=="mRNA"{$5='$endcoord'} 1' OFS="\t" > temp3; \
output=$(awk '$2 == "Inferred"' temp3 | awk -v maxintron="$maxintron" '$5-$4 <= maxintron || maxintron == 0 {print}'); if [ -n "$output" ]; \
then cat temp3 >> temp_all; echo -e "$(cat temp_input)" | sed 's/\t/,/g' | awk -v count="${CurrentGeneID}" '{print count"\t"$1}' >> sourcegene; echo -e "$(cat temp_input)" >> sourcegene2; ((StitchedGeneCount++)); \
else grep -w -f <(tr '\t' '\n' < temp_input) ../annotation/augustus.gff3 >> temp_all; \
fi; \
rm temp; \
rm temp2; \
rm temp3; done

rm temp_input

tr '\t' '\n' < sourcegene2 | sort | uniq > sourcegene3

###################################################################################
### Merging stitched genes into the original annotation

# Get rid of all those single gene candidates before merging together
tr '\t' '\n' < intersect_and_stitch/intersect6_pos.tsv > single.lst
tr '\t' '\n' < intersect_and_stitch/intersect6_neg.tsv >> single.lst
sort single.lst | uniq > single2.lst
grep -v -w -f single2.lst ../annotation/augustus.gff3 > augustus2.gff3

# Now we merge augustus2.gff3 (which does not have those single genes), and temp_all (which contains those single genes but stitched together), then sort them with genometools
gt gff3 -retainids -tidy <(cat augustus2.gff3 temp_all) > augustus_stitched.gff3

### Generating protein fasta for the stitched annotation
echo -e "\n Generating protein fasta file from final annotation." >> ${log}

gffread -y temporary.fa -g ${targetgenome} augustus_stitched.gff3
seqtk seq -l 60 temporary.fa | seqkit sort -N > augustus_stitched_protein.fa

# gffread uses full stop to represent ambiguous amino acids, which is not ideal, changing it to the standard "X"
sed -i '/^[^>]/s/\./X/g' augustus_stitched_protein.fa
rm temporary.fa

# copying this final protein fasta file to workingdir
echo -e " Copying protein fasta to working directory." >> ${log}

cp ${workingdir}/Augustus/annotation_stitch/augustus_stitched_protein.fa ${workingdir}/Augustus_stitched_proteins.fasta

###################################################################################
### Structural identification of the predicted genes by similarity search to uniprot proteins
echo -e "\n Blasting Augustus proteins to uniprot database." >> ${log}

# Blasting augustus proteins to uniprot, printing the full fasta header in uniprot DB
cd ${workingdir}/Augustus/annotation_functional
blastp -query ${workingdir}/Augustus/annotation_stitch/augustus_stitched_protein.fa \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-evalue 1e-3 \
-outfmt "6 std salltitles" \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads ${PBS_NCPUS} \
-out ${workingdir}/Augustus/annotation_functional/blast.tsv

seqtk subseq -l 60 ${workingdir}/Augustus/annotation/augustus.aa <(cat ${workingdir}/Augustus/annotation_stitch/sourcegene3 | awk '{print $0".t1"}') > sourceProt.aa
blastp -query ${workingdir}/Augustus/annotation_functional/sourceProt.aa \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-evalue 1e-3 \
-outfmt "6 std salltitles" \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads ${PBS_NCPUS} \
-out ${workingdir}/Augustus/annotation_functional/blast2.tsv

### concatenating the blast output to the final annotation as two extra attributes on 9th column of "gene" features, separated by ";"

# if a protein does not have gene name in the header then it prints unknown instead
echo -e " Concatenating blast results into annotation gff3 file and copying to working directory." >> ${log}

join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(grep -w "gene" ${workingdir}/Augustus/annotation_stitch/augustus_stitched.gff3 | cut -f9 | sed 's/ID=//g' | awk '{print $0".t1"}' | sort -k1,1) <(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); print $1"\t"$2"\t" b[1]} else {print $1"\t"$2"\tunknown"}}' ${workingdir}/Augustus/annotation_functional/blast.tsv | sort -k1,1) | sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' > geneID_to_uniprotID.tabular
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9";uniprot_ID="a[$9]}1' geneID_to_uniprotID.tabular ${workingdir}/Augustus/annotation_stitch/augustus_stitched.gff3 > ${workingdir}/Augustus_stitched_annotation.gff3

join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(grep -w -f ${workingdir}/Augustus/annotation_stitch/sourcegene3 ${workingdir}/Augustus/annotation/augustus.gff3 | grep -w "gene" | cut -f9 | sed 's/ID=//g' | sed 's/;//g' | awk '{print $0".t1"}' | sort -k1,1) <(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); print $1"\t"$2"\t" b[1]} else {print $1"\t"$2"\tunknown"}}' ${workingdir}/Augustus/annotation_functional/blast2.tsv | sort -k1,1) | sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' > sourceID_to_uniprotID.tabular
gt gff3 -retainids -tidy <(grep -w -f ${workingdir}/Augustus/annotation_stitch/sourcegene3 ${workingdir}/Augustus/annotation/augustus.gff3) > sourceGene.gff3
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9";uniprot_ID="a[$9]}1' sourceID_to_uniprotID.tabular ${workingdir}/Augustus/annotation_functional/sourceGene.gff3 > ${workingdir}/Augustus_source_for_stitching.lst

### generating a tabular file with gene info
echo -e "\n Generating a tabular file with some basic gene information." >> ${log}

cd ${workingdir}
echo -e "Contig\tMiddle position\tGene length\tStrand\tGene ID\tUniprot accession number\tGene name" > ${workingdir}/Augustus_gene_table.tabular

# b[1] is geneid, d[1] is uniprot AC, f[1] is gene name
awk -F'\t' '$3 == "gene" {midpoint = int(($4 + $5) / 2); split($9,a,"ID="); split(a[2],b,";"); split($9,c,"uniprot_ID="); split(c[2],d,";"); split($9,e,"gene_name="); split(e[2],f," "); print $1"\t"midpoint"\t"$5-$4"\t"$7"\t"b[1]"\t"d[1]"\t"f[1]}' ${workingdir}/Augustus_stitched_annotation.gff3 >> ${workingdir}/Augustus_gene_table.tabular

### generating a list of stitched genes and their source genes
cat ${workingdir}/Augustus/annotation_stitch/sourcegene | cat - ${workingdir}/Augustus_source_for_stitching.lst > temp && mv temp ${workingdir}/Augustus_source_for_stitching.lst
echo -e "Stitched gene ID\tSource genes (comma separated)" | cat - ${workingdir}/Augustus_source_for_stitching.lst > temp && mv temp ${workingdir}/Augustus_source_for_stitching.lst
cp ${workingdir}/Augustus/annotation_functional/sourceProt.aa ${workingdir}/Augustus_source_for_stitching.fasta

###################################################################################
### Annotation complete
time=$(date)
echo -e "\n #################### ${time}\n Part FiveA completed.\n\n Annotation pipeline completed, you can find your output files in:\n\t${workingdir}" >> ${log}
