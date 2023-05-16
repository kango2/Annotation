#!/bin/bash
#PBS -N OptionA_Stitch
#PBS -l ncpus=16,walltime=32:00:00,storage=gdata/if89+gdata/xl04,mem=80GB
#PBS -j oe

###################################################################################
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 blast/2.11.0 genometools/1.6.2 seqkit/2.3.1 gffread/0.12.7 bedtools/2.28.0 seqtk/1.3 exonerate/2.2.0

cd ${workingdir}/Scripts
perl processexonerate.pl -inputfasta ${trinity_out} -outputdir ${NEWworkingdir}/Exonerate -targetgenome ${NEWtargetgenome} -querychunktotal 720 -hintpriority 4 -hintsource E
# removing entries with multiple alignments due to identical best alignment score in the hints file
export TRINITY_BASE=$(basename ${trinity_out} | sed 's![^.]*$!!')
cd ${NEWworkingdir}/Exonerate
grep "gene" ${TRINITY_BASE}exonerate.target.genes.gff3 | cut -f9 | cut -f4 -d ";" | sed 's/from=//g' | sort | uniq -d | sort > cDNA_duplicated
grep -v -w -f cDNA_duplicated ${NEWworkingdir}/Exonerate/${TRINITY_BASE}exonerate.target.hints.gff3 > ${NEWworkingdir}/Exonerate/${TRINITY_BASE}exonerate.target.hints_unique.gff3

# Making the Augustus directories
cd ${NEWworkingdir}
mkdir Augustus
mkdir Augustus/annotation
mkdir Augustus/annotation_functional
mkdir Augustus/annotation_stitch
mkdir Augustus/annotation_stitch/blast
mkdir Augustus/annotation_stitch/blast/database
mkdir Augustus/annotation_stitch/intersect_and_stitch
mkdir Augustus/split_genome
mkdir Augustus/hints

# moving hints file to a more identifiable place
mv ${NEWworkingdir}/Exonerate/${TRINITY_BASE}exonerate.target.hints_unique.gff3 ${NEWworkingdir}/Augustus/hints/exonerate_hints.gff3

# splitting genomes into ${PBS_NCPUS} chunks (16 in this case) for parallel augustus execution
fastasplit -f ${NEWtargetgenome} -o ${NEWworkingdir}/Augustus/split_genome -c ${PBS_NCPUS}

# running augustus gene prediction
cd ${NEWworkingdir}/Augustus/split_genome
export GENOME_BASE=$(basename ${NEWtargetgenome})
for i in $(ls -1v ${GENOME_BASE}_chunk* | sed "s/${GENOME_BASE}_chunk_//g"); do echo augustus --species=${species} --softmasking=1 --UTR=on --print_utr=on --hintsfile=${NEWworkingdir}/Augustus/hints/exonerate_hints.gff3 --extrinsicCfgFile=${workingdir}/Augustus/config/extrinsic/extrinsic.MPE.cfg --exonnames=off --stopCodonExcludedFromCDS=t --AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config ${NEWworkingdir}/Augustus/split_genome/${GENOME_BASE}_chunk_"$i" ">" ${NEWworkingdir}/Augustus/annotation/aug"$i".out; done | parallel --jobs ${PBS_NCPUS} {}

for i in $(ls -1v ${GENOME_BASE}_chunk* | sed "s/${GENOME_BASE}_chunk_//g"); do cat ${NEWworkingdir}/Augustus/annotation/aug"$i".out >> ${NEWworkingdir}/Augustus/annotation/aug.out; done

cd ${NEWworkingdir}/Augustus/annotation
cat aug.out | join_aug_pred.pl > augustus.gff
getAnnoFasta.pl augustus.gff
grep -v "#" augustus.gff > augustus.gtf
gtf2gff.pl < augustus.gtf --out=augustus.gff3 --gff3
grep -w 'CDS\|three_prime_utr\|five_prime_utr' augustus.gff3 > augustus_for_bedintersect.gff3

# blast related species cDNA to new genome
cd ${NEWworkingdir}/Augustus/annotation_stitch/blast/database
makeblastdb -in ${NEWtargetgenome} \
-parse_seqids \
-title "genome" -dbtype nucl \
-out ${NEWworkingdir}/Augustus/annotation_stitch/blast/database/genome

blastn -query ${Related_species} \
-db ${NEWworkingdir}/Augustus/annotation_stitch/blast/database/genome \
-max_target_seqs 1 \
-evalue 1e-6 \
-parse_deflines \
-outfmt 6 \
-num_threads ${PBS_NCPUS} \
-out ${NEWworkingdir}/Augustus/annotation_stitch/blast/related_cDNA_to_genome.tsv

cd ${NEWworkingdir}/Augustus/annotation_stitch/blast
join <(sort related_cDNA_to_genome.tsv) <(grep ">" ${Related_species} | cut -f1,4 -d' ' | sed 's/ /\t/g' | sed 's/>//g' | sort) -t $'\t' -o1.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14 | perl -pe 's/\t/|/' | perl -lne '@a=split("\t", $_); print "$a[1]\t". (($a[8] < $a[9]) ? $a[8] : $a[9] )."\t".(($a[8] < $a[9]) ? $a[9] : $a[8])."\t$a[0]\t$a[11]\t".(($a[8] < $a[9]) ? "+" : "-" )' > related_cDNA_to_genome.bed
mv related_cDNA_to_genome.bed ${NEWworkingdir}/Augustus/annotation_stitch/intersect_and_stitch/related_cDNA_to_genome.bed

# generate list of genes to be stitched together, (stitching candidates)
cd ${NEWworkingdir}/Augustus/annotation_stitch/intersect_and_stitch
awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' ../../annotation/augustus_for_bedintersect.gff3 > augustus.bed

bedtools intersect -a related_cDNA_to_genome.bed -b augustus.bed -wa -wb -s > intersect
cut -f4,10 intersect | sed 's/;/\t/g' | cut -f1,3 | sed 's/Parent=//g' | sed 's/\.t1//g' | sort | uniq | sort | awk '{ printf "%s", (NR==1 || pre!=$1? (NR>1? ORS:"")$1: "") OFS $2; pre=$1 }
    END  { print "" }' | sed 's/ /\t/g' | awk 'NF>=3' > intersect2
cut -f2- -d "|" intersect2 | sort | uniq > intersect3.tsv

grep gene ../../annotation/augustus.gff3 | grep "+" | cut -f9 | sed 's/;//g' | sed 's/ID=//g' | paste -s -d "\t" > consecutive_positive.tsv
grep gene ../../annotation/augustus.gff3 | grep "-" | cut -f9 | sed 's/;//g' | sed 's/ID=//g' | paste -s -d "\t" > consecutive_negative.tsv

cut -f2- intersect3.tsv | while read a; do if grep -q "$a" consecutive_positive.tsv; then echo -e yes'\t'"$a"; else echo -e no'\t'"$a"; fi done > intersect4_pos.tsv
grep yes intersect4_pos.tsv | cut -f2- | sort | uniq > intersect5_pos.tsv
cat intersect5_pos.tsv | while read -r a; do OCCUR=$(grep "${a}" intersect5_pos.tsv | wc -l); echo -e $OCCUR'\t'"${a}"; done | awk '$1 == 1 {print ;}' | cut -f2- > intersect6_pos.tsv

cut -f2- intersect3.tsv | while read a; do if grep -q "$a" consecutive_negative.tsv; then echo -e yes'\t'"$a"; else echo -e no'\t'"$a"; fi done > intersect4_neg.tsv
grep yes intersect4_neg.tsv | cut -f2- | sort | uniq > intersect5_neg.tsv
cat intersect5_neg.tsv | while read -r a; do OCCUR=$(grep "${a}" intersect5_neg.tsv | wc -l); echo -e $OCCUR'\t'"${a}"; done | awk '$1 == 1 {print ;}' | cut -f2- > intersect6_neg.tsv

# stitching all the clusters on positive strand
StitchedGeneCount=1
cd ${NEWworkingdir}/Augustus/annotation_stitch
cat intersect_and_stitch/intersect6_pos.tsv | while read a; do last=$(echo $a | wc -w); one=1; lastone=$[last - one]; echo "$a" > temp_input; perl -i -pe "chomp if eof" temp_input; readarray -d $'\t' -t a < temp_input; for index in "${!a[@]}"; do if [ $index == 0 ]; then    \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'gene\|mRNA\|transcription_start_site\|five_prime_utr\|start_codon\|CDS\|intron' > temp); \
elif [ $index == $lastone ]; then \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'transcription_end_site\|stop_codon\|three_prime_utr\|CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
else \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
fi done; gt gff3 -sort -tidy -retainids -addintrons <(sed 's/CDS/exon/g' temp) > temp2; grep "$(cut -f3,4,5 temp2 | grep intron | sort | uniq -u)" temp2 >> temp; gt gff3 -sort -tidy -retainids temp | grep -v "#" > temp2; endcoord=$(sort -n -r -k 5 temp2 | cut -f5 | head -n 1); CurrentGeneID="StitchedPg"${StitchedGeneCount}; \
perl ${workingdir}/Scripts/AppendStitchedGene.pl ${NEWworkingdir}/Augustus/annotation_stitch/temp ${CurrentGeneID} | awk '$2=="."{$2="Inferred"} 1' OFS="\t" | awk '$3=="gene"{$5='$endcoord'} 1' OFS="\t" | awk '$3=="mRNA"{$5='$endcoord'} 1' OFS="\t" > temp3; \
output=$(awk '$2 == "Inferred"' temp3 | awk -v maxintron="$maxintron" '$5-$4 <= maxintron || maxintron == 0 {print}'); if [ -n "$output" ]; \
then cat temp3 >> temp_all; echo -e "$(cat temp_input)" | sed 's/\t/,/g' | awk -v count="${CurrentGeneID}" '{print count"\t"$1}' >> sourcegene; echo -e "$(cat temp_input)" >> sourcegene2; ((StitchedGeneCount++)); \
else grep -w -f <(tr '\t' '\n' < temp_input) ../annotation/augustus.gff3 >> temp_all; \
fi; \
rm temp; \
rm temp2; \
rm temp3; done

# stitching all the clusters on negative strand
StitchedGeneCount=1
cd ${NEWworkingdir}/Augustus/annotation_stitch
cat intersect_and_stitch/intersect6_neg.tsv | while read a; do last=$(echo $a | wc -w); one=1; lastone=$[last - one]; echo "$a" > temp_input; perl -i -pe "chomp if eof" temp_input; readarray -d $'\t' -t a < temp_input; for index in "${!a[@]}"; do if [ $index == 0 ]; then    \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'gene\|mRNA\|transcription_end_site\|three_prime_utr\|stop_codon\|CDS\|intron' > temp); \
elif [ $index == $lastone ]; then \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'transcription_start_site\|start_codon\|five_prime_utr\|CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
else \
$(grep -w "${a[$index]}" ../annotation/augustus.gff3 | grep 'CDS\|intron' | awk -v pattern="${a[$index]}" -v replacement="${a[0]}" 'BEGIN{FS="\t"; OFS="\t"} {gsub(pattern, replacement, $9)} 1' >> temp); \
fi done; gt gff3 -sort -tidy -retainids -addintrons <(sed 's/CDS/exon/g' temp) > temp2; grep "$(cut -f3,4,5 temp2 | grep intron | sort | uniq -u)" temp2 >> temp; gt gff3 -sort -tidy -retainids temp | grep -v "#" > temp2; endcoord=$(sort -n -r -k 5 temp2 | cut -f5 | head -n 1); CurrentGeneID="StitchedNg"${StitchedGeneCount}; \
perl ${workingdir}/Scripts/AppendStitchedGene.pl ${NEWworkingdir}/Augustus/annotation_stitch/temp ${CurrentGeneID} | awk '$2=="."{$2="Inferred"} 1' OFS="\t" | awk '$3=="gene"{$5='$endcoord'} 1' OFS="\t" | awk '$3=="mRNA"{$5='$endcoord'} 1' OFS="\t" > temp3; \
output=$(awk '$2 == "Inferred"' temp3 | awk -v maxintron="$maxintron" '$5-$4 <= maxintron || maxintron == 0 {print}'); if [ -n "$output" ]; \
then cat temp3 >> temp_all; echo -e "$(cat temp_input)" | sed 's/\t/,/g' | awk -v count="${CurrentGeneID}" '{print count"\t"$1}' >> sourcegene; echo -e "$(cat temp_input)" >> sourcegene2; ((StitchedGeneCount++)); \
else grep -w -f <(tr '\t' '\n' < temp_input) ../annotation/augustus.gff3 >> temp_all; \
fi; \
rm temp; \
rm temp2; \
rm temp3; done

rm temp_input

tr '\t' '\n' < sourcegene2 | sort | uniq > sourcegene3

tr '\t' '\n' < intersect_and_stitch/intersect6_pos.tsv > single.lst
tr '\t' '\n' < intersect_and_stitch/intersect6_neg.tsv >> single.lst
sort single.lst | uniq > single2.lst
grep -v -w -f single2.lst ../annotation/augustus.gff3 > augustus2.gff3

gt gff3 -retainids -tidy <(cat augustus2.gff3 temp_all) > augustus_stitched.gff3

gffread -y temporary.fa -g ${NEWtargetgenome} augustus_stitched.gff3
seqtk seq -l 60 temporary.fa | seqkit sort -N > augustus_stitched_protein.fa
sed -i '/^[^>]/s/\./X/g' augustus_stitched_protein.fa
rm temporary.fa

cp ${NEWworkingdir}/Augustus/annotation_stitch/augustus_stitched_protein.fa ${NEWworkingdir}/Augustus_stitched_proteins.fasta

# Blasting augustus proteins to uniprot, printing the full fasta header in uniprot DB
cd ${NEWworkingdir}/Augustus/annotation_functional
blastp -query ${NEWworkingdir}/Augustus/annotation_stitch/augustus_stitched_protein.fa \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-evalue 1e-3 \
-outfmt "6 std salltitles" \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads ${PBS_NCPUS} \
-out ${NEWworkingdir}/Augustus/annotation_functional/blast.tsv

seqtk subseq -l 60 ${NEWworkingdir}/Augustus/annotation/augustus.aa <(cat ${NEWworkingdir}/Augustus/annotation_stitch/sourcegene3 | awk '{print $0".t1"}') > sourceProt.aa
blastp -query ${NEWworkingdir}/Augustus/annotation_functional/sourceProt.aa \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-evalue 1e-3 \
-outfmt "6 std salltitles" \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads ${PBS_NCPUS} \
-out ${NEWworkingdir}/Augustus/annotation_functional/blast2.tsv

join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(grep -w "gene" ${NEWworkingdir}/Augustus/annotation_stitch/augustus_stitched.gff3 | cut -f9 | sed 's/ID=//g' | awk '{print $0".t1"}' | sort -k1,1) <(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); print $1"\t"$2"\t" b[1]} else {print $1"\t"$2"\tunknown"}}' ${NEWworkingdir}/Augustus/annotation_functional/blast.tsv | sort -k1,1) | sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' > geneID_to_uniprotID.tabular
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9";uniprot_ID="a[$9]}1' geneID_to_uniprotID.tabular ${NEWworkingdir}/Augustus/annotation_stitch/augustus_stitched.gff3 > ${NEWworkingdir}/Augustus_stitched_annotation.gff3

join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(grep -w -f ${NEWworkingdir}/Augustus/annotation_stitch/sourcegene3 ${NEWworkingdir}/Augustus/annotation/augustus.gff3 | grep -w "gene" | cut -f9 | sed 's/ID=//g' | sed 's/;//g' | awk '{print $0".t1"}' | sort -k1,1) <(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); print $1"\t"$2"\t" b[1]} else {print $1"\t"$2"\tunknown"}}' ${NEWworkingdir}/Augustus/annotation_functional/blast2.tsv | sort -k1,1) | sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' > sourceID_to_uniprotID.tabular
gt gff3 -retainids -tidy <(grep -w -f ${NEWworkingdir}/Augustus/annotation_stitch/sourcegene3 ${NEWworkingdir}/Augustus/annotation/augustus.gff3) > sourceGene.gff3
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9";uniprot_ID="a[$9]}1' sourceID_to_uniprotID.tabular ${NEWworkingdir}/Augustus/annotation_functional/sourceGene.gff3 > ${NEWworkingdir}/Augustus_source_for_stitching.lst

cd ${NEWworkingdir}
echo -e "Contig\tMiddle position\tGene length\tStrand\tGene ID\tUniprot accession number\tGene name" > ${NEWworkingdir}/Augustus_gene_table.tabular
awk -F'\t' '$3 == "gene" {midpoint = int(($4 + $5) / 2); split($9,a,"ID="); split(a[2],b,";"); split($9,c,"uniprot_ID="); split(c[2],d,";"); split($9,e,"gene_name="); split(e[2],f," "); print $1"\t"midpoint"\t"$5-$4"\t"$7"\t"b[1]"\t"d[1]"\t"f[1]}' ${NEWworkingdir}/Augustus_stitched_annotation.gff3 >> ${NEWworkingdir}/Augustus_gene_table.tabular

cat ${NEWworkingdir}/Augustus/annotation_stitch/sourcegene | cat - ${NEWworkingdir}/Augustus_source_for_stitching.lst > temp && mv temp ${NEWworkingdir}/Augustus_source_for_stitching.lst
echo -e "Stitched gene ID\tSource genes (comma separated)" | cat - ${NEWworkingdir}/Augustus_source_for_stitching.lst > temp && mv temp ${NEWworkingdir}/Augustus_source_for_stitching.lst
cp ${NEWworkingdir}/Augustus/annotation_functional/sourceProt.aa ${NEWworkingdir}/Augustus_source_for_stitching.fasta
