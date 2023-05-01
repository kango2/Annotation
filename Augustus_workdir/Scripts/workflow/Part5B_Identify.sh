#!/bin/bash
#PBS -N Part5B_Identify
#PBS -l ncpus=32,walltime=12:00:00,storage=gdata/if89+gdata/xl04,mem=80GB,jobfs=80GB
#PBS -j oe

###################################################################################
### Part five B of Genome Annotation pipeline script

### Setting up the log file
export log=${workingdir}/annotation.log
echo -e "\n\n\n" >> ${log}
date | tee >> ${log}
echo -e "Part FiveB of annotation pipeline: Similarity search to Uniprot database (BLAST)\n====================\n Copying protein fasta to working directory." >> ${log}

### Loading modules
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 blast/2.11.0 genometools/1.6.2 seqkit/2.3.1

###################################################################################
### Structural identification of the predicted genes by similarity search to uniprot proteins
echo -e "\n Blasting Augustus proteins to uniprot database." >> ${log}

# copying protein fasta file to workingdir
cp ${workingdir}/Augustus/annotation/augustus.aa ${workingdir}/Augustus_proteins.fasta
cd ${workingdir}/Augustus/annotation_functional

# Blasting augustus proteins to uniprot, printing the full fasta header in uniprot DB
blastp -query ${workingdir}/Augustus/annotation/augustus.aa \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-evalue 1e-3 \
-outfmt "6 std salltitles" \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads ${PBS_NCPUS} \
-out ${workingdir}/Augustus/annotation_functional/blast.tsv

### concatenating the blast output to the final annotation as two extra attributes on 9th column of "gene" features, separated by ";"

# if a protein does not have gene name in the header then it prints unknown instead
echo -e " Concatenating blast results into annotation gff3 file and copying to working directory." >> ${log}

join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(grep -w "gene" ${workingdir}/Augustus/annotation/augustus.gff3 | cut -f9 | sed 's/ID=//g' | sed 's/;//g' | awk '{print $0".t1"}' | sort -k1,1) <(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); print $1"\t"$2"\t" b[1]} else {print $1"\t"$2"\tunknown"}}' ${workingdir}/Augustus/annotation_functional/blast.tsv | sort -k1,1) | sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' | awk '{print $1";\t"$2"\t"$3}' > geneID_to_uniprotID.tabular
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9"uniprot_ID="a[$9]";"}1' geneID_to_uniprotID.tabular ${workingdir}/Augustus/annotation/augustus.gff3 > ${workingdir}/Augustus_annotation.gff3

### generating a tabular file with gene info
echo -e "\n Generating a tabular file with some basic gene information." >> ${log}

cd ${workingdir}
echo -e "Contig\tMiddle position\tGene length\tStrand\tGene ID\tUniprot accession number\tGene name" > ${workingdir}/Augustus_gene_table.tabular

# b[1] is geneid, d[1] is uniprot AC, f[1] is gene name
awk -F'\t' '$3 == "gene" {midpoint = int(($4 + $5) / 2); split($9,a,"ID="); split(a[2],b,";"); split($9,c,"uniprot_ID="); split(c[2],d,";"); split($9,e,"gene_name="); split(e[2],f,";"); print $1"\t"midpoint"\t"$5-$4"\t"$7"\t"b[1]"\t"d[1]"\t"f[1]}' ${workingdir}/Augustus_annotation.gff3 >> ${workingdir}/Augustus_gene_table.tabular

###################################################################################
### Annotation complete
time=$(date)
echo -e "\n #################### ${time}\n Part FiveB completed.\n\n Annotation pipeline completed, you can find your output files in:\n\t${workingdir}" >> ${log}
