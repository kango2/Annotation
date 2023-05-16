#!/bin/bash
#PBS -N OptionB_NoStitch
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


cp ${NEWworkingdir}/Augustus/annotation/augustus.aa ${NEWworkingdir}/Augustus_proteins.fasta
cd ${NEWworkingdir}/Augustus/annotation_functional

blastp -query ${NEWworkingdir}/Augustus/annotation/augustus.aa \
-db ${workingdir}/Trinity_filter/first_filter/database/uniprot_sprot \
-parse_deflines \
-evalue 1e-3 \
-outfmt "6 std salltitles" \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads ${PBS_NCPUS} \
-out ${NEWworkingdir}/Augustus/annotation_functional/blast.tsv

join -1 1 -2 1 -t $'\t' -a 1 -e unknown -o1.1,2.2,2.3 <(grep -w "gene" ${NEWworkingdir}/Augustus/annotation/augustus.gff3 | cut -f9 | sed 's/ID=//g' | sed 's/;//g' | awk '{print $0".t1"}' | sort -k1,1) <(awk -F'\t' '{if ($13 ~ /GN=/) {split($13,a,"GN="); split(a[2],b," "); print $1"\t"$2"\t" b[1]} else {print $1"\t"$2"\tunknown"}}' ${NEWworkingdir}/Augustus/annotation_functional/blast.tsv | sort -k1,1) | sort -k2 -V | awk -F'\t' -v OFS='\t' '{sub(/\.t1$/,"",$1); print "ID="$0}' | awk '{print $1";\t"$2"\t"$3}' > geneID_to_uniprotID.tabular
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2";gene_name="$3; next} $9 in a {$9=$9"uniprot_ID="a[$9]";"}1' geneID_to_uniprotID.tabular ${NEWworkingdir}/Augustus/annotation/augustus.gff3 > ${NEWworkingdir}/Augustus_annotation.gff3

cd ${NEWworkingdir}
echo -e "Contig\tMiddle position\tGene length\tStrand\tGene ID\tUniprot accession number\tGene name" > ${NEWworkingdir}/Augustus_gene_table.tabular

awk -F'\t' '$3 == "gene" {midpoint = int(($4 + $5) / 2); split($9,a,"ID="); split(a[2],b,";"); split($9,c,"uniprot_ID="); split(c[2],d,";"); split($9,e,"gene_name="); split(e[2],f,";"); print $1"\t"midpoint"\t"$5-$4"\t"$7"\t"b[1]"\t"d[1]"\t"f[1]}' ${NEWworkingdir}/Augustus_annotation.gff3 >> ${NEWworkingdir}/Augustus_gene_table.tabular
