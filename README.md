# Genome annotation pipeline walkthrough
This pipeline trains and runs Augustus gene prediction for your new eukaryotic species.\
The codes are intended to be running on NCI Gadi, it will need adjusting if running on other servers.\
The entire pipeline is split into 5 parts, each as a PBS script.\
NOTE that all paths needs to be absolute path when setting up the shell variables below.
### Required inputs:
- Trinity assembly assembled from RNA-seq of your species (fasta file)
- Soft-masked genome assembly of your species (fasta file)
### Outputs:
- Gene annotation (gff3)
- Predicted peptides (fasta)
- Optimised Augustus species-specific parameters for your species
## Running the pipeline
**Step 1.**  Create a working directory
```
export workingdir=/path/to/working_directory
mkdir ${workingdir}
```
**Step 2.** Download the latest Uniprot database fasta file, and the directory sturcture in your working directory
```
cd ${workingdir}
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

wget https://github.com/kango2/Annotation/blob/main/Augustus_workingdir.tar.gz
tar -xzvf Augustus_workingdir.tar.gz --strip-components=1
```
**Step 3.** Setting up required shell variables
- use underscores instead of space characters for the name of your species.
- use `echo ${PROJECT}` if you are unsure of what your project code is, this will be needed for submitting PBS jobs.
```
export targetgenome=/path/to/your/softmasked_genome.fa
export trinityout=/path/to/your/trinity.fa
export species=Species_name
export project_code=xl04
```
**Step 4 (Optional).** Add email notification for your PBS jobs
```
cd ${workingdir}/Scripts/workflow
export email="your email"
sed -i "5i\\#PBS -M ${email}" Part*
sed -i '6i\\#PBS -m ae' Part*
```
**Step 5.** Run the first script
```
cd ${workingdir}/Scripts/workflow
qsub -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part1_trainingset.sh
```
**Step 6.** Check log for further instruction
- A log file will be generated in `${workingdir}annotation.log` this file reports some basic information as the pipeline progresses, it will also automatically generate the `qsub` commands for you to execute and start the next script.
- As a general rule, lines beginning with `[!]` in the log file will require attention, these are the instruction and commands needed run the next part.
