# Genome annotation pipeline walkthrough
This pipeline trains and runs Augustus gene prediction for your new eukaryotic species using TRINITY de novo transcriptomic assembly.\
The codes are optimised to run on **NCI Gadi**, it will need adjusting if running on other servers.\
The entire pipeline is split into 5 parts, each as a PBS script.\
**NOTE** that all paths needs to be **absolute path** when setting up the shell variables below.
### Required inputs:
- Trinity assembly assembled from RNA-seq of your species (fasta file)
- Soft-masked genome assembly of your species (fasta file)
- Uniprot protein database (download instruction provided below, needs to be in `${workingdir}`)
### Outputs:
- Gene annotation (gff3)
- Predicted peptides (fasta)
- Optimised Augustus species-specific parameters for your species
## Running the pipeline (AUTO)
Added a new script to start the pipeline automatically, and submits all five jobs in the pipeline at once with starting dependencies. (Basically equivalent to running the entire next section but easier and simpler, and does all the job submitting automatically for you).
- Download the script and make it executable
```
wget https://github.com/kango2/Annotation/raw/main/PipelineSettings_StartPipeline.sh
chmod u+x PipelineSettings_StartPipeline.sh
```
- Edit first section in any text editor to set up your variables (all lines starting with "export")
- Run the script
```
./PipelineSettings_StartPipeline.sh
```
## Running the pipeline
**Step 1.**  Create a working directory
```
export workingdir=/path/to/working_directory
mkdir ${workingdir}
```
**Step 2.** Download the latest Uniprot database fasta file and the directory structure into your working directory
```
cd ${workingdir}
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

wget https://github.com/kango2/Annotation/raw/main/Augustus_workingdir.tar.gz
tar -xzvf Augustus_workingdir.tar.gz --strip-components=1
```
**Step 3.** Setting up required shell variables
- use underscores instead of space characters for the name of your species.
- use `printenv PROJECT` if you are unsure of what your project code is, this will be needed for submitting PBS jobs.
```
export targetgenome=/path/to/your/softmasked_genome.fa
export trinity_out=/path/to/your/trinity.fa
export species=Species_name
export project_code=xl04
```
**Step 4 (Optional).** Add email notification for your PBS jobs
- Replace "your email" in second command with your actual email.
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
- A log file will be generated in `${workingdir}/annotation.log` this file reports some basic information as the pipeline progresses, it will also automatically generate the `qsub` commands for you to execute and start the next script.
- As a general rule, lines beginning with `[!]` in the log file will require attention, these are the instruction and commands needed to run the next part.
## Part 5 of pipeline
There are two options to choose from when running Part 5 (5A or 5B). Instructions for both options will be generated in the log file, **choose and run only one**.\
\
**Option A** will stitch together neighbouring genes (no intervening genes) based on their overlapping with cDNA alignments to the genome (cDNA from a closely related species), then do a similarity blast search to uniprot for gene identification
- Two additional shell variables will be needed to run option A, instructions on setting them up will be in the log file. They are `${Related_species}` and `${maxintron}`
- If the stitched gene contains an inferred intron with length longer than `${maxintron}` then it will not stitch them together, it will instead use the original source genes in the final annotation. Setting it to `0` will stitch together all candidates regardless of the inferred intron lengths. I recommend setting this to a reasonable length, maybe `200000`.
- Note that currently this stitching script only works with cDNA fasta file from ENSEMBL.
<!-- end of the list -->
**Option B** will just do a similarity blast search to uniprot for gene identification

# Annotating new genome for same species using existing training parameters (AUTO)
Similar to the **AUTO** section above but for running with existing training parameters.
- Download the script and make it executable
```
wget https://github.com/kango2/Annotation/raw/main/PipelineSettings_StartPipeline2.sh
chmod u+x PipelineSettings_StartPipeline2.sh
```
- Edit first section in any text editor to set up your variables (all lines starting with "export")
- Run the script
```
./PipelineSettings_StartPipeline2.sh
```
# Annotating new genome for same species using existing training parameters
**Step 1.**  Create a NEW working directory
```
export NEWworkingdir=/path/to/new_working_directory
mkdir ${NEWworkingdir}
```
**Step 2.** Download script
```
cd ${NEWworkingdir}
wget https://github.com/kango2/Annotation/raw/main/OptionA_Stitch.sh
wget https://github.com/kango2/Annotation/raw/main/OptionB_NoStitch.sh
```
**Step 3.** Setting up required shell variables, variable 2-5 should be exactly same as running it the first time above
- use underscores instead of space characters for the name of your species.
- use `printenv PROJECT` if you are unsure of what your project code is, this will be needed for submitting PBS jobs.
```
export NEWtargetgenome=/path/to/your/new_softmasked_genome.fa
export workingdir=/path/to/working_directory
export trinity_out=/path/to/your/trinity.fa
export species=Species_name
export project_code=xl04
```
- If running OptionA_Stitch.sh (Similar to Part5A), set up these two additional variables, see **Part 5 of pipeline** above for explanation.
```
export Related_species=/path/to/related_species_cDNA.fa
export maxintron=0
```
**Step 4.** Generate Exonerate hints
```
mkdir ${NEWworkingdir}/Exonerate
for i in ${trinity_out}; do inputfasta=$(realpath $i); for c in 1 241 481; do qsub -P ${project_code} -o ${NEWworkingdir} -v querychunktotal=720,querychunkstart=$c,querychunkend=$((c+239)),outputdir=${NEWworkingdir}/Exonerate,inputfasta=${inputfasta},targetgenome=${NEWtargetgenome} ${workingdir}/Scripts/runexonerate.sh; done; done
```
**Step 5.** Run OptionA_Stitch.sh or OptionB_NoStitch.sh

Make sure to wait until the 3 exonerate jobs finish before running this step
- Option A
```
cd ${NEWworkingdir}
qsub -P ${project_code} -v workingdir=${workingdir},NEWworkingdir=${NEWworkingdir},NEWtargetgenome=${NEWtargetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code},Related_species=${Related_species},maxintron=${maxintron} OptionA_Stitch.sh
```
- Option B
```
cd ${NEWworkingdir}
qsub -P ${project_code} -v workingdir=${workingdir},NEWworkingdir=${NEWworkingdir},NEWtargetgenome=${NEWtargetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} OptionB_NoStitch.sh
```
