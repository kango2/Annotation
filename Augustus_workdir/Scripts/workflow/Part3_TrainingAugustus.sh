#!/bin/bash
#PBS -N Part3_TrainingAugustus
#PBS -l ncpus=16,walltime=30:00:00,storage=gdata/if89+gdata/xl04,mem=50GB,jobfs=50GB
#PBS -j oe

###################################################################################
### Part three of Genome Annotation pipeline script

### Setting up the log file
export log=${workingdir}/annotation.log
echo -e "\n\n\n" >> ${log}
date | tee >> ${log}
echo -e "Part Three of annotation pipeline: Training Augustus for your species\n====================\n Initial training with genes not chosen for optimising Augustus." >> ${log}

### Loading modules
module use /g/data/if89/apps/modulefiles
module load Augustus/3.4.0 perllib/v5.26.3 blat/37 RepeatMasker/4.1.2-p1 scipio/1.4 pblat/2.5 pslCDnaFilter/0 parallel/20191022 exonerate/2.2.0

### Copying augustus config files to a writable space and setting up path
cd ${workingdir}
rsync -a $AUGUSTUS_CONFIG_PATH/ Augustus/config/
export AUGUSTUS_CONFIG_PATH=${workingdir}/Augustus/config

###################################################################################
### Optimising exon parameters for your species

# First we make parameters files for the new species
new_species.pl --species=${species}

# an initial training using genes not selected for optimising, this is fast
cd ${workingdir}/Augustus/training
etraining --species=${species} ${workingdir}/Exonerate/annotation/CDS_only/training.gb.train
augustus --species=${species} ${workingdir}/Exonerate/annotation/CDS_only/training.gb.test | tee first_evaluation.out
grep -A 22 Evaluation first_evaluation.out > first_evaluation.report

# Now we optimize with 500 genes, 200 for evaluation and all 500 for training, the max 5 rounds has been chosen but it will finish earlier if no improvement are found
echo -e "\n Optimising species-specific parameters with 500 genes, you can check progress in the optimising log file.\n [LOG FILE]\t${workingdir}/Augustus/training/optimize.out" >> ${log}

optimize_augustus.pl \
--species=${species} \
--cpus=${PBS_NCPUS} \
--rounds=5 \
${workingdir}/Exonerate/annotation/CDS_only/training.gb.test.test \
--onlytrain=${workingdir}/Exonerate/annotation/CDS_only/training.gb.test.train \
> optimize.out

# Retrain after optimization
echo -e " Retraining after optimisation." >> ${log}

etraining --species=${species} ${workingdir}/Exonerate/annotation/CDS_only/training.gb.test

###################################################################################
### Optimising UTR parameters for your species

# we now optimize the utr parameters
echo -e "\n Now optimising species-specific UTR parameters, you can check progress in the optimising log file.\n [LOG FILE]\t${workingdir}/Augustus/utr_training/optimize.out" >> ${log}

cd ${workingdir}/Augustus/utr_training
optimize_augustus.pl \
--species=${species} \
--cpus=${PBS_NCPUS} \
--rounds=5 \
${workingdir}/Exonerate/annotation/CDS_with_UTR/training_with_utr.gb.test.test \
--onlytrain=${workingdir}/Exonerate/annotation/CDS_with_UTR/training_with_utr.gb.test.train \
--UTR=on \
--metapars=${workingdir}/Augustus/config/species/${species}/${species}_metapars.utr.cfg \
--trainOnlyUtr=1 > optimize_utr.out

# retrain after optimization
echo -e " Retraining after optimisation." >> ${log}

etraining --species=${species} --UTR=on ${workingdir}/Exonerate/annotation/CDS_with_UTR/training_with_utr.gb.test

# final evaluation
augustus --species=${species} --UTR=on ${workingdir}/Exonerate/annotation/CDS_with_UTR/training_with_utr.gb.test | tee UTR_evaluation.out
grep -A 22 Evaluation UTR_evaluation.out > UTR_evaluation.report

###################################################################################
### splitting genome into chunks in preparation for parallel augustus execution
echo -e "\n Splitting genome into chunks in preparation for running Augustus in parallel." >> ${log}

fastasplit -f ${workingdir}/Simple_header/genome/Target_genome_clean.fa -o ${workingdir}/Augustus/split_genome -c ${PBS_NCPUS}

###################################################################################
### Setting up commands to run for part 4 in the log file
time=$(date)
echo -e "\n #################### ${time}\n Part Three completed, run part Four with the following commands.\n CMDs to run:\n [!]\tcd ${workingdir}/Scripts/workflow\n [!]\tqsub -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part4_RunningAugustus.sh" >> ${log}
