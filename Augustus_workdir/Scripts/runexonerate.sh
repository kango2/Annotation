#!/bin/bash
#PBS -N exonerate
#PBS -l ncpus=48,walltime=10:00:00,storage=gdata/if89+gdata/xl04,mem=190GB,jobfs=400GB
#PBS -j oe
module use /g/data/if89/apps/modulefiles
module load exonerate/2.2.0 parallel

#output directory?
#cd /srv/scratch/waters/LIS9486_LIS9514/Exonerate_output/testingfolder
set -ex
rsync -a ${targetgenome} ${PBS_JOBFS}/
rsync -a ${inputfasta} ${PBS_JOBFS}/

inputfasta=${PBS_JOBFS}/$(basename ${inputfasta})
targetgenome=${PBS_JOBFS}/$(basename ${targetgenome})

seq ${querychunkstart} ${querychunkend} | parallel --resume-failed --results ${outputdir}/$(basename ${inputfasta} .fasta)_{}.exonerate.out \
--joblog ${outputdir}/$(basename ${inputfasta} .fasta)_vs_$(basename ${targetgenome}).C${querychunkstart}-C${querychunkend}.exonerate.parallel.log --jobs ${PBS_NCPUS} \
exonerate \
--querytype dna \
--targettype dna \
--model est2genome \
--hspfilter 100 \
--dpmemory 256 \
--seedrepeat 3 \
--wordjump 5 \
--fsmmemory 3000 \
--softmasktarget TRUE \
--saturatethreshold 50 \
--dnawordlen 17 \
--bestn 1 \
--percent 90 \
--showtargetgff yes \
--showalignment no \
--showcigar yes \
--showquerygff yes \
--showvulgar no \
--querychunkid {} \
--querychunktotal ${querychunktotal} \
--query ${inputfasta} \
--target ${targetgenome}