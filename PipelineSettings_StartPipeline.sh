# set up your variables and paths
## Required
export workingdir=""
export targetgenome=""
export trinity_out=""
export species=""
export project_code=""
## do you want email notification as the pipeline runs? Good for keeping track of progression
## CAPITAL YES or NO
## if YES, set up $email
export emailOPT="NO"
export email=""
## do you want to run Part5A or Part5B? A for stitching and B for no stitching
## CAPITAL A or B
## if A, set up $Related_species and $maxintron
export Part5="B"
export Related_species=""
export maxintron=""


# =================================================================
# below are script to start pipeline, no need to change
mkdir -p ${workingdir}
cd ${workingdir}
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget https://github.com/kango2/Annotation/raw/main/Augustus_workingdirExpress.tar.gz
tar -xzvf Augustus_workingdirExpress.tar.gz --strip-components=1

if [ "$emailOPT" = "YES" ]; then
    cd ${workingdir}/Scripts/workflow
    echo -e "====================\nemailOPT set to: YES"
    echo -e "Adding your email to scripts..."
    sed -i "5i\\#PBS -M ${email}" Part*
    sed -i '6i\\#PBS -m ae' Part*
    echo -e "Done"
elif [ "$emailOPT" = "NO" ]; then
    echo -e "====================\nemailOPT set to: NO"
else
    echo -e "Unexpected value for emailOPT, please select either YES or NO. Stopping script."
    exit 1
fi


cd ${workingdir}/Scripts/workflow

if [ "Part5" = "A" ]; then
    echo -e "Option A selected for Part 5, now submitting jobs to start pipeline..."
    DEPEND_FOR_4=$(qsub -W depend=on:1 -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code},Related_species=${Related_species},maxintron=${maxintron} Part5A_StitchThenIdentify.sh)
elif [ "Part5" = "B" ]; then
    echo -e "Option B selected for Part 5, now submitting jobs to start pipeline..."
    DEPEND_FOR_4=$(qsub -W depend=on:1 -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part5B_Identify.sh)
else
    echo -e "Unexpected value for Part5, please select either A or B. Stopping script."
    exit 1
fi
DEPEND_FOR_3=$(qsub -W depend=on:1,beforeok=${DEPEND_FOR_4} -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part4_RunningAugustus.sh)
DEPEND_FOR_2=$(qsub -W depend=on:1,beforeok=${DEPEND_FOR_3} -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part3_TrainingAugustus.sh)
DEPEND_FOR_1=$(qsub -W depend=on:7,beforeok=${DEPEND_for_2} -P ${project_code} -v workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part2_trainingset.sh)
qsub -W depend=beforeok:${DEPEND_FOR_1} -P ${project_code} -v DEPEND_FOR_1=${DEPEND_FOR_1},workingdir=${workingdir},targetgenome=${targetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} Part1_trainingset.sh

echo -e "Done, pipeline successfully started"
