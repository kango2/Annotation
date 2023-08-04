# set up your variables and paths
## Required
export NEWworkingdir=""
export NEWtargetgenome=""
export workingdir=""
export trinity_out=""
export species=""
export project_code=""
## do you want email notification as the pipeline runs? Good for keeping track of progression
## CAPITAL YES or NO
## if YES, set up $email
export emailOPT="NO"
export email=""
## do you want to run Option A or B? A for stitching and B for no stitching
## CAPITAL A or B
## if A, set up $Related_species and $maxintron
export Option="B"
export Related_species=""
export maxintron=""


# =================================================================
# below are script to start pipeline, no need to change
mkdir -p ${NEWworkingdir}
cd ${NEWworkingdir}
wget https://github.com/kango2/Annotation/raw/main/OptionA_Stitch.sh
wget https://github.com/kango2/Annotation/raw/main/OptionB_NoStitch.sh

if [ "$emailOPT" = "YES" ]; then
    cd ${NEWworkingdir}
    echo -e "====================\nemailOPT set to: YES"
    echo -e "Adding your email to scripts..."
    sed -i "5i\\#PBS -M ${email}" Option*
    sed -i '6i\\#PBS -m ae' Option*
    echo -e "Done"
elif [ "$emailOPT" = "NO" ]; then
    echo -e "====================\nemailOPT set to: NO"
else
    echo -e "Unexpected value for emailOPT, please select either YES or NO. Stopping script."
    exit 1
fi


cd ${NEWworkingdir}

if [ "Option" = "A" ]; then
    echo -e "Option A selected, now submitting jobs to start pipeline..."
    DEPEND_FOR_EXONERATE=$(qsub -W depend=on:3 -P ${project_code} -v workingdir=${workingdir},NEWworkingdir=${NEWworkingdir},NEWtargetgenome=${NEWtargetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code},Related_species=${Related_species},maxintron=${maxintron} OptionA_Stitch.sh)
elif [ "Option" = "B" ]; then
    echo -e "Option B selected, now submitting jobs to start pipeline..."
    DEPEND_FOR_EXONERATE=$(qsub -W depend=on:3 -P ${project_code} -v workingdir=${workingdir},NEWworkingdir=${NEWworkingdir},NEWtargetgenome=${NEWtargetgenome},trinity_out=${trinity_out},species=${species},project_code=${project_code} OptionB_NoStitch.s)
else
    echo -e "Unexpected value for Option, please select either A or B. Stopping script."
    exit 1
fi

mkdir ${NEWworkingdir}/Exonerate
for i in ${trinity_out}; do inputfasta=$(realpath $i); for c in 1 241 481; do qsub -W depend=beforeok:${DEPEND_FOR_EXONERATE} -P ${project_code} -o ${NEWworkingdir} -v querychunktotal=720,querychunkstart=$c,querychunkend=$((c+239)),outputdir=${NEWworkingdir}/Exonerate,inputfasta=${inputfasta},targetgenome=${NEWtargetgenome} ${workingdir}/Scripts/runexonerate.sh; done; done

echo -e "Done, pipeline successfully started"
