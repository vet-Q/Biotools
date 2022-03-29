#/bin/bash

#$ -N {name}
#$ -P hendry.prjc -q {queue}
#$ -o logs/stdout_{log} -e logs/stderr_{log} -j y
#$ -cwd -V
{extra_qsub_args}

echo "**********************************************************************"
echo "ARRAY JOB SUBMISSION OF {script}"
echo "----------------------------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
SECONDS=0
echo "**********************************************************************"

{script} {script_args} {extra_script_args}

echo "**********************************************************************"
echo "Finished at: "`date`
printf 'Time elapsed: %dh:%dm:%ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60))
echo "**********************************************************************"
echo ""
echo ""
echo ""
