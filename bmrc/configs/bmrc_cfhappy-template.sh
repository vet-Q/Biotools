#/bin/bash

#$ -N cfhappy
#$ -P hendry.prjc -q short.qc
#$ -o logs/stdout_cfhappy -e logs/stderr_cfhappy -j y
#$ -cwd -V

echo "**********************************************************************"
echo "BMRC CLUSTER JOB SUBMISSION OF HAP.PY COMPARISON"
echo "----------------------------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
SECONDS=0
echo "**********************************************************************"

args=$@
truthset cfhappy $args

echo "**********************************************************************"
echo "Finished at: "`date`
printf 'Time elapsed: %dh:%dm:%ds\n' $(($SECONDS/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60))
echo "**********************************************************************"
echo ""
echo ""
echo ""
