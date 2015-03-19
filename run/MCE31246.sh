#$ -cwd -V -m e
#$ -l h_vmem=2G
#$ -t 1-1
cd /home/ds/phy4cs/Dropbox/PhysChem/MCE/EXEC/MCEv2-SB-31246/$SGE_TASK_ID-run/
echo "Running on $HOSTNAME in folder $PWD"
./MCE.exe
