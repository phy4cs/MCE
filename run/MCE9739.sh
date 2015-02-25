#$ -cwd -V -m e
#$ -l h_rt=40:00:00
#$ -l h_vmem=4G
#$ -t 1-1
cd /home/ds/phy4cs/MCE/EXEC/MCEv2-SB-9739/$SGE_TASK_ID-run/
echo "Running on $HOSTNAME in folder $PWD"
./MCE.exe
