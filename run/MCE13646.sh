#$ -cwd -V -m e
#$ -pe smp 2
#$ -l h_rt=40:00:00
#$ -l h_vmem=4G
#$ -t 1-50
cd /nobackup/phy4cs/MCEv2-SB-13646/$SGE_TASK_ID-run/
echo "Running on $HOSTNAME in folder $PWD"
module load mkl
./MCE.exe
