#!/bin/bash
#SBATCH --partition generalsky
##SBATCH --partition=EpycPriority
##SBATCH --account=zhz18039
#SBATCH --ntasks 1
#SBATCH --array 1-150
#SBATCH --mail-type=END                       
#SBATCH --mail-user=xiucheng.yang@uconn.edu  

module load matlab/2019b
# Replace with your application's commands
cd /home/xiy19029/DECODE_1.1/
matlab -nodisplay -nosplash -singleCompThread -r 'batchDECODE($SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_TASK_MAX,"TileName")'\;exit\; 
echo 'Finished Matlab Code at '
