#!/bin/bash
#SBATCH --partition generalsky
##SBATCH --partition=EpycPriority
##SBATCH --account=zhz18039
##SBATCH --partition=OSGPriority
##SBATCH --account=osgusers
#SBATCH --ntasks 2
#SBATCH --array 1-20
#SBATCH --mail-type=END                       
#SBATCH --mail-user=xiucheng.yang@uconn.edu  

module load matlab/2019b
# Replace with your application's commands
cd /home/xiy19029/DECODE_1.1/
matlab -nodisplay -nosplash -singleCompThread -r 'batchCoverChangeMaps($SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_TASK_MAX)'\;exit\; 
echo 'Finished Matlab Code at '