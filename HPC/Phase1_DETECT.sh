#!/bin/bash
##SBATCH --partition general

#SBATCH --account=zhz18039
#SBATCH --partition=priority
#SBATCH --qos=zhz18039epyc
#SBATCH --constraint=epyc128

#SBATCH --ntasks 1 
#SBATCH --array 1-100
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=END                       
#SBATCH --mail-user=xiucheng.yang@uconn.edu  

echo $SLURMD_NODENAME  # display the node name


module load matlab
# Replace with your application's commands
cd /home/xiy19029/DECODE_v2_Share/
matlab -nodisplay -nosplash -singleCompThread -r "batchDECODE_Phase1_DETECT($SLURM_ARRAY_TASK_ID, $SLURM_ARRAY_TASK_MAX,'h025v017');exit"; 
echo 'Finished Matlab Code at '

