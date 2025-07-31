#!/bin/bash -x
#SBATCH --job-name=u1potential2p1d
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j-error.out


export OMP_NUM_THREAD=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY=balanced,granularity=fine,verbose

date

## load the modules needed for su2 library, mostly GCC and GCCcore
## replace by your installation path
srun /hiskp4/gross/masterthesis/codenew/software-stack/install/su2-su3_heatbath_overrelaxation/bin/u1-main -f $1

date


