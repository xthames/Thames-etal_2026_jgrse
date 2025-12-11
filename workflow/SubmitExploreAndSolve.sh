#!/bin/bash
#SBATCH --mail-user=USER@domain.ext
#SBATCH --mail-type=ALL
#SBATCH --partition=open
#SBATCH --job-name=varfRt_ExploreSOWs
#SBATCH --output=/storage/home/USER/work/research/VarRegassing/jobs/output/%x.out
#SBATCH --error=/storage/home/USER/work/research/VarRegassing/jobs/error/%x.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8


# (0a) load in the appropriate modules, switch to the right directory
module load matlab/R2023a
module load parallel
cd /storage/home/USER/work/research/VarRegassing


# (0b) environment variable names: ORDER IMPORTANT!
## $1 = script to use (general or target)
## $2 = timeDirection
## $3 = viscDep
## $4 = viscStrength
## $5 = modelType
## $6 = varyfR (optional, if modelType = 'var')
## $7 = varyphiRum (optional, if modelType = 'var')
if [ "$1" == "targeted" ]
then
	MATLAB_SCRIPT="ExploreAndSolveTargetedTXum_2025XI"
else
	MATLAB_SCRIPT="ExploreAndSolveGeneralTXum_2025XI"
fi
IS_SAVED_STR="MODEL NOT SAVED..."


# (1) loop the exploration of initial SOWs for Earth with varfRt until we get some successes 
while [ "$IS_SAVED_STR" == "MODEL NOT SAVED..." ]
do
	> jobs/output/varfRt_ExploreSOWs.out	
	SRUN_OPTS="--nodes=1 --ntasks=1 --cpus-per-task=$SLURM_CPUS_PER_TASK"
	srun $SRUN_OPTS matlab -batch "timeDirection='$2'; viscDep='$3'; viscStrength='$4'; modelType='$5'; varyfR='$6'; varyphiRum='$7'; $MATLAB_SCRIPT" &
	wait
	IS_SAVED_STR=$(cat jobs/output/varfRt_ExploreSOWs.out | tail -n 2 | head -n 1)
done

