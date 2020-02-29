#!/bin/bash
#
# SLURM resource specifications
# (use an extra '#' in front of SBATCH to comment-out any unused options)
#

#SBATCH --job-name=ballistic_trajectory   # shows up in the output of 'squeue'
#SBATCH -o bt_out_%J         # create outfile
#SBATCH -e bt_err_%J         # create errorfile
#SBATCH --time=1-00:00:00        # specify the requested wall-time
#SBATCH --partition=dark         # specify the partition to run on
#SBATCH --nodes=1                # number of nodes allocated for this job
#SBATCH --ntasks-per-node=1      # number of MPI ranks per node, 20 for astro_short and astro_long
#SBATCH --cpus-per-task=1        # number of OpenMP threads per MPI rank
##SBATCH --exclude=<node list>   # avoid nodes (e.g. --exclude=node786)

#SBATCH --mail-user=cliu205@ucsc.edu # please put your own email, I don't want your error messages
#SBATCH --mail-type=ALL

# Load default settings for environment variables
source /users/software/astro/startup.d/modules.sh

# If required, replace specific modules
# module unload intelmpi
# module load mvapich2

# When compiling remember to use the same environment and modules

### set executable
EXE="Three_Body"
### your parameterfile
ARGS=""
### INFILE = 1 when restarting gadget from restart files, nothing when you start with IC, 2 when restarting from a snapshot
INFILE=""

# Execute the code
srun --mpi=pmi2 $EXE
