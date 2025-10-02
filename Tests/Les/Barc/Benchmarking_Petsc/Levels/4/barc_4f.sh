#!/bin/bash

#SBATCH --clusters=merlin6
#SBATCH --job-name=Barc_4F
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=out_barc_4f
#SBATCH --mail-type=ALL                  # mail events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --mail-user=bojan.niceno@psi.ch  # where to send mail

module purge
module load gcc/10.3.0 openmpi/4.0.5-2_slurm

source ~/.bashrc

cd /data/user/niceno/Barc/T-Flows/Tests/Les/Barc/Benchmarking_Petsc/Levels/4/Run_F

mpirun -np 16 ./Process

