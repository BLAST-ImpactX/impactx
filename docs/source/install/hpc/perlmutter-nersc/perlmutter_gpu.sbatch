#!/bin/bash -l

# Copyright 2021-2023 Axel Huebl, Kevin Gott
#
# This file is part of ImpactX.
#
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH -J ImpactX
#    note: <proj> must end on _g
#SBATCH -A <proj>
#SBATCH -q regular
# A100 40GB (most nodes)
#SBATCH -C gpu
# A100 80GB (256 nodes)
#S BATCH -C gpu&hbm80g
#SBATCH --exclusive
# ideally single:1, but NERSC cgroups issue
#SBATCH --gpu-bind=none
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH -o ImpactX.o%j
#SBATCH -e ImpactX.e%j

# python interpreter & script here
EXE=python3
INPUTS=run_fodo.py
# or executable & inputs file
#EXE=./impactx
#INPUTS=input_fodo.in

# enviroment setup
if [[ -z "${MY_PROFILE}" ]]; then
    echo "WARNING: FORGOT TO"
    echo "   source $HOME/perlmutter_gpu_impactx.profile"
    echo "before submission. Doing that now."

    source $HOME/perlmutter_gpu_impactx.profile
fi

# pin to closest NIC to GPU
export MPICH_OFI_NIC_POLICY=GPU

# threads for OpenMP and threaded compressors per MPI rank
#   note: 16 avoids hyperthreading (32 virtual cores, 16 physical)
export SRUN_CPUS_PER_TASK=16

# GPU-aware MPI optimizations
export AMREX_DEFAULT_INIT="amrex.use_gpu_aware_mpi=1"

# CUDA visible devices are ordered inverse to local task IDs
#   Reference: nvidia-smi topo -m
srun --cpu-bind=cores bash -c "
    export CUDA_VISIBLE_DEVICES=\$((3-SLURM_LOCALID));
    ${EXE} ${INPUTS}" \
  > output.txt
