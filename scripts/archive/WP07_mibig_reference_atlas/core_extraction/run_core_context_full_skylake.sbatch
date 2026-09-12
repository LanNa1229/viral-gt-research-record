#!/bin/bash
#SBATCH --job-name=mibig_core_full
#SBATCH --account=nknlab
#SBATCH --partition=batch
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --constraint=Skylake
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/ln72030/BGC_Prophet/core_context_full_20260910/full_retry_%j.out
#SBATCH --error=/scratch/ln72030/BGC_Prophet/core_context_full_20260910/full_retry_%j.err
set -euo pipefail
cd /scratch/ln72030/BGC_Prophet
/scratch/ln72030/conda_envs/bgcprophet/bin/python -B /scratch/ln72030/BGC_Prophet/core_context_full_20260910/run_hpc_core_context_full.py --work /scratch/ln72030/BGC_Prophet --out /scratch/ln72030/BGC_Prophet/core_context_full_20260910/results_validatednode

