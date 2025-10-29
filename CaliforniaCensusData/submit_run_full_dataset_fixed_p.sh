#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=24:mem=64gb
#PBS -l walltime=240:00:00
#PBS -N SPMIX-CaliforniaCensusData-run_full_dataset_fixed_p
#PBS -o log/SPMIX-CaliforniaCensusData-run_full_dataset_fixed_p.out
#PBS -e log/SPMIX-CaliforniaCensusData-run_full_dataset_fixed_p.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Load required modules
source /opt/mox/spack/share/spack/setup-env.sh
spack load parallel

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
export SPMIX=/u/gianella/containers/spmix.sif
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir $SPMIX $@
}
export -f in-apptainer

# Define graph_sparsity values
SPARSITIES=(0.1 0.2 0.3 0.4 0.5)

# Generate log directories
for SPARSITY in "${SPARSITIES[@]}"; do
	mkdir -p log/HRJ/rho0.95/alpha24_beta22/p_${SPARSITY}
done

# Parallel execution in containerized environment
parallel -j 0 'in-apptainer cook exec run_full_dataset_fixed_p_{1} &> log/HRJ/rho0.95/alpha24_beta22/p_{1}/run_full_dataset_fixed_p.log' ::: "${SPARSITIES[@]}"


