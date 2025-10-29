#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=24:mem=64gb
#PBS -l walltime=240:00:00
#PBS -N SPMIX-CaliforniaCensusData-run_full_dataset
#PBS -o log/SPMIX-CaliforniaCensusData-run_full_dataset.out
#PBS -e log/SPMIX-CaliforniaCensusData-run_full_dataset.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Create log directory
mkdir -p log/HRJ/rho0.95/alpha24_beta22/a2_b93/

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
export SPMIX=/u/gianella/containers/spmix.sif
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir $SPMIX $@
}
export -f in-apptainer

# Execution in containerized environment
in-apptainer cook exec run_full_dataset &> log/HRJ/rho0.95/alpha24_beta22/a2_b93/run_full_dataset.log
