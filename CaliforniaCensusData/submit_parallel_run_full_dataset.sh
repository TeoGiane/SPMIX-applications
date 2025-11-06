#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=24:mem=64gb
#PBS -l walltime=240:00:00
#PBS -N SPMIX-CaliforniaCensusData-parallel_run_full_dataset
#PBS -o log/SPMIX-CaliforniaCensusData-parallel_run_full_dataset.out
#PBS -e log/SPMIX-CaliforniaCensusData-parallel_run_full_dataset.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Create log directory
mkdir -p log

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
export SPMIX=/u/gianella/containers/spmix_dev.sif
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir $SPMIX $@
}
export -f in-apptainer

# Execution in containerized environment
in-apptainer cook exec --jobs 24 parallel_run_full_dataset &> log/parallel_run_full_dataset.log

