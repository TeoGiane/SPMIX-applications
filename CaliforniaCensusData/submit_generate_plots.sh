#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -l walltime=144:00:00
#PBS -N SPMIX-CaliforniaCensusData-run_full_dataset
#PBS -o log/SPMIX-CaliforniaCensusData-run_full_dataset.out
#PBS -e log/SPMIX-CaliforniaCensusData-run_full_dataset.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir spmix.sif $@
}
export -f in-apptainer

# Execution in containerized environment
in-apptainer cook exec generate_plot &> generate_plot.log
