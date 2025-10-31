#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=1:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N SPMIX-CaliforniaCensusData-generate_plots
#PBS -o log/SPMIX-CaliforniaCensusData-generate_plots.out
#PBS -e log/SPMIX-CaliforniaCensusData-generate_plots.err

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
in-apptainer cook exec generate_plot &> log/HRJ/rho0.95/alpha24_beta22/a2_b93/generate_plot.log
