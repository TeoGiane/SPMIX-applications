#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=50:mem=128gb
#PBS -l walltime=48:00:00
#PBS -N BD_Scenario2-generate_plots
#PBS -o /u/gianella/log/BD_Scenario2-generate_plots.out
#PBS -e /u/gianella/log/BD_Scenario2-generate_plots.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
export SPMIX=/u/gianella/containers/spmix_dev.sif
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir $SPMIX $@
}
export -f in-apptainer

# Create log directory if it doesn't exists
mkdir -p log

# Execution in containerized environment
in-apptainer cook exec bd_highdim:generate_plots &> log/bd_highdim-generate_plots.log
