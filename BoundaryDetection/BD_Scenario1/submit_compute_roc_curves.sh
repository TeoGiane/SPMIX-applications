#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=50:mem=128gb
#PBS -l walltime=48:00:00
#PBS -N BD_Scenario1-compute_roc_curves
#PBS -o /u/gianella/log/BD_Scenario1-compute_roc_curves.out
#PBS -e /u/gianella/log/BD_Scenario1-compute_roc_curves.err

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
in-apptainer cook exec simulation_study:compute_roc_curves &> log/compute_roc_curves.log
