#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=64:mem=128gb
#PBS -l walltime=48:00:00
#PBS -N BD_Scenario1-simulation_study:run
#PBS -o /u/gianella/log/BD_Scenario1-simulation_study:run.out
#PBS -e /u/gianella/log/BD_Scenario1-simulation_study:run.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
export SPMIX=/u/gianella/containers/spmix.sif
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir $SPMIX $@
}
export -f in-apptainer

# Execution in containerized environment
in-apptainer cook exec --jobs 96 simulation_study:run &> log/simulation_study:run.log
