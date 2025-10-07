#!/bin/bash

# PBS Settings
#PBS -S /bin/bash
#PBS -l select=1:ncpus=96
#PBS -l walltime=504:00:00
#PBS -N SyntheticData-comparative_study
#PBS -o log/SyntheticData-comparative_study.out
#PBS -e log/SyntheticData-comparative_study.err

# Set current working directory
cd ${PBS_O_WORKDIR}

# Load required modules
source /opt/mox/spack/share/spack/setup-env.sh
spack load parallel

# Create alias for apptainer execution
export APPTAINER=/opt/mox/apptainer/bin/apptainer
in-apptainer () {
	$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir spmix.sif $@
}
export -f in-apptainer

# Define models to compare
MODELS=(SPMIX CARBayes naiveMCAR SKATER)

# Create log directory if it doesn't exist
mkdir -p log

# Run comparative study for each model
parallel -j 0 'in-apptainer cook exec run_{1} &> log/run_{1}.log' ::: "${MODELS[@]}"
in-apptainer cook exec generate_plot &> log/generate_plot.log
