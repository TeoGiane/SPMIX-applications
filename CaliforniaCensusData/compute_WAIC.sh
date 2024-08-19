#!/bin/bash

# Set values for H
H=(RJ)

# Make log folders if not present
for h in ${H[@]}; do
    mkdir -p log/H_$h
done

# Execute run_sampler.R in parallel via GNU parallel
parallel -j 15 \
    'Rscript --vanilla ./R/compute_WAIC.R -c {1} ./summary/meanWAIC-H_{1}.csv &> ./log/H_{1}/compute_WAIC.log' ::: ${H[@]}

# Finished job notification
echo "JOB 'compute_meanWAIC.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
