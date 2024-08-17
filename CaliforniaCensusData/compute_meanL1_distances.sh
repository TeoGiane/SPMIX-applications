#!/bin/bash

# Set values for H and RHO
H=(RJ)

# Make log folders if not present
for h in ${H[@]}; do
		mkdir -p log/H_$h
done

# Execute run_sampler.R in parallel via GNU parallel
parallel -j 15 \
    'Rscript --vanilla ./R/compare_densities.R -c {1} ./summary/meanL1-H_{1}.csv &> ./log/H_{1}/compare_densities.log' ::: ${H[@]}

# Finished job notification
echo "JOB 'compute_meanL1_distances.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
