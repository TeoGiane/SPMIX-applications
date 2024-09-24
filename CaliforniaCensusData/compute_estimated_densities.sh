#!/bin/bash

# Define number of datasets
NSIM=100

# Set values for H and RHO
H=(RJ)
RHO=(0.95)

# Make log folders if not present
for h in ${H[@]}; do
  for r in ${RHO[@]}; do
    mkdir -p log/H_$h/rho_$r
  done
done

# Execute run_sampler.R in parallel via GNU parallel
parallel -j 0 \
  'printf -v n "%03d" {1}; Rscript --vanilla ./R/compute_estimated_densities.R -c {2} -r {3} -o ./summary/mean_est_dens {1} &> ./log/H_{2}/rho_{3}/compute_estimated_densitites_$n.log' \
  ::: $(eval echo {1..$NSIM}) ::: ${H[@]} ::: ${RHO[@]}

# Finished job notification
echo "JOB 'compute_estimated_densities.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
