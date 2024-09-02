#!/bin/bash

# Define number of datasets
NSIM=10

# Set values for H and RHO
H=(RJ)
RHO=(0 0.9 0.95 0.99)

# Make log folders if not present
for h in ${H[@]}; do
  for r in ${RHO[@]}; do
    mkdir -p log/H_$h/rho_$r
  done
done

# Execute run_sampler.R in parallel via GNU parallel
parallel -j 15 \
  'printf -v n "%03d" {1}; Rscript --vanilla ./R/run_sampler.R -c {2} -r {3} -o ./output/H_{2}/rho_{3}/chain_$n.dat ./data/data_$n.dat &> ./log/H_{2}/rho_{3}/run_sampler_$n.log' \
  ::: $(eval echo {1..$NSIM}) ::: ${H[@]} ::: ${RHO[@]}

# Finished job notification
echo "JOB 'run_sim_study.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
