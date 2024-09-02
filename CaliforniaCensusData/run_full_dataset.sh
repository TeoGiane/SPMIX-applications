#!/bin/bash

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
  '; Rscript --vanilla ./R/run_sampler.R -c {1} -r {2} -o ./output/H_{1}/rho_{2}/chain_full_dataset.dat ./data/full_dataset.dat &> ./log/H_{1}/rho_{2}/run_full_dataset.log' \
  ::: ${H[@]} ::: ${RHO[@]}

# Finished job notification
echo "JOB 'run_full_dataset.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
