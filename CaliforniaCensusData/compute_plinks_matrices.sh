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
  'printf -v n "%03d" {1}; Rscript --vanilla ./R/compute_plinks_matrix.R -c {2} -r {3} -o ./summary/plinks {1} &> ./log/H_{2}/rho_{3}/compute_plinks_matrix_$n.log' \
  ::: $(eval echo {1..$NSIM}) ::: ${H[@]} ::: ${RHO[@]}

# Finished job notification
echo "JOB 'compute_plinks_matrices.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
