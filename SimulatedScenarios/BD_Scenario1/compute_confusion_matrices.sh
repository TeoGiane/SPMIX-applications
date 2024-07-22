#!/bin/bash

# Define number of datasets
NSIM=50

# Set values for H and RHO
H=(2 4 6 8 10 RJ)
RHO=(0 0.5 0.9 0.95 0.99)

# Make log folders if not present
for h in ${H[@]}; do
    for r in ${RHO[@]}; do
        mkdir -p log/H_$h/rho_$r
    done
done

# Execute run_sampler.R in parallel via GNU parallel
parallel -j 6 \
    'Rscript --vanilla compare_boundary_graph.R -d {1} -c {2} -r {3} ./summary/CM-H_{2}-rho_{3}.csv &> ./log/H_{2}/rho_{3}/compare_boundary_graph.log' \
    ::: $NSIM ::: ${H[@]} ::: ${RHO[@]}
