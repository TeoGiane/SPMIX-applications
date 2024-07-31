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
parallel --dry-run -j 6 \
    'Rscript --vanilla run_sampler.R -c {1} -r {2} -o ./output/H_{1}/rho_{2}/chain.dat &> ./log/H_{1}/rho_{2}/run_sampler_$n.log' \
    ::: ${H[@]} ::: ${RHO[@]}
