#!/bin/bash

# Make log folders if not present
mkdir -p log

# Execute run_sampler.R in parallel via GNU parallel
Rscript --vanilla ./R/SPMIXvsMCAR.R -a 2 -b 36 -o "output" &> ./log/SPMIX_vs_MCAR.log

# Finished job notification
echo "JOB 'SPMIX_vs_MCAR.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
