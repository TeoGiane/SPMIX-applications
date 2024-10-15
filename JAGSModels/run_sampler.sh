#!/bin/bash

# Make log folders if not present
mkdir -p log/JAGSModel

# Execute run_sampler.R in parallel via GNU parallel
Rscript --vanilla ./R/newnew_multivariate_CAR.R &> ./log/JAGSModel/newnew_multivariate_CAR.log

# Finished job notification
echo "JOB 'run_sampler.sh' has finished" | mail -s "[Notification] - Job finished" matteo.gianella@polimi.it
