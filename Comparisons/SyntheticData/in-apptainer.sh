#!/bin/bash

APPTAINER=/opt/mox/apptainer/bin/apptainer

$APPTAINER exec --pwd /workdir --bind `pwd`:/workdir spmix.sif $@

