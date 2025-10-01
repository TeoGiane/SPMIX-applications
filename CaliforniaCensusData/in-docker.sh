#!/bin/bash

docker run --rm -it -v `pwd`:/workdir spmix $@
