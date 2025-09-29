from cook import create_task # , Task
# from cook.contexts import create_group
# import itertools as it
# import numpy as np
# import os
# from pathlib import Path
# from typing import Literal

# Create requirements task
create_task("requirements", action="pip-compile -v", targets=["requirements.txt"],
            dependencies=["requirements.in"])

# Create pip-sync task
create_task("pip-sync", action="pip-sync", dependencies=["requirements.txt"])

# Create generate-data task
action = ["Rscript","--vanilla","src/generate_data.R"]
deps = ["src/generate_data.R"]
input_targets = ["input/data.dat", "input/shp/grid.dbf", "input/shp/grid.prj","input/shp/grid.shp", "input/shp/grid.shx"]
create_task(name="generate_data", dependencies=deps, action=action, targets=input_targets)

# Create tasks for all algorithms
MODELS = ["SPMIX", "CARBayes", "naiveMCAR", "SKATER"]
for model in MODELS:
    action = ["Rscript","--vanilla",f"src/run_{model}.R"]
    deps = [f"src/run_{model}.R"] + input_targets
    create_task(name=f"run_{model}", dependencies=deps, action=action)
