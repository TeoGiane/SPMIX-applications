from cook import create_task

# , Task
# from cook.contexts import create_group
# import itertools as it
# import numpy as np
# import os
# from pathlib import Path
# from typing import Literal

# Create generate-data task
generate_data_action = ["Rscript","--vanilla","src/generate_data.R"]
generate_data_deps = ["src/generate_data.R"]
generate_data_targets = ["input/data.dat", "input/shp/grid.dbf", "input/shp/grid.prj","input/shp/grid.shp", "input/shp/grid.shx"]
create_task(name="generate_data", dependencies=generate_data_deps, action=generate_data_action, targets=generate_data_targets)

# Create tasks for all algorithms
MODELS = ["SPMIX", "CARBayes", "naiveMCAR", "SKATER"]
run_model_targets = [f"output/{model}-fit.dat" for model in MODELS]
for model in MODELS:
    run_model_action = ["Rscript", "--vanilla", f"src/run_{model}.R"]
    run_model_deps = [f"src/run_{model}.R"] + generate_data_targets
    create_task(name=f"run_{model}", dependencies=run_model_deps, action=run_model_action, targets=[run_model_targets[MODELS.index(model)]])

# Create generate_plot task
generate_plot_action = ["Rscript","--vanilla","src/generate_plot.R"]
generate_plot_deps = ["src/generate_plot.R"] + run_model_targets
generate_plot_targets = [f"output/plt_boundaries-{model}.pdf" for model in MODELS] + [f"output/plt_plinks-{model}.pdf" for model in MODELS]
create_task(name="generate_plot", dependencies=generate_plot_deps, action=generate_plot_action, targets=generate_plot_targets)