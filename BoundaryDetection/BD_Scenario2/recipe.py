from cook import create_task, Task
from cook.contexts import create_group
# import requests
# import zipfile
# from pathlib import Path

# import re
# import os
# from typing import Optional

# Define quantities for simulation study
num_datasets = 1
ids = [f"{i:03d}" for i in range(1, num_datasets+1)]

# Define generate_data task
generate_data_action = ["Rscript", "src/generate_data.R"] + \
    ["--num-datasets", num_datasets] + \
    ["--dest-dir", "input"]
generate_data_targets = [f"input/data{i:03d}.dat" for i in range(1, num_datasets+1)]
create_task("bd_highdim:generate_data", action = generate_data_action, targets = generate_data_targets)

# Define run_sampler tasks
with create_group("bd_highdim:run_sampler") as run_sampler_group:
    for id in ids:
        run_sampler_action = ["Rscript", "src/run_sampler.R"] + \
            [f"input/data{id}.dat"] + \
            ["--output-file", f"output/chain{id}.dat"]
        create_task(f"_bd_highdim:run_sampler-{id}", action=run_sampler_action,
                    task_dependencies=["bd_highdim:generate_data"],
                    targets=[f"output/chain{id}.dat"])
        
# Define generate_plot tasks
with create_group("bd_highdim:generate_plots") as generate_plots_group:
    for id in ids:
        generate_plot_action = ["Rscript", "src/generate_plots.R"] + \
            ["--data-file", f"input/data{id}.dat"] + \
            ["--sim-file", f"output/chain{id}.dat"] + \
            ["--output-dir", f"plots/chain{id}"]
        create_task(f"_bd_highdim:generate_plots-{id}", action=generate_plot_action,
                    task_dependencies=[f"_bd_highdim:run_sampler-{id}"])
