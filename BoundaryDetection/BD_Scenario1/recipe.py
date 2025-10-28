from cook import create_task, Task
from cook.contexts import create_group

# Define quantities for simulation study
num_datasets = 50
ids = [f"{i:03d}" for i in range(1, num_datasets+1)]
num_components_values = [2, 4, 6, 8, 10, "RJ"]
rho_values = [0, 0.5, 0.9, 0.95, 0.99]

# Define generate_data task
generate_data_action = ["Rscript", "src/generate_data.R"] + \
    ["--num-datasets", num_datasets] + \
    ["--dest-dir", "input"]
generate_data_targets = [f"input/data_{i:03d}.dat" for i in range(1, num_datasets+1)]
create_task("simulation_study:generate_data", action = generate_data_action, targets = generate_data_targets)

# Function to create single run_sampler task in simulation study (hidden)
def create_run_sampler_task(dataset_id: str, rho: float, num_components: int | str, output_file: str) -> Task:
    input_file = f"input/data_{dataset_id}.dat"
    action = ["Rscript", "src/run_sampler.R"] + \
        ["--num-components", num_components] + \
        ["--rho", rho] + \
        ["--output-file", output_file] + \
        [input_file]
    deps = []# [input_file]
    targets = []# [output_file]
    return create_task(name = f"_simulation-study:run-H{num_components}-rho{rho}-replica{dataset_id}",
                       action = action, targets = targets, dependencies = deps)

# Create run_sampler task group
with create_group("simulation_study:run") as simulation_study_group:
    for num_components in num_components_values:
        for rho in rho_values:
            for id in ids:
                output_file = f"output/H{num_components}/rho{rho}/chain_{id}.dat"
                create_run_sampler_task(id, rho, num_components, output_file)

# Function to create compute_confusion_matrices task in simulation study (hidden)
def create_confusion_matrices_task(num_datasets: int, num_components: int | str, rho: float) -> Task:
    output_file = f"summary/confusion_matrices-H{num_components}-rho{rho}.csv"
    action = ["Rscript", "src/compute_confusion_matrices.R"] + \
        ["--num-datasets", num_datasets] + \
        ["--num-components", num_components] + \
        ["--rho", rho] + \
        [output_file]
    deps = []# [input_file]
    targets = []# [output_file]
    return create_task(name = f"_simulation-study:compute_confusion-matrices-H{num_components}-rho{rho}",
                       action = action, targets = targets, dependencies = deps)

# Create compute_confusion_matrices task group
with create_group("simulation_study:compute_confusion_matrices") as confusion_matrices_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_confusion_matrices_task(num_datasets, num_components, rho)

# Function to create compute_mean_L1_distances task in simulation study (hidden)
def create_mean_L1_distances_task(num_datasets: int, num_components: int | str, rho: float) -> Task:
    output_file = f"summary/mean_L1_distances-H{num_components}-rho{rho}.csv"
    action = ["Rscript", "src/compute_mean_L1_distances.R"] + \
        ["--num-datasets", num_datasets] + \
        ["--num-components", num_components] + \
        ["--rho", rho] + \
        [output_file]
    deps = []# [input_file]
    targets = []# [output_file]
    return create_task(name = f"_simulation-study:compute-mean-L1-distances-H{num_components}-rho{rho}",
                       action = action, targets = targets, dependencies = deps)

# Create compute_mean_L1_distances task group
with create_group("simulation_study:compute_mean_L1_distances") as mean_L1_distances_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_mean_L1_distances_task(num_datasets, num_components, rho)

# Function to create compute_WAIC task in simulation study (hidden)
def create_WAIC_task(num_datasets: int, num_components: int | str, rho: float) -> Task:
    output_file = f"summary/WAIC-H{num_components}-rho{rho}.csv"
    action = ["Rscript", "src/compute_WAIC.R"] + \
        ["--num-datasets", num_datasets] + \
        ["--num-components", num_components] + \
        ["--rho", rho] + \
        [output_file]
    deps = []# [input_file]
    targets = []# [output_file]
    return create_task(name = f"_simulation-study:compute-WAIC-H{num_components}-rho{rho}",
                       action = action, targets = targets, dependencies = deps)

# Create compute_WAIC task group
with create_group("simulation_study:compute_WAIC") as WAIC_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_WAIC_task(num_datasets, num_components, rho)
