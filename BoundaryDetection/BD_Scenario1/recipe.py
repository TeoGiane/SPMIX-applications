from cook import create_task, Task
from cook.contexts import create_group

# Define quantities for simulation study
num_datasets = 50
ids = [f"{i:03d}" for i in range(1, num_datasets+1)]
num_components_values = [2, 4, 6, 8, 10] # + "RJ"
poisson_rate_values = [1.0]#, 2.0, 5.0, 10.0]
rho_values = [0, 0.5, 0.9, 0.95, 0.99]

# Define generate_data task
generate_data_action = ["Rscript", "src/generate_data.R"] + \
    ["--num-datasets", num_datasets] + \
    ["--dest-dir", "input"]
generate_data_targets = [f"input/data{i:03d}.dat" for i in range(1, num_datasets+1)]
create_task("simulation_study:generate_data", action = generate_data_action, targets = generate_data_targets)

# Function to create single run_sampler task in simulation study (hidden)
def create_run_sampler_task(dataset_id: str, rho: float, num_components: int | str, poisson_rate: float, output_file: str) -> Task:
    input_file = f"input/data_{dataset_id}.dat"
    if num_components == "RJ":
        task_name = f"_simulation-study:run-HRJ-poisson{poisson_rate}-rho{rho}-replica{dataset_id}"
        action = ["Rscript", "src/run_sampler.R"] + \
            ["--rho", rho] + \
            ["--num-components", "RJ"] + \
            ["--poisson-rate", poisson_rate] + \
            ["--output-file", output_file] + \
            [input_file]
    else:
        task_name = f"_simulation-study:run-H{num_components}-rho{rho}-replica{dataset_id}"  
        action = ["Rscript", "src/run_sampler.R"] + \
            ["--rho", rho] + \
            ["--num-components", num_components] + \
            ["--output-file", output_file] + \
            [input_file]
    deps = [] # [simulation_study:generate_data"]
    targets = [] # [output_file]
    return create_task(name = task_name, action = action, targets = targets, task_dependencies = deps)

# Create run_sampler task group
with create_group("simulation_study:run") as simulation_study_group:
    for num_components in num_components_values:
        for rho in rho_values:
            for id in ids:
                output_file = f"output/H{num_components}/rho{rho}/chain{id}.dat"
                create_run_sampler_task(id, rho, num_components, None, output_file)
    for poisson_rate in poisson_rate_values:
        for rho in rho_values:
            for id in ids:
                output_file = f"output/HRJ/poisson{poisson_rate}/rho{rho}/chain{id}.dat"
                create_run_sampler_task(id, rho, "RJ", poisson_rate, output_file)


# Function to create compute_confusion_matrices task in simulation study (hidden)
def create_confusion_matrices_task(num_datasets: int, num_components: int | str, poisson_rate: float, rho: float) -> Task:
    if num_components == "RJ":
        task_name = f"_simulation-study:compute_confusion_matrices-HRJ-poisson{poisson_rate}-rho{rho}"
        output_file = f"summary/confusion_matrices-HRJ-poisson{poisson_rate}-rho{rho}.csv"
        action = ["Rscript", "src/compute_confusion_matrices.R"] + \
            ["--num-datasets", num_datasets] + \
            ["--num-components", "RJ"] + \
            ["--poisson-rate", poisson_rate] + \
            ["--rho", rho] + \
            [output_file]
    else:
        task_name = f"_simulation-study:compute_confusion_matrices-H{num_components}-rho{rho}"
        output_file = f"summary/confusion_matrices-H{num_components}-rho{rho}.csv"
        action = ["Rscript", "src/compute_confusion_matrices.R"] + \
            ["--num-datasets", num_datasets] + \
            ["--num-components", num_components] + \
            ["--rho", rho] + \
            [output_file]
    deps = [] # [simulation_study:run"]
    targets = [] # [output_file]
    return create_task(name = task_name, action = action, targets = targets, task_dependencies = deps)

# Create compute_confusion_matrices task group
with create_group("simulation_study:compute_confusion_matrices") as confusion_matrices_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_confusion_matrices_task(num_datasets, num_components, None, rho)
    for poisson_rate in poisson_rate_values:
        for rho in rho_values:
            create_confusion_matrices_task(num_datasets, "RJ", poisson_rate, rho)


# Function to create compute_mean_L1_distances task in simulation study (hidden)
def create_mean_L1_distances_task(num_datasets: int, num_components: int | str, poisson_rate: float, rho: float) -> Task:
    if num_components == "RJ":
        task_name = f"_simulation-study:compute_mean_L1_distances-HRJ-poisson{poisson_rate}-rho{rho}"
        output_file = f"summary/mean_L1_distances-HRJ-poisson{poisson_rate}-rho{rho}.csv"
        action = ["Rscript", "src/compute_mean_L1_distances.R"] + \
            ["--num-datasets", num_datasets] + \
            ["--num-components", "RJ"] + \
            ["--poisson-rate", poisson_rate] + \
            ["--rho", rho] + \
            [output_file]
    else:
        task_name = f"_simulation-study:compute_mean_L1_distances-H{num_components}-rho{rho}"
        output_file = f"summary/mean_L1_distances-H{num_components}-rho{rho}.csv"
        action = ["Rscript", "src/compute_mean_L1_distances.R"] + \
            ["--num-datasets", num_datasets] + \
            ["--num-components", num_components] + \
            ["--rho", rho] + \
            [output_file]
    deps = [] # ["simulation_study:run"]
    targets = [] # [output_file]
    return create_task(name = task_name, action = action, targets = targets, task_dependencies = deps)

# Create compute_mean_L1_distances task group
with create_group("simulation_study:compute_mean_L1_distances") as mean_L1_distances_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_mean_L1_distances_task(num_datasets, num_components, None, rho)
    for poisson_rate in poisson_rate_values:
        for rho in rho_values:
            create_mean_L1_distances_task(num_datasets, "RJ", poisson_rate, rho)


# Function to create compute_WAIC task in simulation study (hidden)
def create_WAIC_task(num_datasets: int, num_components: int | str, poisson_rate: float, rho: float) -> Task:
    if num_components == "RJ":
        task_name = f"_simulation-study:compute_WAIC-HRJ-poisson{poisson_rate}-rho{rho}"
        output_file = f"summary/WAIC-HRJ-poisson{poisson_rate}-rho{rho}.csv"
        action = ["Rscript", "src/compute_WAIC.R"] + \
            ["--num-datasets", num_datasets] + \
            ["--num-components", "RJ"] + \
            ["--poisson-rate", poisson_rate] + \
            ["--rho", rho] + \
            [output_file]
    else:
        task_name = f"_simulation-study:compute_WAIC-H{num_components}-rho{rho}"
        output_file = f"summary/WAIC-H{num_components}-rho{rho}.csv"
        action = ["Rscript", "src/compute_WAIC.R"] + \
            ["--num-datasets", num_datasets] + \
            ["--num-components", num_components] + \
            ["--rho", rho] + \
            [output_file]
    deps = [] # [simulation_study:run"]
    targets = [] # [output_file]
    return create_task(name = task_name, action = action, targets = targets, task_dependencies = deps)

# Create compute_WAIC task group
with create_group("simulation_study:compute_WAIC") as WAIC_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_WAIC_task(num_datasets, num_components, None, rho)
    for poisson_rate in poisson_rate_values:
        for rho in rho_values:
            create_WAIC_task(num_datasets, "RJ", poisson_rate, rho)


# Function to create compute_roc_curves task in simulation study (hidden)
def create_roc_curves_task(num_datasets: int, num_components: int | str, rho: float) -> Task:
    output_file = f"summary/ROC_curves-H{num_components}-rho{rho}.csv"
    action = ["Rscript", "src/compute_roc_curves.R"] + \
        ["--num-datasets", num_datasets] + \
        ["--num-components", num_components] + \
        ["--rho", rho] + \
        [output_file]
    deps = [] # [simulation_study:run"]
    targets = [] # [output_file]
    return create_task(name = f"_simulation-study:compute_roc_curves-H{num_components}-rho{rho}",
                       action = action, targets = targets, task_dependencies = deps)

# Create compute_roc_curves task group
with create_group("simulation_study:compute_roc_curves") as roc_curves_group:
    for num_components in num_components_values:
        for rho in rho_values:
            create_roc_curves_task(num_datasets, num_components, rho)

# Create generate_tables task
generate_tables_action = ["Rscript", "src/generate_tables.R"] + \
    ['--summary-path', 'summary'] + \
    ['--num-components-values', ",".join(str(nc) for nc in num_components_values)] + \
    ['--rho-values', ",".join(str(r) for r in rho_values)] + \
    ['--output-dir', 'tables']
generate_tables_deps = [] # ["simulation_study:compute_confusion_matrices","simulation_study:compute_mean_L1_distances","simulation_study:compute_WAIC"]
generate_tables_targets = []
create_task("simulation_study:generate_tables",
            action = generate_tables_action, targets = generate_tables_targets, task_dependencies = generate_tables_deps)