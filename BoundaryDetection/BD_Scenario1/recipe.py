from cook import create_task, Task
from cook.contexts import create_group

# Define generate_data task
num_datasets = 50
generate_data_action = ["Rscript", "src/generate_data.R"] + \
    ["--num-datasets", num_datasets] + \
    ["--dest-dir", "input"]
generate_data_targets = [f"input/data_{i:03d}.dat" for i in range(1, num_datasets+1)]
create_task("simulation_study:generate_data", action = generate_data_action, targets = generate_data_targets)

# Function to create single run_sampler task in simulation study (hidden)
def create_simulation_study_task(dataset_id: str, rho: float, num_components: int | str, output_file: str) -> Task:
    input_file = f"input/data_{dataset_id}.dat"
    action = ["Rscript", "src/run_sampler.R"] + \
        ["--num-components", num_components] + \
        ["--rho", rho] + \
        ["--output-file", output_file] + \
        [input_file]
    deps = []#[input_file]
    targets = []#[output_file]
    return create_task(name = f"_simulation-study:run-H{num_components}-rho{rho}-replica{dataset_id}",
                       action = action, targets = targets, dependencies = deps)

# Create simulation study task group
with create_group("simulation_study:run") as simulation_study_group:
    num_components_values = [2, 4, 6, 8, 10, "RJ"]
    rho_values = [0, 0.5, 0.9, 0.95, 0.99]
    ids = [f"{i:03d}" for i in range(1, num_datasets+1)]
    for num_components in num_components_values:
        for rho in rho_values:
            for id in ids:
                output_file = f"output/H{num_components}/rho{rho}/chain_{id}.dat"
                create_simulation_study_task(id, rho, num_components, output_file)
