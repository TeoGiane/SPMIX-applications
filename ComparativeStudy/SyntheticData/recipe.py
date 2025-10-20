from cook import create_task

# Create generate-data task
generate_data_deps = ["src/generate_data.R"]
generate_data_targets = ["input/data.dat", "input/shp/grid.dbf", "input/shp/grid.prj","input/shp/grid.shp", "input/shp/grid.shx"]
create_task(name="generate_data", action=["Rscript","src/generate_data.R"],
            dependencies=generate_data_deps,
            targets=generate_data_targets)

# Create tasks for all algorithms
MODELS = ["SPMIX", "CARBayes", "naiveMCAR", "SKATER"]
run_model_targets = [f"output/{model}-fit.dat" for model in MODELS]
for model in MODELS:
    run_model_deps = [f"src/run_{model}.R"] + generate_data_targets
    create_task(name=f"run_{model}", action=["Rscript", f"src/run_{model}.R"],
                dependencies=run_model_deps,
                targets=[run_model_targets[MODELS.index(model)]])

# Create generate_plot task
generate_plot_deps = ["src/generate_plot.R"] + run_model_targets
generate_plot_targets = [f"plots/plt_boundaries-{model}.pdf" for model in MODELS]
create_task(name="generate_plot", action=["Rscript","src/generate_plot.R"],
            dependencies=generate_plot_deps,
            targets=generate_plot_targets)