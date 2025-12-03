from cook import create_task, Task
from cook.contexts import create_group
import requests
import zipfile
from pathlib import Path

import re
import os
from typing import Optional

def _clean_prior_string(prior_str: str, pattern: str) -> str:
    """
    Helper function to extract and clean the core part of a prior string.

    This function removes a given pattern (like 'beta_prior { ... }'),
    replaces whitespace sequences with a single underscore, and removes
    any resulting ':_' sequences.

    Args:
        prior_str (str): The input string to clean.
        pattern (str): The regex pattern for the wrapper to remove.

    Returns:
        str: The cleaned string.
    """
    # Remove the wrapper (e.g., "beta_prior { ... }")
    cleaned_str = re.sub(pattern, "", prior_str, flags=re.DOTALL)
    # Replace one or more whitespace characters with a single underscore
    cleaned_str = re.sub(r'\s+', '_', cleaned_str).strip('_')
    # Remove the pattern ':_' that can result from the previous step
    cleaned_str = cleaned_str.replace(":_", "")
    return cleaned_str

def create_output_path(
    num_components_prior: Optional[str] = None,
    p0_prior: Optional[str] = None,
    rho_prior: Optional[str] = None,
    sigma_prior: Optional[str] = None,
    graph_prior: Optional[str] = None,
) -> str:
    """
    Generate an output path given priors in ASCII protocol buffer format.

    Args:
        num_components_prior (str, optional): Prior on the number of components. Defaults to None.
        p0_prior (str, optional): Prior on P0. Defaults to None.
        rho_prior (str, optional): Prior on rho. Defaults to None.
        sigma_prior (str, optional): Prior on sigma. Defaults to None.
        graph_prior (str, optional): Prior on the graph. Defaults to None.

    Returns:
        str: A formatted output path string.
        
    Raises:
        ValueError: If any prior string has an unrecognized format.
    """
    out_path_parts = []

    # Prior on number of components
    if num_components_prior is not None:
        if "shifted_poisson_prior" in num_components_prior:
            pattern = r"shifted_poisson_prior\s*\{|\}"
            cleaned_str = _clean_prior_string(num_components_prior, pattern)
            out_path_parts.append(f"H-RJ{cleaned_str}")
        elif "fixed:" in num_components_prior:
            pattern = r"fixed:\s*"
            cleaned_str = _clean_prior_string(num_components_prior, pattern)
            out_path_parts.append(f"H-{cleaned_str}")
        else:
            raise ValueError("Wrong format for 'num_components_prior' input.")
    
    # Prior on P0
    if p0_prior is not None:
        if "p0_params" in p0_prior:
            pattern = r"p0_params\s*\{|\}"
            cleaned_str = _clean_prior_string(p0_prior, pattern)
            out_path_parts.append(f"P0-{cleaned_str}")
        else:
            raise ValueError("Wrong format for 'p0_prior' input.")

    # Prior on rho
    if rho_prior is not None:
        if "beta_prior" in rho_prior:
            pattern = r"beta_prior\s*\{|\}"
            cleaned_str = _clean_prior_string(rho_prior, pattern)
            out_path_parts.append(f"rho-{cleaned_str}")
        elif "fixed:" in rho_prior:
            pattern = r"fixed:\s*"
            cleaned_str = _clean_prior_string(rho_prior, pattern)
            out_path_parts.append(f"rho-{cleaned_str}")
        else:
            raise ValueError("Wrong format for 'rho_prior' input.")

    # Prior on sigma
    if sigma_prior is not None:
        if "inv_wishart_prior" in sigma_prior:
            pattern = r"inv_wishart_prior\s*\{|\}"
            cleaned_str = _clean_prior_string(sigma_prior, pattern)
            out_path_parts.append(f"sigma-{cleaned_str}")
        elif "inv_gamma_prior" in sigma_prior:
            pattern = r"inv_gamma_prior\s*\{|\}"
            cleaned_str = _clean_prior_string(sigma_prior, pattern)
            out_path_parts.append(f"sigma-{cleaned_str}")
        elif "fixed:" in sigma_prior:
            pattern = r"fixed:\s*"
            cleaned_str = _clean_prior_string(sigma_prior, pattern)
            out_path_parts.append(f"sigma-{cleaned_str}")
        else:
            raise ValueError("Wrong format for 'sigma_prior' input.")

    # Prior on the graph
    if graph_prior is not None:
        if "beta_prior" in graph_prior:
            pattern = r"beta_prior\s*\{|\}"
            cleaned_str = _clean_prior_string(graph_prior, pattern)
            out_path_parts.append(f"p-{cleaned_str}")
        elif "fixed:" in graph_prior:
            pattern = r"fixed:\s*"
            cleaned_str = _clean_prior_string(graph_prior, pattern)
            out_path_parts.append(f"p-{cleaned_str}")
        else:
            raise ValueError("Wrong format for 'graph_prior' input.")

    # Join all parts into a single path string
    return os.path.join(*out_path_parts)


def download_and_extract(url, dest_folder):
    """
    Downloads a zip file from a URL and extracts it to a specified folder.

    Args:
        url (str): URL of the zip file to download
        dest_folder (str): Name of the folder to extract files to
    """
    # Create destination folder if it doesn't exist
    dest_path = Path(dest_folder)
    dest_path.mkdir(exist_ok=True)
    
    # Download the file
    print(f"Downloading from {url}")
    response = requests.get(url)
    zip_path = dest_path / "temp.zip"
    
    # Save the zip file
    with open(zip_path, 'wb') as f:
        f.write(response.content)
    
    # Extract the contents
    print("Extracting files")
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(dest_path)
    
    # Remove the temporary zip file
    zip_path.unlink()
    print(f"Files extracted to {dest_path}")

# Specify download URLs
urls = ["https://www2.census.gov/programs-surveys/acs/experimental/2020/data/pums/1-Year/csv_pca.zip",
        "https://usa.ipums.org/usa/resources/volii/shapefiles/ipums_puma_2010_tl20.zip"]

# Define the action function for the download_data task
def download_data(_ : Task) -> None:
    for url in urls:
        download_and_extract(url, dest_folder="raw")

# Define download_data task
download_data_targets = [
    "raw/psam_p06.csv",
    "raw/ipums_puma_2010_tl20.shp"
]
create_task("download_data", action=download_data, targets=download_data_targets)

# Define generate_data task
num_datasets = 10
subsample_size = 100
generate_data_action = ["Rscript", "src/generate_data.R",
                        "--num-datasets", num_datasets,
                        "--subsample-size", subsample_size]
generate_data_targets = ["input/counties-pumas/counties-pumas.shp"] + \
    [f"input/data_{i:03d}.dat" for i in range(1, num_datasets+1)] + \
    ["input/adj_matrix.dat"]
create_task("generate_data", action = generate_data_action, dependencies = download_data_targets, targets = generate_data_targets)

# Define run_full_dataset priors to explore
num_components_priors = ["shifted_poisson_prior { rate: 1.0 }"]
rho_priors = ["fixed: 0.95"]
sigma_priors = ["inv_gamma_prior { alpha: 6 beta: 4 }"]
graph_priors = ["beta_prior { a: 2 b: 93 }"]

# Define run_full_dataset tasks
with create_group("run_full_dataset") as parallel_runs_full_dataset_group:
    for num_components_prior in num_components_priors:
        for rho_prior in rho_priors:
            for sigma_prior in sigma_priors:
                for graph_prior in graph_priors:
                    output_path = create_output_path(num_components_prior, None, rho_prior, sigma_prior, graph_prior)
                    run_full_dataset_action = ["Rscript", "src/run_sampler_new.R",
                                               "--num-components-prior", num_components_prior,
                                               "--rho-prior", rho_prior,
                                               "--sigma-prior", sigma_prior,
                                               "--graph-prior", graph_prior,
                                               "--output-file", f"output-new/{output_path}/full_dataset_chain.dat",
                                               "input/full_dataset.dat"]
                    # run_full_dataset_targets = [f"output/{output_path}/full_dataset_chain.dat"]
                    create_task(f"_run_full_dataset-{output_path}", action=run_full_dataset_action)#,
                                # dependencies=["input/full_dataset.dat"],
                                # targets=run_full_dataset_targets)

with create_group("generate_plots") as parallel_generate_plots_group:
    for num_components_prior in num_components_priors:
        for rho_prior in rho_priors:
            for sigma_prior in sigma_priors:
                for graph_prior in graph_priors:
                    output_path = create_output_path(num_components_prior, None, rho_prior, sigma_prior, graph_prior)
                    generate_plot_action = ["Rscript", "src/generate_plot.R",
                                            "--data-file", "input/full_dataset.dat",
                                            "--sim-file", f"output-new/{output_path}/full_dataset_chain.dat",
                                            "--output-dir", f"plots-new/{output_path}/full_dataset"]
                    create_task(f"_generate_plots-{output_path}", action=generate_plot_action)

# Define generate_plot task
generate_plot_action = ["Rscript", "src/generate_plot.R",
                        "--data-file", "input/full_dataset.dat",
                        "--sim-file", "dump/output/HRJ/rho0.95/alpha6_beta4/a2_b93/full_dataset_chain.dat",
                        "--output-dir", "plots-def/HRJ/rho0.95/alpha6_beta4/a2_b93/full_dataset"]
create_task("generate_plot", action=generate_plot_action)#,
            # dependencies=["input/full_dataset.dat", "input/counties-pumas/counties-pumas.shp"])


# Define generate_explain_boundaries_shapefiles task
generate_explain_boundaries_shapefiles_action = ["Rscript", "src/generate_explain_boundaries_shapefiles.R",
                                                 "--output-dir", "input/explain_boundaries"]
generate_explain_boundaries_shapefiles_targets = ["input/explain_boundaries/crimes_sf.dat", "input/explain_boundaries/health_insurance_sf.dat"]
create_task("generate_explain_boundaries_shapefiles", action=generate_explain_boundaries_shapefiles_action,
            targets=generate_explain_boundaries_shapefiles_targets)