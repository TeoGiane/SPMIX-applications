from cook import create_task, Task
# from cook.contexts import create_group
import requests
import zipfile
from pathlib import Path

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
create_task("download_data", action=download_data,
            targets=download_data_targets)

# Define generate_data task
num_datasets = 10
subsample_size = 100
generate_data_action = ["Rscript", "src/generate_data.R",
                        "--num-datasets", num_datasets,
                        "--subsample-size", subsample_size]
create_task("generate_data", action=generate_data_action)#,
            # dependencies=download_data_targets,
            # targets=["input/counties-pumas/counties-pumas.shp"] +
            #         [f"input/data_{i:03d}.dat" for i in range(1, num_datasets+1)] +
            #         ["input/adj_matrix.dat"])

# Define run_full_dataset task
run_full_dataset_action = ["Rscript", "src/run_sampler.R",
                           "--num-components", "RJ",
                           "--rho", 0.95,
                           "--output-file", "output/H_RJ/rho_0.95/strong_p/full_dataset_chain.dat",
                           "input/full_dataset.dat"]
run_full_dataset_targets = ["output/H_RJ/rho_0.95/strong_p/full_dataset_chain.dat"]
create_task("run_full_dataset", action=run_full_dataset_action)#,
            # dependencies=["input/full_dataset.dat"],
            # targets=run_full_dataset_targets)

# Define generate_plot task
generate_plot_action = ["Rscript", "src/generate_plot.R",
                        "--data-file", "input/full_dataset.dat",
                        "--sim-file", "output/H_RJ/rho_0.95/strong_p/full_dataset_chain.dat",
                        "--output-dir", "plots/H_RJ/rho_0.95/strong_p/full_dataset"]
create_task("generate_plot", action=generate_plot_action)#,
            # dependencies=["input/full_dataset.dat", "input/counties-pumas/counties-pumas.shp"])