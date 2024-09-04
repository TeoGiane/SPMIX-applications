#! /bin/bash

# Define useful variables
NUM_DATASETS=100
SUBSAMPLE_SIZE=100

# Clean up data/ folder if already present
echo "Clean up"
rm -rf data/

# Download raw.zip file
echo "Getting raw data"
wget -nv 'https://api.onedrive.com/v1.0/shares/u!aHR0cHM6Ly8xZHJ2Lm1zL3UvcyFBbEhCVllHSnJYREJqZDhzR0Q4MmRpbTZERl9KTEE_ZT1iQ1pZTGw/root/content' -O raw.zip

# Unzip folder
echo "Unzip raw.zip file into raw/"
unzip -q raw.zip && rm raw.zip

# Generate the data folder via R script
echo "Run generate_data.R script"
Rscript --vanilla ./R/generate_data.R -n $NUM_DATASETS -s $SUBSAMPLE_SIZE

# Remove the raw/ folder
echo "Remove raw/ folder"
rm -rf raw/
