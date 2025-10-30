# Reference image
FROM python:3.12

# Update packages
RUN apt-get update && apt-get upgrade -y

# Install required libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    jags \
    protobuf-compiler \
    libprotobuf-dev \
    libprotoc-dev \
    libgdal-dev \
    libudunits2-dev \
    libgsl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Install R and dependencies
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /workdir

# Install R dependencies
COPY requirements.R .
RUN Rscript requirements.R

# Install SPMIX R package from GitHub
RUN Rscript -e 'devtools::install_github("TeoGiane/SPMIX", ref = "dev")'

# Install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
