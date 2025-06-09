# Use satijalab/seurat:5.0.0 as the base image
# This image should include R, Seurat v5.0.0, and many common dependencies,
# including Seurat, SeuratDisk, and glue (as you've confirmed).
FROM satijalab/seurat:5.0.0

# --- System Dependencies ---
# Install system libraries that might be needed for Python.
# We're relying on the base image for most R-related system dependencies now.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    # Python is essential for your workflow
    python3-dev \
    # Clean up APT cache to reduce image size
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# --- R Package Installation ---
# No R packages to install, as confirmed by interactive testing.

# --- Python Package Installation ---
# Install Python packages required for AnnData/h5ad format.
RUN pip3 install --no-cache-dir \
    anndata \
    h5py \
    pandas \
    numpy \
    scipy \
    scanpy