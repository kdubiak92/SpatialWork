FROM condaforge/miniforge3

# Create environment
RUN conda create -n seurat_env -c conda-forge -c bioconda \
    r-base=4.3 \
    r-seurat \
    r-hdf5r \
    r-remotes \
    r-optparse \
    r-future \
    r-future.apply \
    r-rmarkdown \
    pandoc \
    libxml2 \
    libcurl \
    openssl \
    git \
    && conda clean -afy

# Activate environment and set as default
ENV PATH=/opt/conda/envs/seurat_env/bin:$PATH
ENV CONDA_DEFAULT_ENV=seurat_env

# Install additional R packages
RUN Rscript -e "install.packages('harmony', repos='https://cloud.r-project.org/')"
# Install SeuratDisk from GitHub (in case conda version lags)
RUN Rscript -e "remotes::install_github('mojaveazure/seurat-disk')"