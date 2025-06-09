#!/usr/bin/env Rscript

# Example usage:
# Rscript render-sample-anlys.R \
#   -x sample-anlys.Rmd \
#   -i /path/to/spaceranger/outs \
#   -s SAMPLE-ID \
#   -b .008um \
#   -o output_directory
#
# Optional parameters with defaults:
#   --detected_genes_low_threshold=0.01 \
#   --detected_genes_high_threshold=0.99 \
#   --umi_count_low_threshold=0.01 \
#   --umi_count_high_threshold=0.99 \
#   --ribosomal_contents_threshold=25 \
#   --mitochondrial_contents_threshold=50 \
#   --norm_method=LogNormalize \
#   --varFeatures=2000 \
#   --remove_mito_genes=FALSE \
#   --remove_ribo_genes=FALSE \
#   --npcs=20

library(optparse)
library(rmarkdown)

# Define command-line options
option_list <- list(
  make_option(c("-x", "--rmd_file"),
    type = "character", default = "sample-anlys.Rmd",
    help = "R Markdown File [default: %default]", metavar = "FILE"
  ),
  make_option(c("-i", "--input_dir"),
    type = "character", default = "data/space_ranger_out",
    help = "Input directory for 10X spatial data [default: %default]", metavar = "DIRECTORY"
  ),
  make_option(c("-s", "--sample_id"),
    type = "character", default = "sample",
    help = "Sample ID [default: %default]", metavar = "character"
  ),
  make_option(c("-b", "--target_bin"),
    type = "character", default = ".008um",
    help = "Target bin to use; choose among '.002um', '.008um', or '.016um' [default: %default]",
    metavar = "BIN"
  ),
  make_option(c("-o", "--output_directory"),
    type = "character", default = "output_directory",
    help = "Output directory to save HTML and Seurat objects [default: %default]",
    metavar = "DIRECTORY"
  ),
  # QC threshold options
  make_option("--detected_genes_low_threshold",
    type = "numeric", default = 0.01,
    help = "Lower quantile threshold for detected genes [default: %default]"
  ),
  make_option("--detected_genes_high_threshold",
    type = "numeric", default = 0.99,
    help = "Upper quantile threshold for detected genes [default: %default]"
  ),
  make_option("--umi_count_low_threshold",
    type = "numeric", default = 0.01,
    help = "Lower quantile threshold for UMI counts [default: %default]"
  ),
  make_option("--umi_count_high_threshold",
    type = "numeric", default = 0.99,
    help = "Upper quantile threshold for UMI counts [default: %default]"
  ),
  make_option("--ribosomal_contents_threshold",
    type = "numeric", default = 25,
    help = "Maximum percentage of ribosomal content [default: %default]"
  ),
  make_option("--mitochondrial_contents_threshold",
    type = "numeric", default = 50,
    help = "Maximum percentage of mitochondrial content [default: %default]"
  ),
  # Analysis options
  make_option("--norm_method",
    type = "character", default = "LogNormalize",
    help = "Normalization method [default: %default]"
  ),
  make_option("--varFeatures",
    type = "integer", default = 2000,
    help = "Number of variable features to select [default: %default]"
  ),
  make_option("--remove_mito_genes",
    type = "logical", default = FALSE,
    help = "Remove mitochondrial genes [default: %default]"
  ),
  make_option("--remove_ribo_genes",
    type = "logical", default = FALSE,
    help = "Remove ribosomal genes [default: %default]"
  ),
  make_option("--npcs",
    type = "integer", default = 20,
    help = "Number of principal components to use [default: %default]"
  )
)

# Parse the command-line options
opt <- parse_args(OptionParser(option_list = option_list))

# Keep existing target_bin validation and mapping logic
if (!(opt$target_bin %in% c(".002um", ".008um", ".016um"))) {
  stop(
    "Invalid target_bin value: ", opt$target_bin,
    ". Allowed values are: .002um, .008um, .016um"
  )
}

# Map target_bin to the correct format
bin_mapping <- list(
  ".002um" = "Spatial.002um",
  ".008um" = "Spatial.008um",
  ".016um" = "Spatial.016um"
)

# Update target_bin to the correct format
opt$target_bin <- bin_mapping[[opt$target_bin]]

# Check if the required options are provided
if (is.null(opt$rmd_file) || is.null(opt$input_dir) ||
  is.null(opt$sample_id) || is.null(opt$output_directory)) {
  stop("Please provide the Rmd file, input directory, sample ID, and output directory.")
}

# Ensure the output directory exists
dir.create(opt$output_directory, showWarnings = FALSE, recursive = TRUE)

# Validate input directory
if (!dir.exists(opt$input_dir)) {
  stop("Input directory does not exist: ", opt$input_dir)
}

# Parameters to pass to the Rmd file
params <- list(
  input_dir = opt$input_dir,
  sample_id = opt$sample_id,
  target_bin = opt$target_bin,
  output_directory = opt$output_directory,
  # QC parameters
  detected_genes_low_threshold = opt$detected_genes_low_threshold,
  detected_genes_high_threshold = opt$detected_genes_high_threshold,
  umi_count_low_threshold = opt$umi_count_low_threshold,
  umi_count_high_threshold = opt$umi_count_high_threshold,
  ribosomal_contents_threshold = opt$ribosomal_contents_threshold,
  mitochondrial_contents_threshold = opt$mitochondrial_contents_threshold,
  # Analysis parameters
  norm_method = opt$norm_method,
  varFeatures = opt$varFeatures,
  remove_mito_genes = opt$remove_mito_genes,
  remove_ribo_genes = opt$remove_ribo_genes,
  npcs = opt$npcs
)

# Define output file name
default_output_file <- file.path(
  opt$output_directory,
  paste0(opt$sample_id, "_preprocessing_report.html")
)

# Render the R Markdown file
tryCatch(
  {
    rmarkdown::render(
      input = opt$rmd_file,
      output_format = "html_document",
      params = params,
      output_file = default_output_file,
      envir = new.env(parent = globalenv())
    )
    cat("Rendering completed. Output saved to:", default_output_file, "\n")
  },
  error = function(e) {
    cat("Error during rendering:", e$message, "\n")
    quit(status = 1)
  }
)