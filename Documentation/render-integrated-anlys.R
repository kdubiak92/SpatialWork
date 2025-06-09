#!/usr/bin/env Rscript

# Example Command:
# Rscript render-integrated-anlys.R \
# -r "output_directory/SD-A1_processed.rds,output_directory/SD-D1_processed.rds" \
# -t .008um \
# -n LogNormalize \
# -p 20 \
# -i integrated-anlys-visiumhd.Rmd \
# -o output_integration

# Load necessary libraries
library(optparse)
library(rmarkdown)

# Define command-line options
option_list <- list(
  make_option(c("-r", "--rds_paths"),
    type = "character", default = NULL,
    help = "Comma-separated list of RDS files paths containing processed Seurat objects",
    metavar = "RDS_FILES"
  ),
  make_option(c("-t", "--target_assay"),
    type = "character", default = ".008um",
    help = "Target bin to use; choose among '.002um', '.008um', or '.016um'",
    metavar = "BIN"
  ),
  make_option(c("-n", "--normalization_method"),
    type = "character", default = "LogNormalize",
    help = "Method used for data normalization", metavar = "METHOD"
  ),
  make_option(c("-p", "--npcs"),
    type = "integer", default = 20,
    help = "Number of principal components to use for integration", metavar = "NUMBER"
  ),
  make_option(c("-i", "--input"),
    type = "character", default = NULL,
    help = "Path to the input R Markdown (.Rmd) file", metavar = "FILE"
  ),
  make_option(c("-o", "--output_directory"),
    type = "character", default = NULL,
    help = "Output directory for results and HTML report",
    metavar = "DIRECTORY"
  )
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate and map target_assay
if (!(opt$target_assay %in% c(".002um", ".008um", ".016um"))) {
  stop(
    "Invalid target_assay value: ", opt$target_assay,
    ". Allowed values are: .002um, .008um, .016um"
  )
}

# Map target_assay to the correct format
bin_mapping <- list(
  ".002um" = "Spatial.002um",
  ".008um" = "Spatial.008um",
  ".016um" = "Spatial.016um"
)

# Update target_assay to the correct format
opt$target_assay <- bin_mapping[[opt$target_assay]]

# Check for required arguments
if (is.null(opt$rds_paths)) {
  stop("The --rds_paths argument is required.")
}
if (is.null(opt$input)) {
  stop("The --input argument is required.")
}
if (is.null(opt$output_directory)) {
  stop("The --output_directory argument is required.")
}

# Process rds_paths (trim whitespace from file paths)
rds_files <- trimws(unlist(strsplit(opt$rds_paths, ",")))
for (rds_file in rds_files) {
  if (!file.exists(rds_file)) {
    stop(paste("RDS file not found:", rds_file))
  }
}

# Create output directory if it doesn't exist
dir.create(opt$output_directory, showWarnings = FALSE, recursive = TRUE)

# Construct output HTML filename
output_html <- file.path(opt$output_directory, "integration_analysis.html")

cat("RDS files to process:\n")
cat(paste(rds_files, collapse = "\n"), "\n")

# Render the R Markdown file
cat("Rendering R Markdown file...\n")
rmarkdown::render(
  input = opt$input,
  output_file = output_html,
  params = list(
    rds_paths = rds_files,
    target_assay = opt$target_assay,
    output_directory = opt$output_directory,
    normalization_method = opt$normalization_method,
    npcs = opt$npcs
  )
)

cat("Analysis complete. Results saved to:", opt$output_directory, "\n")
cat("HTML report:", output_html, "\n")
