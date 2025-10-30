# # ---- BD SCENARIO 1 - GENERATE TABLES ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_tables", hide.opts = TRUE,
                         description = "Generate summary tables for the simulated scenario")
opt_parser <- add_argument(opt_parser, arg = "--summary-path", type = "character",
                           help = "Relative path to the summary folder")
opt_parser <- add_argument(opt_parser, arg = "--num-components-values", type = "character",
                           help = "Comma-separated values for the number of components considered")
opt_parser <- add_argument(opt_parser, arg = "--poisson-rate-values", type = "character", default = NULL,
                           help = "Comma-separated values for the Poisson rate parameter considered (for RJMCMC)")
opt_parser <- add_argument(opt_parser, arg = "--rho-values", type = "character",
                           help = "Comma-separated values for the rho parameter considered")
opt_parser <- add_argument(opt_parser, arg = "--output-dir", type = "character",
                           help = "Relative path to the output folder")
extra_args <- parse_args(opt_parser)


# Preliminary checks ------------------------------------------------------

# Stop if summary folder does not exist
if(is.na(extra_args$summary_path)) {
  print(opt_parser)
  stop("Input parameter '--summary-path' not specified")
}

# Stop if number of components values not specified
if(is.na(extra_args$num_components_values)) {
  print(opt_parser)
  stop("Input parameter '--num-components-values' not specified")
}

# Stop if rho values not specified
if(is.na(extra_args$rho_values)) {
  print(opt_parser)
  stop("Input parameter '--rho-values' not specified")
}

# Stop if output directory not specified
if(is.na(extra_args$output_dir)) {
  print(opt_parser)
  stop("Input parameter '--output-dir' not specified")
}

# Set working directory relative to the script location
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  # Running in RStudio
  setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))
} else {
  # Running from command line
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  setwd(dirname(dirname(script.name)))
}
cat("Setting working directory to: ", getwd(), "\n")

# Check if summary path is specified
summary_path <- file.path(getwd(), extra_args$summary_path)
if(!dir.exists(summary_path)){
  stop(sprintf("'%s' does not exists", summary_path))
}
cat(sprintf("Summary Folder: %s\n", summary_path)) # Log

# Create output directory if does not exist
output_dir <- file.path(getwd(), extra_args$output_dir)
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
cat(sprintf("Output Folder: %s\n", output_dir)) # Log

# Parse H and rho values
H <- unlist(strsplit(extra_args$num_components_values, split = ","))
if("RJ" %in% H && is.null(extra_args$poisson_rate_values)){
  stop("Please provide values for --poisson-rate-values when using RJMCMC")
}
if("RJ" %in% H) {
  poisson_rate_values <- unlist(strsplit(extra_args$poisson_rate_values, split = ","))
  H <- c(setdiff(H, "RJ"), sapply(poisson_rate_values, function(pr){sprintf("RJ-poisson%s", pr)}))
}
rho <- unlist(strsplit(extra_args$rho_values, split = ","))


# Main code ---------------------------------------------------------------

# Create table buffers
precisions <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
sensitivities <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
specificities <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
meanL1 <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
meanWAIC <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
meanAUC <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))

# Generate tables
for (i in 1:length(H)) {
  for (j in 1:length(rho)) {
    # Parse confusion_matrix csv file
    file_name <- file.path(summary_path, sprintf("confusion_matrices-H%s-rho%s.csv", H[i], rho[j]))
    if(!file.exists(file_name)){
      stop(sprintf("File '%s' not found", file_name))
    } else {
      cat(sprintf("Parsing file: %s\n", file_name)) # Log
      CM <- read.csv(file_name)
    }
    # Compute precision
    prec_vec <- sapply(1:nrow(CM), function(e){CM[e,"TP"]/(CM[e,"TP"]+CM[e,"FP"])}); prec_vec[which(is.nan(prec_vec))] <- 0
    precisions[i,j] <- sprintf("%g (%g)", round(mean(prec_vec),3), round(sd(prec_vec), 3))
    # Compute sensitivity
    sens_vec <- sapply(1:nrow(CM), function(e){CM[e,"TP"]/(CM[e,"TP"]+CM[e,"FN"])})
    sensitivities[i,j] <- sprintf("%g (%g)", round(mean(sens_vec),3), round(sd(sens_vec), 3))
    # Compute specificity
    spec_vec <- sapply(1:nrow(CM), function(e){CM[e,"TN"]/(CM[e,"TN"]+CM[e,"FP"])})
    specificities[i,j] <- sprintf("%g (%g)", round(mean(spec_vec),3), round(sd(spec_vec),3))

    # Parse mean L1 csv file
    file_name <- file.path(summary_path, sprintf("mean_L1_distances-H%s-rho%s.csv", H[i], rho[j]))
    if(!file.exists(file_name)){
      stop(sprintf("File '%s' not found", file_name))
    } else {
      cat(sprintf("Parsing file: %s\n", file_name)) # Log
      L1_df <- read.csv(file_name)
    }
    # Compute overall mean L1 distance
    meanL1[i,j] <- sprintf("%g (%g)", round(mean(L1_df[,"MeanL1"]),3), round(sd(L1_df[,"MeanL1"]), 3))

    # Parse WAIC csv file
    file_name <- file.path(summary_path, sprintf("WAIC-H%s-rho%s.csv", H[i], rho[j]))
    if(!file.exists(file_name)){
      stop(sprintf("File '%s' not found", file_name))
    } else {
      cat(sprintf("Parsing file: %s\n", file_name)) # Log
      WAIC_df <- read.csv(file_name)
    }
    # Compute mean WAIC
    meanWAIC[i,j] <- sprintf("%g (%g)", round(mean(WAIC_df[,"WAIC"]),1), round(sd(WAIC_df[,"WAIC"]), 1))

    # Parse ROC curves csv file
    file_name <- file.path(summary_path, sprintf("ROC_curves-H%s-rho%s.csv", H[i], rho[j]))
    if(!file.exists(file_name)){
      stop(sprintf("File '%s' not found", file_name))
    } else {
      cat(sprintf("Parsing file: %s\n", file_name)) # Log
      ROC_df <- read.csv(file_name)
    }
    # Compute mean AUC
    meanAUC[i,j] <- sprintf("%g (%g)", round(mean(ROC_df[,"auc"]),3), round(sd(ROC_df[,"auc"]), 3))
  }
}

# Save tables
write.csv(precisions, file = file.path(output_dir, "precision_table.csv"), row.names = TRUE)
write.csv(sensitivities, file = file.path(output_dir, "sensitivity_table.csv"), row.names = TRUE)
write.csv(specificities, file = file.path(output_dir, "specificity_table.csv"), row.names = TRUE)
write.csv(meanL1, file = file.path(output_dir, "mean_L1_table.csv"), row.names = TRUE)
write.csv(meanWAIC, file = file.path(output_dir, "mean_WAIC_table.csv"), row.names = TRUE)
write.csv(meanAUC, file = file.path(output_dir, "mean_AUC_table.csv"), row.names = TRUE)

# # Table of Elapsed Times --------------------------------------------------

# parseTimeDF <- function(h){
#   file_name <- file.path(getwd(), "summary", alphabeta, sprintf("elapsed_times-H_%s.txt", h))
#   out <- data.frame(read.table(file_name, header = F))
#   names(out) <- h
#   return(out)
# }
# TimesDF <- do.call(cbind, lapply(H, parseTimeDF))

