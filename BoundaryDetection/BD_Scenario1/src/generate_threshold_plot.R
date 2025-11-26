# # ---- BD SCENARIO 1 - GENERATE THRESHOLD PLOTS ---- # #

# Command line input options via argparser
suppressMessages(library("argparser"))
opt_parser <- arg_parser(name = "generate_threshold_plot", hide.opts = TRUE,
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

suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))

# Generate ROC_df
ROC_df <- list()
for (i in 1:length(H)) {
  for (j in 1:length(rho)) {
    file_name <- file.path(summary_path, sprintf("ROC_curves-H%s-rho%s.csv", H[i], rho[j]))
    if(!file.exists(file_name)){
      stop(sprintf("File '%s' not found", file_name))
    } else {
      cat(sprintf("Parsing file: %s\n", file_name)) # Log
      tmp_df <- read.csv(file_name) %>%
        group_by(threshold) %>% 
        mutate("mean_fnr" = mean(fn/(tp+fn)),
               "mean_fpr" = mean(fp/(fp+tn)))
      tmp_df$model <- rep(sprintf("H = %s, rho = %s",H[i],rho[j]), nrow(tmp_df))
      ROC_df <- append(ROC_df, list(tmp_df))
    }
  }
}
ROC_df <- do.call(rbind, ROC_df)
ROC_df$model <- factor(ROC_df$model, levels = unique(ROC_df$model))

# PLOT - FPR and FNR for all models in facet_wrap
plt_threshold_facet <- ggplot(data = ROC_df) +
  geom_line(aes(x = threshold, y = fp/(fp+tn), group=dataset_id), color='steelblue1', alpha = 0.3, linetype = 2, linewidth = 0.5) +
  geom_line(aes(x = threshold, y = mean_fpr), color='steelblue4', linewidth=1) +
  geom_line(aes(x = threshold, y = fn/(tp+fn), group=dataset_id), color='darkorange1', alpha = 0.3, linetype = 2, linewidth = 0.5) +
  geom_line(aes(x = threshold, y = mean_fnr), color='darkorange4', linewidth=1) +
  geom_vline(xintercept = 0.5, color='darkred', linewidth=1, linetype = 4) +
  facet_wrap(~model, ncol = 5)
pdf(file.path(output_dir, "plt_threshold_facet.pdf"), height = 16, width = 16); print(plt_threshold_facet); dev.off()


# # ---- END OF SCRIPT ---- # #
