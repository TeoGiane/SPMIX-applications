# Table for CM ------------------------------------------------------------

# Define values for H and rho
H <- c(2,4,6,8,10,'RJ')
rho <- c(0,0.5,0.9,0.95,0.99)

# Generate CM tables
precisions <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
sensitivities <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
specificities <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
for (i in 1:length(H)) {
  for (j in 1:length(rho)) {
    # Parse csv file
    file_name <- file.path(getwd(), "summary", sprintf("CM-H_%s-rho_%s.csv", H[i], rho[j]))
    CM <- read.csv(file_name)
    # Compute precision
    prec_vec <- sapply(1:nrow(CM), function(e){CM[e,"TP"]/(CM[e,"TP"]+CM[e,"FP"])}); prec_vec[which(is.nan(prec_vec))] <- 0
    precisions[i,j] <- sprintf("%g (%g)", round(mean(prec_vec),3), round(sd(prec_vec), 3))
    # Compute sensitivity
    sens_vec <- sapply(1:nrow(CM), function(e){CM[e,"TP"]/(CM[e,"TP"]+CM[e,"FN"])})
    sensitivities[i,j] <- sprintf("%g (%g)", round(mean(sens_vec),3), round(sd(sens_vec), 3))
    # Compute specificity
    spec_vec <- sapply(1:nrow(CM), function(e){CM[e,"TN"]/(CM[e,"TN"]+CM[e,"FP"])})
    specificities[i,j] <- sprintf("%g (%g)", round(mean(spec_vec),3), round(sd(spec_vec),3))
  }
}

# Table for meanL1 --------------------------------------------------------

# Define values for H and rho
H <- c(2,4,6,8,10,'RJ')
rho <- c(0,0.5,0.9,0.95,0.99)

# Generate CM tables
meanL1 <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
for (i in 1:length(H)) {
  for (j in 1:length(rho)) {
    # Parse csv file
    file_name <- file.path(getwd(), "summary", sprintf("meanL1-H_%s-rho_%s.csv", H[i], rho[j]))
    L1_df <- read.csv(file_name)
    # Compute overall mean L1 distance
    meanL1[i,j] <- sprintf("%g (%g)", round(mean(L1_df[,"MeanL1"]),3), round(sd(L1_df[,"MeanL1"]), 3))
  }
}


# Table for WAIC ----------------------------------------------------------

# Define values for H and rho
H <- c(2,4,6,8,10,'RJ')
rho <- c(0,0.5,0.9,0.95,0.99)

# Generate CM tables
meanWAIC <- matrix(NA, nrow = length(H), ncol = length(rho), dimnames = list(H,rho))
for (i in 1:length(H)) {
  for (j in 1:length(rho)) {
    # Parse csv file
    file_name <- file.path(getwd(), "summary", sprintf("meanWAIC-H_%s-rho_%s.csv", H[i], rho[j]))
    WAIC_df <- read.csv(file_name)
    # Compute overall mean L1 distance
    meanWAIC[i,j] <- sprintf("%g (%g)", round(mean(L1_df[,"MeanWAIC"]),3), round(sd(L1_df[,"MeanWAIC"]), 3))
  }
}
