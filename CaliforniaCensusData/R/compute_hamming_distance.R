# Required libraries
library("SPMIX")
library("sf")
library("spdep")

# FUNCTION - get boundary graph from serialized chain file
get_Gb <- function(chain_file, adj_matrix) {
  # Load serialized chain
  load(chain_file)
  # Compute admissible edges
  Eadj <- which(adj_matrix == 1, arr.ind = TRUE)
  # Deserialize chain and compute plinks
  chain <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  G_chain <- lapply(chain, function(x) matrix(x$G$data,x$G$rows,x$G$cols))
  plinks <- Reduce('+', G_chain)/length(G_chain)
  # Compute boundary graph (whole matrix)
  Gb <- matrix(0L, nrow(plinks), ncol(plinks))
  Gb[Eadj] <- ifelse(plinks[Eadj] < 0.5, 1L, 0L)
  # Return upper triangular matrix
  return(Gb[upper.tri(Gb)])
}

# FUNCTION - Compute standardized hamming distance
compute_std_hamming_distance <- function(Gb_lhs, Gb_rhs) {
  # Check dimension match
  if(length(Gb_lhs) != length(Gb_rhs)){
    stop("Dimensions does not match")
  }
  # Compute standarized hamming distance
  return(sum(Gb_lhs != Gb_rhs) / length(Gb_lhs))
}


# Main code ---------------------------------------------------------------

# Import shapefile
sf_counties <- read_sf("data/counties-pumas/counties-pumas.shp")

# Compute adjacency list and matrix
adj_list <- poly2nb(sf_counties, queen = FALSE)
W <- nb2mat(adj_list, style = "B")

# List all available chain files
chain_files <- list.files("output/H_RJ/rho_0.95", recursive = F, full.names = T)

# Compute std hamming distances between all pairs
std_hamming_distances <- vector(length = choose(length(chain_files),2)); k <- 1
for (i in 1:(length(chain_files)-1)) {
  Gb_lhs <- get_Gb(chain_files[i],W)
  for (j in (i+1):length(chain_files)) {
    std_hamming_distances[k] <- compute_std_hamming_distance(Gb_lhs, get_Gb(chain_files[j],W))
    k <- k+1
  }
}

# Return summary table
summary(std_hamming_distances)
