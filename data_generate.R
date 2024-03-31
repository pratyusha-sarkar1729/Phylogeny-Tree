rm(list = ls())

# Define the base sequence
base_sequence <- "GCAT"

# Define a mapping from nucleotides to numeric values
nucleotide_mapping <- c(A = 1, T = 2, C = 3, G = 4)

# Convert the base sequence to numeric representation
numeric_sequence <- sapply(strsplit(base_sequence, "")[[1]], function(nucleotide) nucleotide_mapping[nucleotide])

# Define the number of sequences to generate
num_sequences <- 3

# Define the Jukes-Cantor model parameter
mu <- 5

# Define the alphabet of nucleotides
nucleotides <- c("A", "T", "C", "G")

# Function to calculate Jukes-Cantor model probabilities with custom transition probability matrix
calculate_transition_probs <- function(branch_length, mu) {
  # Initialize an empty list to store transition probability matrices
  transition_matrices <- list()
  
  # Loop over each branch length*
  for (bl in branch_length) {
    # Compute the probabilities for same and different nucleotides
    p_same <- 1 - 3/4 * exp(-4/3 * mu * bl)
    p_diff <- 1/4 * exp(-4/3 * mu * bl)
    
    # Create the transition probability matrix with p_diff
    transition_matrix <- matrix(rep(p_diff, 16), nrow = 4, byrow = TRUE)
    
    # Change the diagonal elements to p_same
    diag(transition_matrix) <- rep(p_same, 4)
    
    # Append the transition matrix to the list
    transition_matrices[[length(transition_matrices) + 1]] <- transition_matrix
  }
  
  return(transition_matrices)
}

# Function to generate a new sequence using Jukes-Cantor model
generate_new_sequence <- function(base_sequence, mu, branch_length) {
  # Calculate transition probabilities
  transition_probs <- calculate_transition_probs(branch_length, mu)
  
  # Initialize new sequence
  new_sequence <- ""
  
  # Convert the base sequence to numeric representation
  numeric_sequence <- sapply(strsplit(base_sequence, "")[[1]], function(nucleotide) nucleotide_mapping[nucleotide])
  
  # Loop through each base in the base sequence
  for (base in numeric_sequence) {
    # Determine the transition probabilities for this base
    prob_vector <- transition_probs[[length(new_sequence)]][base, ]
    
    # Generate new base according to Jukes-Cantor model probabilities
    new_base <- sample(nucleotides, size = 1, prob = prob_vector)
    
    # Append new base to the new sequence
    new_sequence <- paste0(new_sequence, new_base)
  }
  
  return(new_sequence)
}



library(ape)


# Set the number of tips (leaves) in the tree
num_tips <- 6

# Create a star tree with the specified number of tips
star_tree <- rcoal(num_tips,br = sort(rexp(num_tips-1, rate = 10),decreasing = FALSE))

# Plot the ladder-like tree
tree_structure <- plot(star_tree, cex = 0.6)

branch_length_new = sort(star_tree$edge.length,decreasing = TRUE)
branch1 <- generate_new_sequence(base_sequence, mu, branch_length_new[1])
branch2 <- generate_new_sequence(base_sequence, mu, branch_length_new[2])
branch3 <- generate_new_sequence(branch2, mu, branch_length_new[3])
branch4 <- generate_new_sequence(branch2, mu, branch_length_new[4])
branch5 <- generate_new_sequence(branch4, mu, branch_length_new[5])
branch6 <- generate_new_sequence(branch4, mu, branch_length_new[8])
branch7 <- generate_new_sequence(branch5, mu, branch_length_new[9])
branch8 <- generate_new_sequence(branch5, mu, branch_length_new[10])
branch9 <- generate_new_sequence(branch6, mu, branch_length_new[6])
branch10 <- generate_new_sequence(branch6, mu, branch_length_new[7])



star_tree$tip.label = c(branch3,branch9,branch10,branch7,branch8,branch1)

tree_structure_final <- plot(star_tree, cex = 0.6)

# Save the final tree plot as an image file
# Save the tree_structure_plot object
save(tree_structure_final, file = "tree_structure_plot.RData")
