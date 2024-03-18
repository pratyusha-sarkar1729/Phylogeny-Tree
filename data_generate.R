rm(list = ls())

# Define the base sequence
base_sequence <- "GCAT"

# Define a mapping from nucleotides to numeric values
nucleotide_mapping <- c(A = 1, T = 2, C = 3, G = 4)

# Convert the base sequence to numeric representation
numeric_sequence <- sapply(strsplit(base_sequence, "")[[1]], function(nucleotide) nucleotide_mapping[nucleotide])

# Define the number of sequences to generate
num_sequences <- 6

# Define the Jukes-Cantor model parameter
mu <- 5

# Define the alphabet of nucleotides
nucleotides <- c("A", "T", "C", "G")

# Function to generate random branch lengths from exponential distribution with parameter 10
generate_branch_lengths <- function(num_branches) {
  rexp(num_branches, rate = 10)
}

# Generate random branch lengths
branch_length <- generate_branch_lengths(num_sequences - 2)

# Function to calculate Jukes-Cantor model probabilities with custom transition probability matrix
calculate_transition_probs <- function(branch_length, mu) {
  # Initialize an empty list to store transition probability matrices
  transition_matrices <- list()
  
  # Loop over each branch length
  for (bl in branch_length) {
    # Compute the probabilities for same and different nucleotides
    p_same <- 1 - 3/4 * exp(-4/3 * mu * bl)
    p_diff <- 1/4 * exp(-4/3 * mu * bl)
    
    # Create the transition probability matrix directly
    transition_matrix <- matrix(c(p_same, p_diff, p_diff, p_diff,
                                  p_diff, p_same, p_diff, p_diff,
                                  p_diff, p_diff, p_same, p_diff,
                                  p_diff, p_diff, p_diff, p_same), 
                                nrow = 4, byrow = TRUE)
    
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

# Generate new sequence using the previous sequence
new_sequence <- generate_new_sequence(previous_sequence, mu, branch_length)
new_sequence

library(ape)

# Generate sequences for the branches
branch1 <- generate_new_sequence(base_sequence, mu, branch_length)
branch2 <- generate_new_sequence(base_sequence, mu, branch_length)
branch3 <- generate_new_sequence(base_sequence, mu, branch_length)
branch4 <- generate_new_sequence(base_sequence, mu, branch_length)
branch5 <- generate_new_sequence(base_sequence, mu, branch_length)
branch6 <- generate_new_sequence(base_sequence, mu, branch_length)


# Define the Newick tree string with branch lengths
# Concatenate branch lengths to the Newick tree string
tree_newick <- paste0("(((", branch1, ":",1- branch_length[1], ",", branch2, ":",1- branch_length[1], "):", branch_length[2], ",", 
                      "(", branch3, ":",1 -  branch_length[4], ",", branch4, ":",1- branch_length[4], "):", branch_length[4], "):", branch_length[2], ",", 
                      "(", branch5, ":", 1 -branch_length[3], ",", branch6, ":", 1- branch_length[3], "):", branch_length[4], ");")

# Create the phylogenetic tree object
tree <- read.tree(text = tree_newick)

# Plot the tree with adjusted parameters
plot(tree, cex = 0.8, label.offset = 0.01)

tree$edge.length
nodeHeights(tree)
tree$edge



