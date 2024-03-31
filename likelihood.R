library(phangorn)



# Define the sequence names
sequence_names <- c("A", "B", "C", "D", "E", "F")

# Extract sequences from the tree
leaf_nodes <- star_tree$tip.label

# Prepare the output
output <- paste0(">", sequence_names, "\n", leaf_nodes, "\n")

# Write sequences to a temporary FASTA file
writeLines(output, "temp_alignment.fasta")

# Read the FASTA file into a DNAbin object
D <- read.dna("temp_alignment.fasta", format = "fasta")

# Check the resulting DNAbin object
D

D2 = as.phyDat(D)

# Define tip labels
tip_labels <- c("A", "B", "C", "D", "E", "F")

# Assign tip labels to the tree
star_tree$tip.label <- tip_labels

# Apply PML estimation
result <- pml(star_tree, data = D2,k = 4)
result

#############function for likelihood ##########
###############################################

calculate_likelihood <- function(tree, sequence_names) {
  # Extract sequences from the tree
  leaf_nodes <- tree$tip.label
  
  # Prepare the output
  output <- paste0(">", sequence_names, "\n", leaf_nodes, "\n")
  
  # Write sequences to a temporary FASTA file
  writeLines(output, "temp_alignment.fasta")
  
  # Read the FASTA file into a DNAbin object
  D <- read.dna("temp_alignment.fasta", format = "fasta")
  
  # Convert DNAbin object to phyDat object
  D2 <- as.phyDat(D)
  
  # Define tip labels
  tip_labels <- sequence_names
  
  # Assign tip labels to the tree
  tree$tip.label <- tip_labels
  
  # Apply PML estimation
  result <- pml(tree, data = D2, k = 4)
  
  # Return the likelihood estimation result
  return(result)
}


# Define sequence names based on the number of nucleotides
sequence_names <- c("A", "B", "C", "D", "E", "F")

# Call the function with your tree object and sequence names
likelihood_result <- calculate_likelihood(star_tree, sequence_names)

# View the likelihood estimation result
likelihood_result
