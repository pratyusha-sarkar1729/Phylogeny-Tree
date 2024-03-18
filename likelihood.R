library(phangorn)



# Define the sequence names
sequence_names <- c("A", "B", "C", "D", "E", "F")

# Extract sequences from the tree
leaf_nodes <- tree$tip.label

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
tree$tip.label <- tip_labels

# Apply PML estimation
result <- pml(tree, data = D2,k = 4)
result
