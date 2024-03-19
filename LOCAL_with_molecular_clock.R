library(igraph)
library(phytools)


# Define a function to propose a new phylogeny tree
propose_new_tree <- function(current_tree, sequence_names) {
  # Copy the current tree to the proposed tree
  proposed_tree <- current_tree
  
  # Step 1: Choose a random internal branch
  # Get the edge matrix directly from the tree object
  edge_matrix <- tree$edge
  
  # Filter out internal edges
  internal_edges <- edge_matrix[edge_matrix[, 2] > Ntip(tree), ]
  
  # Sample a random row from internal_edges
  random_row <- sample(1:nrow(internal_edges), 1)
  
  # Get the randomly selected internal edge
  selected_internal_edge <- internal_edges[random_row, ]
  
  # Extract node IDs from the selected internal edge
  u <- selected_internal_edge[[2]]
  v <- selected_internal_edge[[1]]
  
  # Step 2: Modify the tree topology
  # Get children of nodes u and v
  children_u <- tree$edge[which(tree$edge[, 1] == u),][,2]
  children_v <- tree$edge[which(tree$edge[, 1] == v & tree$edge[, 2] != u), , drop = FALSE][,2]
  children_all <- c(children_u,children_v)
  # Get the heights of nodes a, b, and c above v
  heights <- nodeHeights(proposed_tree)
  heights_u <- (heights[which(tree$edge[, 1] == u),][,2])
  heights_v <- heights[which(tree$edge[, 1] == v & tree$edge[, 2] != u), , drop = FALSE][, 2]

  height_children <- c(heights_u,heights_v)
  # Filter heights greater than 0
  filtered_heights <- height_children[height_children > 0]
  # Sort the filtered heights
  sorted_heights <- unique(sort(filtered_heights))
  h1 <- sorted_heights[1]
  h2 <- sorted_heights[2]
  h3 <- sorted_heights[3]
  
  # Compute the proposed distance between v and the nearest child
  proposed_distance <- runif(1) * h1
  
  # Compute the proposed location of u
  proposed_height_u <- runif(1) * h2
  
  # Check if the proposed position of u is closer to v
  if (proposed_height_u < proposed_distance) {
    # Connect one of the three children nodes to v and connect the others to u
    # Randomly choose one child to connect to v
    child_to_connect <- sample(children_v, 1)
    # Connect the other children to u
    children_to_connect_to_u <- setdiff(children_v, child_to_connect)
    
    # Connect nodes
    proposed_tree <- drop.tip(proposed_tree, children_to_connect_to_u)
    proposed_tree <- drop.tip(proposed_tree, v)
    proposed_tree <- bind.tree(proposed_tree, u, v)
    proposed_tree <- bind.tree(proposed_tree, children_to_connect_to_u, u)
  } else {
    # Force the tree topology and connect the lowest child to v
    child_to_connect <- children_all[which.min((height_children))]
    # Connect the other children to u
    children_to_connect_to_u <- setdiff(children_all, child_to_connect)
    
    # Connect nodes
    proposed_tree <- add.edge(proposed_tree, c(u, child_to_connect))
    proposed_tree <- add_edges(proposed_tree, c(u, children_to_connect_to_u))
    proposed_tree <- drop.tip(proposed_tree, node = v)
  }
  
  # Step 3: Acceptance/Rejection stage (Metropolis-Hastings)
  # Calculate the Metropolis-Hastings acceptance probability based on the likelihood of the proposed tree modification
  # Generate a random number from a uniform distribution
  # Accept or reject the proposed tree modification based on the acceptance probability
  
  # Calculate the likelihood of the proposed tree modification using the PML method
  likelihood_proposed <- calculate_likelihood(proposed_tree, sequence_names)
  
  # Calculate the likelihood of the current tree
  likelihood_current <- calculate_likelihood(current_tree, sequence_names)
  
  # Calculate the acceptance probability
  acceptance_prob <- exp(likelihood_proposed - likelihood_current)
  
  # Generate a random number from a uniform distribution
  random_number <- runif(1)
  
  # Accept or reject the proposed tree modification based on the acceptance probability
  if (random_number <= acceptance_prob) {
    # Accept the proposed tree modification
    current_tree <- proposed_tree
  }
  
  
  # Return the proposed tree
  return(proposed_tree)
}















## load tree
data(vertebrate.tree)
vertebrate.tree$edge
## compute height of all nodes
H<-nodeHeights(vertebrate.tree)
print(H)
## compute total tree depth
max(H)
