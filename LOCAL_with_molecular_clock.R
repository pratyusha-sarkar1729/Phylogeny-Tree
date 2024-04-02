library(igraph)
library(phytools)


propose_new_tree <- function(tree, sequence_names) {
  proposed_tree <- tree
  current_tree <- tree
  # Choose a random internal branch
  # Get the edge matrix directly from the tree object
  edge_matrix <- current_tree$edge

  # Get the leaves
  leaves <- setdiff(unique(edge_matrix), edge_matrix[, 1])
  
  
  # Filter out internal edges
  internal_edges <- edge_matrix[edge_matrix[, 2] > Ntip(tree), ]
  
  # Sample a random row from internal_edges
  random_row <- sample(1:nrow(internal_edges), 1)
  
  # Get the randomly selected internal edge
  selected_internal_edge <- internal_edges[random_row, ]
  
  # Extract node IDs from the selected internal edge
  u <- selected_internal_edge[[2]]
  v <- selected_internal_edge[[1]]
  
  # Modify the tree topology
  # Get children of nodes u and v
  children_u <- tree$edge[which(tree$edge[, 1] == u),][,2]
  children_v <- tree$edge[which(tree$edge[, 1] == v & tree$edge[, 2] != u), , drop = FALSE][,2]
  children_all <- c(children_u,children_v)
  
  # Get the heights of nodes a, b, and c above v
  heights <- nodeHeights(proposed_tree)
  heights_u <- (heights[which(tree$edge[, 1] == u),][,2]) - (heights[which(tree$edge[, 1] == u),][,1])
  heights_v <- heights[which(tree$edge[, 1] == v & tree$edge[, 2] != u), , drop = FALSE][, 2] - heights[which(tree$edge[, 1] == v & tree$edge[, 2] != u), , drop = FALSE][, 1]

  height_children <- c(heights_u,heights_v)

  filtered_heights <- height_children[height_children > 0]
  
  # Sort the filtered heights
  sorted_heights <- unique(sort(filtered_heights))
  
  h1 <- sorted_heights[1]
  h2 <- sorted_heights[2]
  h3 <- sorted_heights[3]
  
  x <- runif(1) * h1
  y <- runif(1) * h2
  
  proposed_height_u = max(x,y)
  proposed_height_v = min(x,y)
  
  
  # Convert phylogenetic tree to graph object
  proposed_tree_graph <- as.igraph(proposed_tree)

    # Check if the proposed position of u is closer to v
  if (proposed_height_u < h1) {
    # Connect one of the three children nodes to v and connect the others to u
    # Randomly choose one child to connect to v
    child_to_connect <- sample(children_all, 1)
    # Connect the other children to u
    children_to_connect_to_u <- setdiff(children_all, child_to_connect)
    
    # Connect nodes
    clade_tree_u <- extract.clade(proposed_tree,u)
    clade_tree_v <- extract.clade(proposed_tree,children_v)
    clade_tree_u$root.edge <- clade_tree_v$root.edge <- a <- .2
    z = bind.tree(clade_tree_u, clade_tree_v,po = 0.2)
    a = setdiff(star_tree$tip.label, c(clade_tree_u$tip.label, clade_tree_v$tip.label))
    plot(z)
    plot(bind.tip(z,a,edge.length = branch_length_new[1],po = .02))
  } else {
    # Force the tree topology and connect the lowest child to v
    child_to_connect <- children_all[which.min((height_children))]
    # Connect the other children to u
    children_to_connect_to_u <- setdiff(children_all, child_to_connect)
    
    # Connect nodes
    proposed_tree_graph <- delete_edges(proposed_tree_graph, which(proposed_tree_graph$edges[,1] == v & proposed_tree_graph$edges[,2] == child_to_connect))
    proposed_tree_graph <- delete_edges(proposed_tree_graph, which(proposed_tree_graph$edges[,1] == u & proposed_tree_graph$edges[,2] == child_to_connect))
    proposed_tree_graph <- add_edges(proposed_tree_graph, c(v, child_to_connect))
    proposed_tree_graph <- add_edges(proposed_tree_graph, c(u, children_to_connect_to_u))
  }
  
  # Convert back to phylogenetic tree object
  proposed_tree <- igraph_to_phylo(proposed_tree_graph)
  
  # Acceptance/Rejection stage 
  
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




