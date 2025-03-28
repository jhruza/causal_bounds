#generate random but valid probabilities based on nvals and names of a graph
valid_p_sample <- function(graph, return_joint=FALSE){
    nvals <- V(graph)$nvals
    names <- V(graph)$name

    prob <- list()
    for (i in 1:vcount(graph)){
      incoming_edges <- incident (graph, names[i], mode="in")
      incoming_names <- ends(graph, incoming_edges)[,1]
      incoming_positions <- which(names %in% incoming_names)

      col_names <- c("outcome", incoming_names ,"p")
      df <- data.frame(matrix(ncol = length(col_names), nrow = 0))
      colnames(df) <- col_names

      seq_list <- lapply(nvals[incoming_positions], function(x) seq(0, x-1))
      grid <- do.call(expand.grid, seq_list)
      
      if (nrow(grid) == 0){
        random_numbers <- runif(nvals[i])
        random_numbers <- random_numbers / sum(random_numbers) 
        df <- data.frame(cbind(0:(nvals[i]-1), random_numbers))
        colnames(df)<-col_names
      }

      for (row in seq_len(nrow(grid))){   
        random_numbers <- runif(nvals[i])
        random_numbers <- random_numbers / sum(random_numbers) 
        cond_df<-cbind(0:(nvals[i]-1), grid[row,], random_numbers)
        colnames(cond_df) <- col_names
        df<- rbind(df, cond_df)
        }
      prob[[names[i]]] <- df
    }
if (return_joint){
  seq_list <- lapply(nvals, function(x) seq(0, x - 1))
  grid <- do.call(expand.grid, seq_list)
  colnames(grid) <- names
  grid$p <- NA
  for (row in seq_len(nrow(grid))){
    p <- 1
    for (node in names){
      node_df <- prob[[node]]
      # Filter conditions - start with matching the outcome
      conditions <- node_df$outcome == grid[row, node]

      # Add conditions for parent values
      parent_cols <- setdiff(colnames(node_df), c("outcome", "p"))
      for (parent in parent_cols) {
        parent_val <- grid[row, parent]
        conditions <- conditions & node_df[[parent]] == parent_val
      }

      # Get probability from the matching row
      matched_prob <- node_df$p[conditions]

      # Multiply to get joint probability
      if (length(matched_prob) == 1) {
        p <- p * matched_prob
      } else {
        stop(paste("Error finding unique probability for node", node))
      }
    }
    grid$p[row] <- p
  }
  return(grid)
}
else {
  return(prob)
}
}

