#generate random but valid probabilities based on nvals and names
valid_p_sample <- function(obj){
    nvals <- V(obj$graph)$nvals
    names <- V(obj$graph)$name
    parameters <- obj$parameters

    prob <- list()
    for (i in 1:vcount(obj$graph)){
      incoming_edges <- incident (obj$graph, names[i], mode="in")
      incoming_names <- ends(obj$graph, incoming_edges)[,1]
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
return(prob)
}