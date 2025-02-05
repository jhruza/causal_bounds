library(causaloptim)

# generate graph_full
graph_full <- initialize_graph(igraph::graph_from_literal(Z -+ X, S-+X, X-+ Y, 
                                    Ur -+ X, Ur -+ Y))
# print(V(graph))
# plot(graph)
# Assuming V(graph)$name = ["Z", "X", "S", "Y", "Ur"]
V(graph_full)$leftside <- c(1, 0, 1, 0, 0)
V(graph_full)$latent   <- c(0, 0, 0, 0, 1)
V(graph_full)$nvals    <- c(2, 2, 2, 2, 2)
# E(graph_full)$rlconnect     <- c(0, 0, 0, 0, 0)
# E(graph_full)$edge.monotone <- c(0, 0, 0, 0, 0)

to_bound <- "p{Y(X = 0) = 1}"
obj_full <- analyze_graph(graph_full, constraints = NULL, effectt = to_bound)
bounds_full <- optimize_effect_2(obj_full)
cat(latex_bounds(bounds = bounds_full$bounds, parameters = obj_full$parameters, prob.sym = "P"))



# generate margenalized graph
graph_marg <- initialize_graph(igraph::graph_from_literal(Z -+ X, X-+ Y, 
                                    Ur -+ X, Ur -+ Y))


# Assuming V(graph)$name = ["Z", "X", "Y", "Ur"]
V(graph_marg)$leftside <- c(1, 0, 0, 0)
V(graph_marg)$latent   <- c(0, 0, 0, 1)
V(graph_marg)$nvals    <- c(2, 2, 2, 2)
# E(graph_full)$rlconnect     <- c(0, 0, 0, 0)
# E(graph_full)$edge.monotone <- c(0, 0, 0, 0)
obj_marg <- analyze_graph(graph_marg, constraints = NULL, effectt = to_bound)
bounds_marg <- optimize_effect_2(obj_marg)



valid_p_sample <- function(obj){
    nvals <- V(obj$graph)$nvals
    name <- V(obj$graph)$name
    parameters <- obj$parameters
    #create named list with NA as entries
    prob_list <- setNames(as.list(rep(NA, length(parameters))), parameters)
    rightvars <- attr(parameters, "rightvars")
    condvars <-attr(parameters, "condvars")
    positions_rightvars <- which(name %in% rightvars)
    total_vals_rightvars <- prod(nvals[positions_rightvars]) 

    # looop list with NA entries and fill with valid probabilities
    for (param in parameters){
      if(is.na(prob_list[param])){
        random_numbers <- runif(total_vals_rightvars)
        random_numbers <- random_numbers / sum(random_numbers) 

        cond_instance <- strsplit(param, "_")[[1]][2]
        selected_elements <- grep(paste0(cond_instance, "$"), parameters)
        prob_list[selected_elements] <- random_numbers
        }
    }
return(prob_list)
}




p <- valid_p_sample(obj_full)
attr(obj_full$parameters,"condvars")
nval <-2
random_s <- runif(2)
random_s <- random_s / sum(random_s) 

calc_bounds_marg<-c(lower = 0, upper =0)
for (val in 1:nval) {
    index <- grep(paste0(val-1, "$"), names(p))
    marg_p<-p[index]
    names(marg_p) <- sapply(names(marg_p), function(x) substr(x, 1, nchar(x) - 1))
    calc_bounds_marg <- calc_bounds_marg + random_s[val] * do.call(bounds_marg$bounds_function, marg_p)
}

calc_bounds_full <- do.call(bounds_full$bounds_function, p)

calc_bounds_marg
calc_bounds_full
calc_bounds_marg[1]<calc_bounds_full[1] & calc_bounds_marg[2]>calc_bounds_full[2]


