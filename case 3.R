library(causaloptim)
source("helper_function.R")

# generate graph_full
graph_full <- initialize_graph(igraph::graph_from_literal(Z -+ X, S -+X, S -+Y , X- + Y, 
                                    Ur -+ S, Ur -+ X, Ur -+ Y))
# print(V(graph))
plot(graph_full)
# Assuming V(graph)$name = ["Z", "X", "S", "Y", "Ur"]
V(graph_full)$leftside <- c(1, 0, 0, 0, 0)
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
cat(latex_bounds(bounds = bounds_marg$bounds, parameters = obj_marg$parameters, prob.sym = "P"))

# generate valid but random conditinal probabilities based on graph 
p <- valid_p_sample(obj_full)


rightvars <- attr(obj_full$parameters, "rightvars")
condvars <-attr(obj_full$parameters, "condvars")


prob_list <- setNames(as.list(rep(NA, length(obj_full$parameters))), obj_full$parameters)

# create list that matches arguments of bounds_function
for (param in obj_full$parameters){
    x<- as.numeric(substring(param, 2, 2))
    s<- as.numeric(substring(param, 3, 3))
    y<- as.numeric(substring(param, 4, 4))
    z<- as.numeric(substring(param, 6, 6))
    # get joint probability and margenalize out Ur
    joint_p<-sum(subset(p$X, outcome == x & S == s & Z == z)[ ,"p"]*subset(p$Y, outcome == y & X == x & S == s)[ ,"p"]*subset(p$S, outcome == s)[ ,"p"]*subset(p$Z, outcome == z)[ ,"p"]*subset(p$Ur)[ ,"p"]) # nolint
    # condition on Z
    cond_p<-joint_p/subset(p$Z, outcome == z)[ ,"p"]
    prob_list[param] <- cond_p
} 

# nvals of S
nval<-V(graph_full)$nvals[3]
result_marg<-c(lower = 0, upper =0)
for (val in 1:nval) {
    # filter values condition S=s
    index <- grep(paste0("^p.",val-1, "._."), names(prob_list))
    marg_p<-prob_list[index]

    # renaming list so parameters match that of bounds_function
    names(marg_p) <- sapply(names(marg_p), function(x) paste0(substr(x, 1, 2), substr(x, 4, nchar(x))))
    
    # result_marg + P(s)*Lowerbound(s)
    result_marg <- result_marg + sum(subset(p$Ur)[ ,"p"]*subset(p$S, outcome == val-1)[ ,"p"]) * do.call(bounds_marg$bounds_function, marg_p)
}

result_full<- do.call(bounds_full$bounds_function, prob_list)

print(result_full)
print(result_marg)

