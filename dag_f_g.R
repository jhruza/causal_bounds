library(causaloptim)

# generate graph_full
graph_no_iv <- initialize_graph(igraph::graph_from_literal(A-+Y, Ur-+A, Ur-+ B, Ur-+ G, Ur-+ Y, G -+ A))
plot(graph_no_iv)
V(graph_no_iv)$name
# Assuming V(graph_full)$name = "A"  "Y"  "Ur" "B" "G" 
V(graph_no_iv)$name
V(graph_no_iv)$leftside <- c(0, 0, 0, 0, 0)
V(graph_no_iv)$latent   <- c(0, 0, 1, 0, 0)
# V(graph_no_iv)$nvals    <- c(2, 2, 2, 2, 2)
V(graph_no_iv)$nvals    <- c(3, 2, 2, 3,3)
# V(graph_no_iv)$nvals    <- c(4, 2, 2, 4)


# E(graph_full)$rlconnect     <- c(0, 0, 0, 0, 0)
# E(graph_full)$edge.monotone <- c(0, 0, 0, 0, 0)

# to_bound <- "p{Y(A = 0) = 1; B=0}+p{Y(A = 1) = 1; B=1} - p{Y(A=0) = 1; G=0} + p{Y(A = 1) = 1; G=1}"
to_bound <- "p{Y(A = 0) = 1; B=0}+p{Y(A = 1) = 1; B=1} + p{Y(A = 2) = 1; B=2} - p{Y(A=0) = 1; G=0} - p{Y(A = 1) = 1; G=1} - p{Y(A=2) = 1; G=2}"
# to_bound <- "p{Y(A = 0) = 1; B=0}+p{Y(B = 1) = 1; B=1}+p{Y(B = 2) = 1; B=2} + p{Y(B = 3) = 1; B=3}"
# to_bound <- "p{Y(A = 0) = 1; B=0}"
obj_no_iv <- analyze_graph(graph_no_iv, constraints = NULL, effectt = to_bound)
bounds_no_iv <- optimize_effect_2(obj_no_iv)
cat(latex_bounds(bounds = bounds_no_iv$bounds, parameters = obj_no_iv$parameters, prob.sym = "P"))



# IV setting
graph_iv <- initialize_graph(igraph::graph_from_literal(Z-+A, A-+Y, Ur-+A, Ur-+ B, Ur-+ G, Ur-+ Y, G -+ A))
plot(graph_iv)
V(graph_iv)$name
# Assuming V(graph)$name = "Z"  "A"  "Y"  "Ur" "B" "G"
V(graph_iv)$name
V(graph_iv)$leftside <- c(1, 0, 0, 0, 0, 0)
V(graph_iv)$latent   <- c(0, 0, 0, 1, 0, 0)
# V(graph_iv)$nvals    <- c(2, 2, 2, 2, 2, 2)
V(graph_iv)$nvals    <- c(2, 3, 2, 2, 3, 3)


obj_iv <- analyze_graph(graph_iv, constraints = NULL, effectt = to_bound)
bounds_iv <- optimize_effect_2(obj_iv)
cat(latex_bounds(bounds = bounds_iv$bounds, parameters = obj_iv$parameters, prob.sym = "P"))
