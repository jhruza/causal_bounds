library(causaloptim)
source("helper_function.R")

f <- function(x, preimage = FALSE) {
  mapping <- list(
    `0` = c(0, 1),
    `1` = c(2, 3),
    `2` = c(4, 5)
  )
  
  if (preimage) {
    for (key in names(mapping)) {
      if (x == as.numeric(key)) {
        return(mapping[[key]])
      }
    }
    warning("Preimage not defined for this scenario")
    return(NA)
  } else {
    for (key in names(mapping)) {
      if (x %in% mapping[[key]]) {
        return(as.numeric(key))
      }
    }
    warning("Decision function not defined for this scenario")
    return(NA)
  }
}

# define full graph
graph_full <- initialize_graph(igraph::graph_from_literal(Z -+ A, X-+A, A-+Y , X-+ Y, 
                                    Ur -+ X, Ur-+A, Ur -+ Y))
# print(V(graph))
plot(graph_full)
V(graph_full)$name
# Assuming V(graph)$name =  "Z"  "A"  "X"  "Y"  "B"  "Ur"
V(graph_full)$leftside <- c(1, 0, 0, 0, 0)
V(graph_full)$latent   <- c(0, 0, 0, 0, 1)
V(graph_full)$nvals    <- c(2, 3, 6, 2, 3)
# E(graph_full)$rlconnect     <- c(0, 0, 0, 0, 0)
# E(graph_full)$edge.monotone <- c(0, 0, 0, 0, 0)


# define conditianal graph
graph_cond <- initialize_graph(igraph::graph_from_literal(Z -+ A, A-+Y, 
                                    Ur-+A, Ur -+ Y))
plot(graph_cond)
V(graph_cond)$name
# Assuming V(graph)$name =  "Z"  "A"  "Y"  "Ur"
V(graph_cond)$leftside <- c(1, 0, 0, 0)
V(graph_cond)$latent   <- c(0, 0, 0, 1)
V(graph_cond)$nvals    <- c(2, 3, 2, 2)


# define reduced graph
graph_reduced <- initialize_graph(igraph::graph_from_literal(Z -+ A, A-+Y, 
                                    Ur -+ B, Ur-+A, Ur -+ Y))
plot(graph_reduced)
V(graph_cond)$name

# Assuming V(graph)$name =  "Z"  "A"  "Y"  "Ur" "B" 
V(graph_reduced)$leftside <- c(1, 0, 0, 0, 0)
V(graph_reduced)$latent   <- c(0, 0, 0, 1, 0)
V(graph_reduced)$nvals    <- c(2, 3, 2, 2, 3)


#sample from full graph
sample <- valid_p_sample(graph_full)
sample_joint <- valid_p_sample(graph_full, return_joint = TRUE)


# calculate boudns using approach where we condition on X
#calculate bounds for each value of A that 
bound_a<-list()
obj_a<-list()
bounds_a<-list()
for(a in 0:(V(graph_full)["A"]$nvals-1)){
    bound_a[paste0(a)]<- paste0("p{Y(A = ", a, ") = 1}")
    obj_a[[paste0(a)]]<- analyze_graph(graph_cond, constraints = NULL, effectt = bound_a[paste0(a)])
    bounds_a[[paste0(a)]] <- optimize_effect_2(obj_a[[paste0(a)]])
}

result_cond <- c(lower = 0, upper = 0)
# marginalize over X
for (x in 0:(V(graph_full)["X"]$nvals-1)){
  #select the correct bounds based on f(i)
  bounds <- bounds_a[[paste0(f(x))]]
  #calculate the distribution for the bounds conditiond on X=x and save in prob_list
  #prob list must have syntax of P(A=a, Y=y | Z=z) which is in reality P(A=a, Y=y | Z=z, X=x)
  prob_list <- setNames(as.list(rep(NA, length(obj_a[[paste0(f(x))]]$parameters))), obj_a[[paste0(f(x))]]$parameters)
  for (param in names(prob_list)){
    a <- as.numeric(substring(param, 2, 2))
    y <- as.numeric(substring(param, 3, 3))
    z <- as.numeric(substring(param, 5, 5))

    #joint observed distribution marged over Ur  P(A=a, Y=y, Z=z, X=x)
    joint_p<- sum(subset(sample_joint, A == a & Y == y & Z == z & X==x)[, "p"])
    
    #joint observed distribution conditioned on Z=z and X=x .
    joint_p <- joint_p / sum(subset(sample_joint, Z==z & X==x)[ , "p"])
    prob_list[param] <- joint_p
  }
  # calculate p(X=x)
  p_x <- sum(subset(sample_joint, X==x)[ , "p"]) # nolint

  result_cond <- result_cond + p_x * do.call(bounds$bounds_function, prob_list)
}

# calculate true estimand E(f(X)) = sum_a{p(Y(A=a)=1, f(X)=a)}
result_estimand <- 0
for (a in 0:(V(graph_full)["A"]$nvals-1)){
  # calculate p(Y(A=a)=1, f(X)=a) = p(Y(A=a) = 1 | f(x)=a) * p(f(x)=a)

  #first term is p(Y(A=a)=1 | f(x)=a) = sum_{{x:f(x)=a}, {u}} {p(x,u) * p( y = 1 | a, x, u)}
  first_term <- 0
  for (x in f(a, preimage = TRUE)){
    for (u in 0:(V(graph_full)["A"]$nvals-1)){
      first_term <- first_term + sum(subset(sample_joint, X == x & Ur == u)[, "p"])*sum(subset(sample_joint, Y == 1 & A == a & X == x & Ur == u)[, "p"])/sum(subset(sample_joint, A == a & X == x & Ur == u)[, "p"]) # nolint
    }
  }
  #second term is p(f(x)=a) = sum_{{x:f(x)=a}, {u}} {p(x,u)}
  second_term <- 0
  for (x in f(a, preimage = TRUE)){
    second_term_partial <- sum(subset(sample_joint, X == x)[, "p"])
    second_term <- second_term + second_term_partial
  }
  result_estimand <- result_estimand + first_term * second_term
}


# calculate bounds by reducing the full graph
bound_reduced <- paste0("p{Y(A = ", 0: (V(graph_full)["A"]$nvals - 1), ") = 1; B = ", 0: (V(graph_full)["A"]$nvals - 1),"}", collapse = " + ") #nolint
obj_reduced <- analyze_graph(graph_reduced, constraints = NULL, effectt = bound_reduced)
bounds_reduced <- optimize_effect_2(obj_reduced)

#we need the prob_list_reduced in the form of P(A = a, Y = a, B = b | Z = z)
prob_list_reduced <- setNames(as.list(rep(NA, length(obj_reduced$parameters))), obj_reduced$parameters)

for (param in names(prob_list_reduced)){
  a <- as.numeric(substring(param, 2, 2))
  y <- as.numeric(substring(param, 3, 3))
  b <- as.numeric(substring(param, 4, 4))
  z <- as.numeric(substring(param, 6, 6))

  # P(A = a, Y = a, B = b | Z = z)= sum_{{x:f(x)=b} {P(A = a, Y = a, X = x | Z = z)}
  prob <- 0
  for (x in f(b, preimage = TRUE)){
    #P(A = a, Y = a, X = x | Z = z))
    prob_partial <- sum(subset(sample_joint, A == a & Y == y & X == x & Z == z)[, "p"])/ sum(subset(sample_joint, Z == z)[, "p"]) # nolint
    prob <- prob + prob_partial
  }
  prob_list_reduced[param] <- prob
}
result_reduced <- do.call(bounds_reduced$bounds_function, prob_list_reduced)

result_estimand
result_reduced
result_cond

