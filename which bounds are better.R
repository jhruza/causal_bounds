library(causaloptim)
source("helper_function.R")

#define parameters
N_A <- 4
N_X <- 10
N_Y <- 2
N_Z <- 2
N_Ur <- 2

simulation_n <- 10000

#ITR function
f <- function(x, preimage = FALSE) {
  mapping <- list(
    `0` = c(0, 1),
    `1` = c(2, 3, 8, 9),
    `2` = c(4, 5, 6, 7)
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
# Assuming V(graph_full)$name =  "Z"  "A"  "X"  "Y" "Ur"
V(graph_full)$leftside <- c(1, 0, 0, 0, 0)
V(graph_full)$latent   <- c(0, 0, 0, 0, 1)
V(graph_full)$nvals    <- c(2, 3, 6, 2, 3)
V(graph_full)$nvals    <- c(N_Z, N_A, N_X, N_Y, N_Ur)

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
V(graph_cond)$nvals    <- c(N_Z, N_A, N_Y, N_Ur)


# define reduced graph
graph_reduced <- initialize_graph(igraph::graph_from_literal(Z -+ A, A-+Y, 
                                    Ur -+ B, Ur-+A, Ur -+ Y))
plot(graph_reduced)
V(graph_cond)$name

# Assuming V(graph)$name =  "Z"  "A"  "Y"  "Ur" "B" 
V(graph_reduced)$leftside <- c(1, 0, 0, 0, 0)
V(graph_reduced)$latent   <- c(0, 0, 0, 1, 0)
V(graph_reduced)$nvals    <- c(N_Z, N_A, N_Y, N_Ur, N_A)



  bound_a<-list()
  obj_a<-list()
  bounds_a<-list()
  for(a in 0:(V(graph_full)["A"]$nvals-1)){
      bound_a[paste0(a)]<- paste0("p{Y(A = ", a, ") = 1}")
      obj_a[[paste0(a)]]<- analyze_graph(graph_cond, constraints = NULL, effectt = bound_a[paste0(a)])
      bounds_a[[paste0(a)]] <- optimize_effect_2(obj_a[[paste0(a)]])
  }

# calculate boudns using approach where we condition on X
#calculate bounds for each value of A that 
calculate_bound_conditional <- function(sample_joint){

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
  return(result_cond)
}

# calculate true estimand E(f(X)) = sum_a{p(Y(A=a)=1, f(X)=a)}
calculate_true_estimand <- function(sample_joint) {
  result_estimand <- 0
  for (a in 0:(V(graph_full)["A"]$nvals - 1)){
    # First calculate p(f(X)=a)
    p_fx_equals_a <- 0
    for (x in f(a, preimage = TRUE)){
      p_fx_equals_a <- p_fx_equals_a + sum(subset(sample_joint, X == x)[, "p"])
    }
    
    # Calculate p(Y(A=a)=1 | f(X)=a)
    first_term <- 0
    for (x in f(a, preimage = TRUE)){
      for (u in 0:(V(graph_full)["Ur"]$nvals-1)){
        # p(X=x, Ur=u | f(X)=a)
        p_xu_given_fx_equals_a <- sum(subset(sample_joint, X == x & Ur == u)[, "p"]) / p_fx_equals_a

        # p(Y=1 | do(A=a), X=x, Ur=u) = p(Y=1 | A=a, X=x, Ur=u)
        p_y1_given_do_a_xu <- sum(subset(sample_joint, Y == 1 & A == a & X == x & Ur == u)[, "p"]) / 
                              sum(subset(sample_joint, A == a & X == x & Ur == u)[, "p"])

        first_term <- first_term + p_xu_given_fx_equals_a * p_y1_given_do_a_xu
      }
    }
    
    # p(Y(A=a)=1, f(X)=a) = p(Y(A=a)=1 | f(X)=a) * p(f(X)=a)
    result_estimand <- result_estimand + first_term * p_fx_equals_a
  }
  return(result_estimand)
}

  bound_reduced <- paste0("p{Y(A = ", 0: (V(graph_full)["A"]$nvals - 1), ") = 1; B = ", 0: (V(graph_full)["A"]$nvals - 1),"}", collapse = " + ") #nolint
  obj_reduced <- analyze_graph(graph_reduced, constraints = NULL, effectt = bound_reduced)
  bounds_reduced <- optimize_effect_2(obj_reduced)
# calculate bounds by reducing the full graph
calculate_bound_reduced <- function(sample_joint){

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
  return(result_reduced)
}



# calculate tight bounds for the full graph
calculate_tight_bounds <- function(sample_joint){
  bound_tight <- paste0("p{Y(A = ", sapply(0:(N_X-1), f), ") = 1; X = ", 0:(N_X-1), "}", collapse = " + ") #nolint
  obj_tight <- analyze_graph(graph_full, constraints = NULL, effectt = bound_tight)
  bounds_tight <- optimize_effect_2(obj_tight)

  #we need the prob_list_tight in the form of P(A = a, X = x, Y = y | Z = z)

  prob_list_tight <- setNames(as.list(rep(NA, length(obj_tight$parameters))), obj_tight$parameters)

  for (param in names(prob_list_tight)){
    a <- as.numeric(substring(param, 2, 2))
    x <- as.numeric(substring(param, 3, 3))
    y <- as.numeric(substring(param, 4, 4))
    z <- as.numeric(substring(param, 6, 6))

    prob_list_tight[param] <- sum(subset(sample_joint, A == a & X == x & Y == y & Z == z)[, "p"]) / 
                              sum(subset(sample_joint, Z == z)[, "p"])

  }
  result_tight <- do.call(bounds_tight$bounds_function, prob_list_tight)
  return(result_tight)
}



#sample from full graph
sample_joint <- valid_p_sample(graph_full, return_joint = TRUE)

print(calculate_true_estimand(sample_joint))
print(calculate_bound_reduced(sample_joint))
print(calculate_bound_conditional(sample_joint))
# print(calculate_tight_bounds(sample_joint))


# SINGLE CORE COMPUTING
# true_estimand<-list()
# bound_reduced<-list()
# bound_conditional<-list()
# pb <- txtProgressBar(min = 0, max = simulation_n, style = 3)
# for (i in 1:simulation_n){
#   setTxtProgressBar(pb, i)
#   sample_joint <- valid_p_sample(graph_full, return_joint = TRUE)
#   true_estimand[[i]] <- calculate_true_estimand(sample_joint)
#   bound_reduced[[i]] <- calculate_bound_reduced(sample_joint)
#   bound_conditional[[i]] <- calculate_bound_conditional(sample_joint)
# }

# true_estimand_values <- sapply(true_estimand, function(x) x)
# bound_reduced_lower <- sapply(bound_reduced, function(x) x$lower)
# bound_reduced_upper <- sapply(bound_reduced, function(x) x$upper)
# bound_conditional_lower <- sapply(bound_conditional, function(x) x$lower)
# bound_conditional_upper <- sapply(bound_conditional, function(x) x$upper)



#parallel computing

library(parallel)
library(pbapply)

# Define the number of cores to use
num_cores <- detectCores() - 1  # Use one less than the total number of cores

# Parallelize the simulation
cl <- makeCluster(num_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, library(causaloptim))  # Load necessary libraries on each worker

results <- pblapply(cl=cl, 1:simulation_n, function(i) {
  sample_joint <- valid_p_sample(graph_full, return_joint = TRUE)
  return(list(
    true_estimand = calculate_true_estimand(sample_joint),
    bound_reduced = calculate_bound_reduced(sample_joint),
    bound_conditional = calculate_bound_conditional(sample_joint)
  ))
})

stopCluster(cl)

#end parallel computing


# Extract results
true_estimand_values <- sapply(results, function(x) x$true_estimand)
bound_reduced_lower <- sapply(results, function(x) x$bound_reduced$lower)
bound_reduced_upper <- sapply(results, function(x) x$bound_reduced$upper)
bound_conditional_lower <- sapply(results, function(x) x$bound_conditional$lower)
bound_conditional_upper <- sapply(results, function(x) x$bound_conditional$upper)


# Compare if true estimand is less than the lower bound
sum(true_estimand_values < bound_reduced_lower & true_estimand_values > bound_reduced_upper)
sum(true_estimand_values < bound_conditional_lower & true_estimand_values > bound_conditional_upper)

min(true_estimand_values-bound_reduced_lower)
min(true_estimand_values-bound_conditional_lower)
min(bound_reduced_upper-true_estimand_values)
min(bound_conditional_upper-true_estimand_values)

# Check if the reduced bounds are better 
sum((bound_reduced_lower - bound_conditional_lower)>1e-10)
sum((bound_conditional_upper - bound_reduced_upper)>1e-10) 

# Save the variables as an RDS file
saveRDS(list(
  true_estimand_values = true_estimand_values,
  bound_reduced_lower = bound_reduced_lower,
  bound_reduced_upper = bound_reduced_upper,
  bound_conditional_lower = bound_conditional_lower,
  bound_conditional_upper = bound_conditional_upper
), file = "results.rds")
