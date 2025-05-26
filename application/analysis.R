library(causaloptim)


data <- readRDS("application/analysis_data.rds")
# Create a sub-dataframe of 'data' and rename some columns
df <- data[, c("trt", "stratum", "consume_group", "outcome")]  # replace with actual column names
# Explicitly map "a" to 1 and "b" to 2 in consume_group
df$consume_group <- ifelse(df$consume_group == "less than .2g", 0,
                           ifelse(df$consume_group == ".2 to 6", 1,
                                  ifelse(df$consume_group == "6 or more", 2, NA)))
# Map IV veriable to 0 and 1, trt is randomizied treatment
df$trt <- ifelse(df$trt == "Peanut Avoidance", 0,
                ifelse(df$trt == "Peanut Consumption", 1, NA))
df$stratum <- df$stratum - 1  # Adjust stratum to be 0-indexed
df$outcome <- ifelse(df$outcome == TRUE, 1,
                     ifelse(df$outcome == FALSE, 0, NA))
colnames(df) <- c("Z", "X", "A", "Y")


#ITR function 
## treatment rule: if stratum = 1 give A = 0.2 to 6
##                 if stratum = 2 give  A = 6 or more
f <- function(x, preimage = FALSE) {
  mapping <- list(
    `1` = c(0,1),
    `2` = c(1)
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
#standard of care
g <- function(x, preimage = FALSE) {
  mapping <- list(
    `0` = c(0,1)
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

# calculate the treatment according to the ITR function
for (i in seq_len(nrow(df))) {
  df$F[i] <- f(df$X[i], preimage = FALSE)
  df$G[i] <- g(df$X[i], preimage = FALSE)
}
# Standard of Care is to give treatment A=0


#define parameters
N_Z <- 2
N_A <- 3
N_X <- 2
N_G <- 3
N_F <- 3
N_Y <- 2

# define full graph
graph_full <- initialize_graph(igraph::graph_from_literal(Z -+ A, X-+A, A -+ Y,  X-+ F, X-+G , G-+A, X-+ Y, 
                                                          Ur -+ X, Ur-+A, Ur -+ Y))

# print(V(graph))
plot(graph_full)
V(graph_full)$name
# Assuming V(graph_full)$name = "X"  "A"  "Y"  "F"  "G"  "Ur"
V(graph_full)$leftside <- c(1, 0, 0, 0, 0, 0, 0)
V(graph_full)$latent   <- c(0, 0, 0, 0, 0, 0, 1)
V(graph_full)$nvals    <- c(N_Z, N_X, N_A, N_Y, N_F, N_G, 2)


#function to calculate probabilities for given probability list from bounds and dataframe 
#conditioned is a named list if there is a variable which should be conditioned on, most differ from keys_string
calculate_prob <- function(prob_list, key_string, df, conditioned = list()) {
  for (param_name in names(prob_list)) {
    # 1. Parse param_name (e.g., "p01_3")
    # It should follow the "p<joint_digits>_<conditional_digits>" format
    if (!startsWith(param_name, "p")) {
      warning(paste("Parameter name '", param_name, "' does not start with 'p'. Assigning NA.", sep = ""))
      prob_list[[param_name]] <- NA
      next
    }

    value_part_string <- sub("^p", "", param_name) # e.g., "01_3"

    underscore_pos <- regexpr("_", value_part_string, fixed = TRUE)

    if (underscore_pos == -1) {
      warning(paste("Parameter name '", param_name, "' does not contain '_' for conditioning. Assigning NA.", sep = ""))
      prob_list[[param_name]] <- NA
      next
    }

    joint_values_str <- substr(value_part_string, 1, underscore_pos - 1)
    conditional_values_str <- substr(value_part_string, underscore_pos + 1, nchar(value_part_string))

    if (nchar(joint_values_str) == 0) {
      warning(paste("Parameter name '", param_name, "' has an empty joint part (before '_'). Assigning NA.", sep = ""))
      prob_list[[param_name]] <- NA
      next
    }

    if (nchar(conditional_values_str) == 0) {
      warning(paste("Parameter name '", param_name, "' has an empty conditional part (after '_'). Assigning NA.", sep = "")) # nolint
      prob_list[[param_name]] <- NA
      next
    }

    # 2. Convert value strings to numeric vectors
    joint_individual_values_str <- strsplit(joint_values_str, "")[[1]]
    joint_numeric_values <- as.numeric(joint_individual_values_str)

    conditional_individual_values_str <- strsplit(conditional_values_str, "")[[1]] # nolint
    conditional_numeric_values <- as.numeric(conditional_individual_values_str)

    if (anyNA(joint_numeric_values) || anyNA(conditional_numeric_values)) {
      warning(paste("Could not parse all numeric values from '", param_name, "'. Assigning NA.", sep = ""))
      prob_list[[param_name]] <- NA
      next
    }

    # 3. Determine corresponding keys from key_string
    num_joint_vars <- length(joint_numeric_values)
    num_conditional_vars <- length(conditional_numeric_values)

    if (nchar(key_string) != (num_joint_vars + num_conditional_vars + 1)) {
      warning(paste("Length of key_string '", key_string,
                    "' (", nchar(key_string), ") does not match total number of variables in param '", param_name, 
                    "' (", num_joint_vars + num_conditional_vars, "). Assigning NA.", sep = ""))
      prob_list[[param_name]] <- NA
      next
    }

    all_keys_vec <- strsplit(key_string, "")[[1]]

    joint_keys <- if (num_joint_vars > 0) all_keys_vec[1:num_joint_vars] else character(0)
    conditional_keys <- if (num_conditional_vars > 0) all_keys_vec[(num_joint_vars + 2):(num_joint_vars + num_conditional_vars+1)] else character(0) # nolint
    
    #3.5 add the values which come from conditioned
    if (length(intersect(names(conditioned), conditional_keys)) > 0) {
      warning(paste("There are shared variables to be conditioned on in '", names(conditioned),
                    "' and '", conditional_keys, sep = ""))
      prob_list[[param_name]] <- NA
      next
    }
    conditional_keys <- c(conditional_keys, names(conditioned))
    conditional_numeric_values <- c(conditional_numeric_values, unlist(unname(conditioned)))

    # 4. Calculate probability P(joint_vars | conditional_vars) from df

    # Build logical vector for the conditional part
    current_selection_cond <- rep(TRUE, nrow(df))
    for (i in seq_along(conditional_keys)) {
      key_col <- conditional_keys[i]
      val <- conditional_numeric_values[i]
      if (!key_col %in% names(df)) {
          stop(paste("Conditional key '", key_col, "' from key_string is not a column in df. Please check df column names.", sep="")) # nolint
      }
      current_selection_cond <- current_selection_cond & (df[[key_col]] == val)
    }
    df_subset_conditional <- df[current_selection_cond, , drop = FALSE]
    count_denominator <- nrow(df_subset_conditional)

    if (count_denominator == 0) {
      prob_list[[param_name]] <- NaN # P(A|B) is undefined if P(B)=0
    } else {
      # Build logical vector for the joint part, applied to the conditionally subsetted df
      current_selection_joint <- rep(TRUE, nrow(df_subset_conditional))
      for (i in seq_along(joint_keys)) {
        key_col <- joint_keys[i]
        val <- joint_numeric_values[i]
        if (!key_col %in% names(df_subset_conditional)) { # Should exist if in df
              stop(paste("Joint key '", key_col, "' from key_string is not a column in df. Please check df column names.", sep="")) # nolint
        }
        current_selection_joint <- current_selection_joint & (df_subset_conditional[[key_col]] == val)
      }
      # Count rows matching joint conditions within the conditional subset
      count_numerator <- sum(current_selection_joint) 

      prob_list[[param_name]] <- count_numerator / count_denominator
    }
  }
  return(prob_list)
}

# # define reduced graph
# graph_reduced <- initialize_graph(igraph::graph_from_literal(Z -+ A, A-+Y, G -+ A, 
#                                                              Ur -+ G, Ur -+ F, Ur-+A, Ur -+ Y))
# plot(graph_reduced)
# V(graph_reduced)$name

# # Assuming V(graph)$name = "Z" "A"  "Y"  "G"  "Ur" "F" 
# V(graph_reduced)$leftside <- c(1, 0, 0, 0, 0, 0)
# V(graph_reduced)$latent   <- c(0, 0, 0, 0, 1, 0)
# V(graph_reduced)$nvals    <- c(N_Z, N_A, N_Y, N_G, 2, N_F)


# objective_F <- paste0("p{Y(A = ", 0: (N_A - 1) , ") = 1; F = ", 0: (N_G - 1),"}", collapse = " + ") #nolint
# objective_G <- paste0("p{Y(A = ", 0: (N_G - 1 ), ") = 1; G = ", 0: (N_G - 1),"}", collapse = " - ") #nolint
# bound_reduced <- paste0(objective_F, " - ", objective_G)
# obj_reduced <- analyze_graph(graph_reduced, constraints = NULL, effectt = bound_reduced)
# bounds_reduced <- optimize_effect_2(obj_reduced)

# #create empty list of arguments needed for caluclate boudns
# prob_list_reduced <- setNames(as.list(rep(NA, length(obj_reduced$parameters))), obj_reduced$parameters)

# key_string <- attr(obj_reduced$parameters, "key")

# prob_list_reduced <- calculate_prob(prob_list_reduced, key_string, df)


# #sanity check
# sum(unlist(prob_list_reduced[seq(1, length(prob_list_reduced), 2)])) # should be 1

# latex_bounds(bounds_reduced$bounds, obj_reduced$parameters)


# define conditianal graph
graph_cond <- initialize_graph(igraph::graph_from_literal(Z-+A, A -+ Y,
                                                          Ur-+A, Ur -+ Y))
plot(graph_cond)

V(graph_cond)$name
# Assuming V(graph)$name =  "Z"  "A"  "Y"  "Ur"
V(graph_cond)$leftside <- c(1, 0, 0, 0)
V(graph_cond)$latent   <- c(0, 0, 0, 1)
V(graph_cond)$nvals    <- c(N_Z, N_A, N_Y, 2)




bound_idx <- list()
obj_idx <- list()
bounds_idx <- list()
for (var_a in 0:(N_A - 1)){
  for (var_g in 0:(N_G - 1)){
    idx <- paste0(var_a, "_", var_g)
    bound_idx[[idx]] <- paste0("p{Y(A = ", var_a, ") = 1} - p{Y(A = ", var_g, ") = 1}")
    obj_idx[[idx]] <- analyze_graph(graph_cond, constraints = NULL, effectt = bound_idx[[idx]])
    bounds_idx[[idx]] <- optimize_effect_2(obj_idx[[idx]])
  }
}

result_cond <- c(lower = 0, upper = 0)
for (x in 0:(N_X-1)){
  idx <- paste0(f(x), "_", g(x))
  bounds <- bounds_idx[[idx]]
  key_string <- attr(obj_idx[[idx]]$parameters, "key")
  prob_list <- setNames(as.list(rep(NA, length(obj_idx[[idx]]$parameters))), obj_idx[[idx]]$parameters)
  prob_list <- calculate_prob(prob_list, key_string, df, conditioned = list(X = x))

  count_values_x <- sum(df$X == x)

  p_x <- count_values_x / nrow(df) # P(X=x)

  if (f(x) == g(x)) {
    # If the treatment rule and standard of care are the same skip as bounds are zero
    next
  }

  result_cond <- result_cond + p_x * do.call(bounds$bounds_function, prob_list)
}

print(result_cond)



# #experimental bounding counterfactuals

# bound_idx <- list()
# obj_idx <- list()
# bounds_idx <- list()
# for (var_a in 0:(N_A - 1)){
#     idx <- paste0(var_a)
#     bound_idx[[idx]] <- paste0("p{Y(A = ", var_a, ") = 1}")
#     obj_idx[[idx]] <- analyze_graph(graph_cond, constraints = NULL, effectt = bound_idx[[idx]])
#     bounds_idx[[idx]] <- optimize_effect_2(obj_idx[[idx]])
# }

# result_cond <- c(lower = 0, upper = 0)
# for (x in 0:(N_X-1)){
#   idx <- paste0(2)
#   bounds <- bounds_idx[[idx]]
#   key_string <- attr(obj_idx[[idx]]$parameters, "key")
#   prob_list <- setNames(as.list(rep(NA, length(obj_idx[[idx]]$parameters))), obj_idx[[idx]]$parameters)
#   prob_list <- calculate_prob(prob_list, key_string, df, conditioned = list(X = x))

#   count_values_x <- sum(df$X == x)

#   p_x <- count_values_x / nrow(df) # P(X=x)

#   result_cond <- result_cond + p_x * do.call(bounds$bounds_function, prob_list)
# }
# print(result_cond)
