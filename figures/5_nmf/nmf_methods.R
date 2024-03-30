library(Matrix)
library('NMF')
library(ggplot2)
library(gtable)
library(data.table)

library(ggsignif)
library('ComplexHeatmap')
library('GEOquery')
library(linseed)
library(dplyr)
library(grid)
library(gridExtra)
library(ggrastr)
library(ggpubr)
library(svglite)
library(hNMF)
library(progress)
library(DualSimplex)

#### Data Generation methods ####
make_syntetic_data <- function(M, N, K ) {
  nmf_object <- syntheticNMF(M, K, N, factors = T, noise=TRUE)
  V <- nmf_object[[1]]
  rownames(V) <- c(paste0("feature", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  nmf_object$H
  true_W <- nmf_object$W
  colnames(true_H) <- colnames(V)
  rownames(true_W) <- rownames(V)
  return (list(V=V, true_H=true_H, true_W=true_W))
  
}

## Simplex sampling method. Can sample points for H
sampleFromSimplexUniformly <- function(n, k=3, M=100000, fixed_proportion=NULL) {
  X <- matrix(0, nrow = k + 1, ncol=n)
  X[k + 1, ] <- M
  
  X[2:k, ] <- replicate(n, sample(1:(M*(1-fixed_proportion)-1), k - 1))
  
  
  X <- apply(X, 2, sort)
  if (!is.null(fixed_proportion)) {
    X[k,] <- M - fixed_proportion*M
  }
  Y <- (X - X[c(k + 1, 1:k), ])[2:(k + 1), ]
  return(Y / M)
}

## Generate all matrices using provides pictures_W matrix. H sampled from simplex.
make_simulated_data<- function(M, N, K, noise_deviation = 0.1, noise_nature="proportional" ) {
  
  H <-  matrix(runif(K*N, min = 0, max = 10),nrow=K )  
  W <-  t(matrix(runif(K*M, min = 0, max = 10),nrow=K ) )
  
  V <- W %*% H
  
  rownames(V) <- c(paste0("feature", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  H
  true_W <- W
  colnames(true_H) <- colnames(V)
  rownames(true_W) <- rownames(V)
  if (noise_nature == "proportional") {
    noise_mask <- 1+ matrix(rnorm(length(V), sd = noise_deviation), nrow = nrow(V), ncol = ncol(V))
    noizy_V <-  V * noise_mask
  } else if (noise_nature == "additive") {
    noise_mask <-  matrix(rnorm(length(V), sd = noise_deviation), nrow = nrow(V), ncol = ncol(V))
    noizy_V <-  V + noise_mask
  }
  noizy_V[noizy_V< 0] <-0
  return (list(V=V, noizy_V= noizy_V, true_H=true_H, true_W=true_W))
}





## Generate all matrices using provides pictures_W matrix. H sampled from simplex.
make_picture_data <- function(M, N, K, pictures_W, noise_deviation = 0.1, noise_nature="proportional" ) {
  
  H <- linseed:::sampleFromSimplexUniformly(N, K, 10000)
  V <- pictures_W %*% H
  
  rownames(V) <- c(paste0("pixels", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  H
  true_W <- pictures_W
  colnames(true_H) <- colnames(V)
  rownames(true_W) <- rownames(V)
  if (noise_nature == "proportional") {
    noise_mask <- 1+ matrix(rnorm(length(V), sd = noise_deviation), nrow = nrow(V), ncol = ncol(V))
    noizy_V <-  V * noise_mask
  } else if (noise_nature == "additive") {
    noise_mask <-  matrix(rnorm(length(V), sd = noise_deviation), nrow = nrow(V), ncol = ncol(V))
    noizy_V <-  V + noise_mask
  }
  noizy_V[noizy_V< 0] <-0
  return (list(V=V, noizy_V= noizy_V, true_H=true_H, true_W=true_W))
  
}

## Generate all matrices using provides pictures_W matrix. H is random.
make_picture_data_random_mixtures <- function(M, N, K, pictures_W,noise_deviation = 0.2 ) {
  
  H <- matrix(runif(K*N),nrow=K)
  V <- pictures_W %*% H
  
  rownames(V) <- c(paste0("pixels", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  H
  true_W <- pictures_W
  colnames(true_H) <- colnames(V)
  rownames(true_W) <- rownames(V)
  noise_mask <- 1 + matrix(rnorm(length(V), sd = noise_deviation), nrow = nrow(V), ncol = ncol(V))
  noizy_V <-  V * noise_mask
  noizy_V[noizy_V< 0] <-0
  
  
  return (list(V=V, noizy_V= noizy_V, true_H=true_H, true_W=true_W))
  
}



#### End section


#### Logic methods ####

# Main realization of LinseedNMF method. Should be passed to nmf() package. Will train linseedv2 object
dualsimplex_nmf_algorithm <- function(x, seed, 
                                    start_with_hinge_W=1, start_with_hinge_H=1,
                                    LR_DECAY_STEPS = 7,
                                    RUNS_EACH_STEP = 500,
                                    PARAMETERS_INCREASE_STEPS = 2,
                                    lr_x =  0.1,
                                    lr_omega = 0.1
                                    ){
  factorization_rank <- nbasis(seed)
  
  
  dso <-  DualSimplexSolver$new()
  
  
  
  dso$set_data(x)
  dso$project(factorization_rank)
  dso$init_solution("random")
  
  lr_decay <- 0.1
  params_increase <- 100
  print(paste("current params", start_with_hinge_H, start_with_hinge_W))
  original_lambda_term <- start_with_hinge_H  #coef_hinge_H
  original_beta_term <- start_with_hinge_W #coef_hinge_W
  
  
  for  (lr_step in 1:LR_DECAY_STEPS) {
    lambda_term <-  original_lambda_term * lr_x * lr_x
    beta_term <- original_beta_term   * lr_omega * lr_omega
    RUNS_EACH_STEP <-  RUNS_EACH_STEP * 0.7
    for (incr_step in 1:PARAMETERS_INCREASE_STEPS) {
      # Main training method, you can just run this
      dso$optim_solution(RUNS_EACH_STEP, 
                         optim_config(coef_hinge_H = lambda_term,
                                      coef_hinge_W =beta_term,
                                      coef_der_X = lr_x, 
                                      coef_der_Omega = lr_omega,
                                      coef_pos_D_w = 0,
                                      coef_pos_D_h = 0
                         ))
      lambda_term <- lambda_term * params_increase
      beta_term <-  beta_term * params_increase
    }
    lr_x <- lr_x * lr_decay
    lr_omega <-  lr_omega * lr_decay
    
  }
  
  # R <-  dso$st$proj$meta$R
  # S <-  dso$st$proj$meta$S
  # X_solution <- dso$st$solution_proj$X
  # Omega_solution <- t(dso$st$solution_proj$Omega)
  # 
  # # Correct X and Omega to produce sum-to-one matrices
  # H_ss <-  X_solution %*% R
  # H_ss[H_ss < 0] <-  0
  # H_ss <-  H_ss / rowSums(H_ss) # should be row normalized
  # W_gs <-  t(S) %*% Omega_solution
  # W_gs[W_gs < 0] <-  0
  # W_gs <-  t(t(W_gs) / rowSums( t(W_gs))) # should be column normalized
  # 
  # X_corrected <- H_ss %*% t(R)  # recalculate  X for sum-to-one matrix
  # Omega_corrected <- S %*% W_gs # recalculate  Omega for sum-to-one matrix
  # 
  # dso$st$solution_proj$X <-  X_corrected
  # dso$st$solution_proj$Omega <-  Omega_corrected
  

  solution <- dso$finalize_solution() # this performs reverse sinkhorn procedure 
  sol_no_corr <-  dso$st$solution_no_corr  # factorization calculated from 1 step of sinkhorn

  # W matrix
  W <- sol_no_corr$Dv_inv_W_row
  W[W < 0] <-  0
  print(dim(W))
  basis(seed) <- W
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  
  H <- sol_no_corr$H_row
  print(dim(H))
  H[H < 0] <-  0
  coef(seed) <- H
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  # return updated data
  seed$obj <- dso
  
  print("Getting distance solution")


  distance <- nmfDistance(objective(seed))(seed, x)
  print(paste("objective is ", distance) )
  if (is.na(distance)) {
    print("Seems like objective is NA. Could happen randomly. Reruning the training")
    seed <- dualsimplex_nmf_algorithm(x, seed, start_with_hinge_W=start_with_hinge_W, start_with_hinge_H=start_with_hinge_H)
  }
  print("return from method")
  return(seed)
}

# HALSacc algorithm
# "Accelerated hierarchical alternating least squares NMF. For a refer-
# ence to the method, see N. Gillis, Nonnegative matrix factorization:
# complexity, algorithms and applications [Section 4.2, Algo. 6], PhD
# thesis, Université catholique de Louvain, February 2011"
halsacc_nmf_algorithm <- function(x, seed){
  print("New run")
  
  factorization_rank <- nbasis(seed)
  
  
  X <- x
  bgImageTensor <- array(0,dim=dim(X))
  selectVect <- array(1,dim=dim(X))
  nmfInput <- NULL
  nmfInput$numRows <- nrow(X)
  nmfInput$numCols <- ncol(X)
  nmfInput$numSlices <- 1
  nmfInput$bgImageTensor <- bgImageTensor
  nmfInput$selectVect <- selectVect
  
  
  solution <- oneLevelNMF(X, rank=factorization_rank, method="HALSacc", nruns=1)

  # W matrix
  W <- solution@W
  print(dim(W))
  W[W < 0] <-  0
  basis(seed) <- W
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  
  H <- solution@H
  H[H < 0] <-  0
  coef(seed) <- H
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  # return updated data
  seed$obj <- solution
  
  
  distance <- nmfDistance(objective(seed))(seed, x)
  print(paste("objective is ", distance) )
  return(seed)
}


# PGNMF algorithm
# "NMF by alternating non-negative least squares using projected gradi-
# ents. For a reference to the method, see C.-J. Lin, "Projected Gradient
# Methods for Non-negative Matrix Factorization, Neural computation
# 19.10 (2007): 2756-2779."


pgnmf_nmf_algorithm <- function(x, seed){
  print("New run")
  
  factorization_rank <- nbasis(seed)
  
  
  X <- x
  bgImageTensor <- array(0,dim=dim(X))
  selectVect <- array(1,dim=dim(X))
  nmfInput <- NULL
  nmfInput$numRows <- nrow(X)
  nmfInput$numCols <- ncol(X)
  nmfInput$numSlices <- 1
  nmfInput$bgImageTensor <- bgImageTensor
  nmfInput$selectVect <- selectVect
  
  
  solution <- oneLevelNMF(X, rank=factorization_rank, method="PGNMF", nruns=1)
  
  # W matrix
  W <- solution@W
  print(dim(W))
  W[W < 0] <-  0
  basis(seed) <- W
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  
  H <- solution@H
  H[H < 0] <-  0
  coef(seed) <- H
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  # return updated data
  seed$obj <- solution
  
  
  distance <- nmfDistance(objective(seed))(seed, x)
  print(paste("objective is ", distance) )
  return(seed)
}

# ALS
#DeBruine, ZJ, Melcher, K, and Triche, TJ. (2021). "High-performance non-negative matrix factorization for large single-cell data." BioRXiv. 
als_nmf_algorithm <- function(x, seed){
  factorization_rank <- nbasis(seed)
  
  solution <- RcppML::nmf(Matrix(x), k = factorization_rank)
  
  # W matrix
  W <- solution$w
  print(dim(W))
  W[W < 0] <-  0
  basis(seed) <- W
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  
  H <- solution$h
  H[H < 0] <-  0
  coef(seed) <- H
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  # return updated data
  seed$obj <- solution
  
  
  distance <- nmfDistance(objective(seed))(seed, x)
  print(paste("objective is ", distance) )
  return(seed)
}


#### End section




#### Plot methods ####

# Metric to compare models
pearson_correlation_function <- function(target, w){
  1-cor(w, target, method ="pearson")^2
       # + eta * sum(w^2) 
       # + beta * sum( colSums( h )^2 )
}

# Metric to compare models
cosine_distance_function <- function(target, w){
  1-lsa::cosine(w, target)
  # + eta * sum(w^2) 
  # + beta * sum( colSums( h )^2 )
}

spearman_correlation_function <- function(target, w){
  1-cor(w, target, method ="spearman")^2
  # + eta * sum(w^2) 
  # + beta * sum( colSums( h )^2 )
}

# Metric to compare models
RMSE_function <- function(target, w){
  
  sqrt((1/(length(w)) * sum( (target - w )^2 )) 
       # + eta * sum(w^2) 
       # + beta * sum( colSums( h )^2 )
  )
}

# Metric to compare models, will scale matrices (since methods get different scales) 
normalized_RMSE_function <- function(target, w){
  w[w < 0] <-  0
  target <- (target - min(target)) / (max(target) - min(target))
  w <- (w - min(w)) / (max(w) - min(w))
  
  sqrt((1/(length(w)) * sum((target - w )^2 )) 
       # + eta * sum(w^2) 
       # + beta * sum( colSums( h )^2 )
  )
}

# Will calculate RMSE distances for each component (column)
get_metric_distances_by_column <- function(original_matrix, curr_matrix, metric="rmse"){
  col_num <-  dim(original_matrix)[[2]]
  distances <-  sapply(c(1:col_num), function(original_component_index){
    candidates <- lapply(c(1:col_num), function(current_component_index){
      # print(paste(original_component_index,current_component_index))
      #return (RMSE_function(original_matrix[,original_component_index], curr_matrix[,current_component_index]))
      result <-  switch(metric, 
                        "rmse" = RMSE_function(original_matrix[,original_component_index], curr_matrix[,current_component_index]),
                        "pearson" = pearson_correlation_function(original_matrix[,original_component_index], curr_matrix[,current_component_index]),
                        "spearman" = spearman_correlation_function(original_matrix[,original_component_index], curr_matrix[,current_component_index]),
                        "cosine" = cosine_distance_function(original_matrix[,original_component_index], curr_matrix[,current_component_index]))
      return(result)
    })
    # print(unlist(candidates))
    return (min(unlist(candidates)))
  })
  return (setNames(data.table(c(1:col_num), unlist(distances)), c("component", metric)))
}

# Will calculate RMSE distances for each component (row)
get_metric_distances_by_row <- function(original_matrix, curr_matrix, metric="rmse"){
  current_K <-  dim(original_matrix)[[1]]
  distances <-  sapply(c(1:current_K), function(original_component_index){
    
    candidates <- lapply(c(1:current_K), function(current_component_index){
      # print(paste(original_component_index,current_component_index))
      result <-  switch(metric, 
                        "rmse" = RMSE_function(
                          original_matrix[original_component_index,], curr_matrix[current_component_index,]),
                        "pearson" = pearson_correlation_function(
                          original_matrix[original_component_index,], curr_matrix[current_component_index,]),
                        "spearman" = spearman_correlation_function(
                          original_matrix[original_component_index,], curr_matrix[current_component_index,]),
                        "cosine" = cosine_distance_function(
                          original_matrix[original_component_index,], curr_matrix[current_component_index,]))
      
      
      return (RMSE_function(original_matrix[original_component_index,], curr_matrix[current_component_index,]))
    })
    # print(unlist(candidates))
    return (min(unlist(candidates)))
  })
  return (setNames(data.table(c(1:current_K), unlist(distances)), c("component", metric)))
}


# Will calculate RMSE distances for each component (row)
get_D_distances_on_diagonal <- function(original_matrix, curr_matrix){
  current_K <-  dim(original_matrix)[[1]]
  distances <-  sapply(c(1:current_K), function(original_component_index){
    candidates <- lapply(c(1:current_K), function(current_component_index){
      # print(paste(original_component_index,current_component_index))
      difference <- abs(original_matrix[original_component_index, original_component_index] - curr_matrix[current_component_index,current_component_index])
      percentage_change <-  difference / original_matrix[original_component_index, original_component_index] 
      
      
      distance_value <- percentage_change
      return (distance_value)
    })
    # print(unlist(candidates))
    return (min(unlist(candidates)))
  })
  return (data.table(component=c(1:current_K), diff=unlist(distances)))
}




## Will calculate RMSE distances for each method and each run. W matrix.
get_W_metric_values <- function(res_methods, data_used, metric="rmse") {
  
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <-  data_used[[current_run]]
      curr_W <- run_result$W
      true_W <- curr_data$true_W
      distances_table <- get_metric_distances_by_column(true_W, curr_W,metric=metric)
      distances_table$run <- current_run
      return(distances_table)
    })
    method_df <- rbindlist(per_run_tables)
    
    method_df$method <- method_name
    return(method_df)
    
  })
  
  total_df <- rbindlist(per_method_tables)
  return(total_df)
}


## Will calculate RMSE distances for each method and each run. H matrix.
get_H_metric_values <- function(res_methods, data_used, metric="rmse") {
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <- data_used[[current_run]]
      curr_H <- run_result$H
      curr_W <- run_result$W

      true_H <- curr_data$true_H
      true_W <- curr_data$true_W

      
      distances_table <- get_metric_distances_by_row(true_H, curr_H, metric=metric)
      distances_table$run <- current_run
      return(distances_table)
    })
    method_df <- rbindlist(per_run_tables)
    
    method_df$method <- method_name
    return(method_df)
    
  })
  
  total_df <- rbindlist(per_method_tables)
  return(total_df)
}


get_V_metric_values <- function(res_methods, data_used, metric="rmse") {
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <-  data_used[[current_run]]
      curr_H <- run_result$H
      curr_W <- run_result$W
      curr_D <- run_result$D_inv
      curr_H[curr_H < 0] <-  0
      curr_W[curr_W < 0] <-  0
      
      true_H <- curr_data$true_H
      true_W <- curr_data$true_W
      true_D <- curr_data$D_inv
      
      curr_V <- curr_W %*% curr_D %*% curr_H
      true_V <- true_W %*% true_D %*% true_H
      distance <-  switch(metric, 
                          "rmse" = RMSE_function( true_V, curr_V),
                          "pearson" = mean(diag(pearson_correlation_function(true_V, curr_V))),
                          "spearman" = mean(diag(spearman_correlation_function(true_V, curr_V))))
      distances_table <-  setNames(data.table(distance), c(metric))
      distances_table$run <- current_run
      return(distances_table)
    })
    method_df <- rbindlist(per_run_tables)
    
    method_df$method <- method_name
    return(method_df)
    
  })
  
  total_df <- rbindlist(per_method_tables)
  return(total_df)
}




# get_V_rmse_values <- function(res_methods, data_used) {
#   
#   
#   per_method_tables <- lapply(names(res_methods), function (method_name){
#     curr_method_result <- res_methods[[method_name]]
#     per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
#       run_result <- curr_method_result[[current_run]]
#       curr_data <-  data_used[[current_run]]
#       curr_H <- run_result$H
#       curr_W <- run_result$W
#       curr_D <- run_result$D_inv
#       curr_H[curr_H < 0] <-  0
#       curr_W[curr_W < 0] <-  0
#       
#       true_H <- curr_data$true_H
#       true_W <- curr_data$true_W
#       true_D <- curr_data$D_inv
#       
#       distance <- RMSE_function(true_W %*% true_D %*% true_H, curr_W %*% curr_D %*% curr_H)
#       
#       distances_table <-  data.table(rmse=distance)
#       distances_table$run <- current_run
#       return(distances_table)
#     })
#     method_df <- rbindlist(per_run_tables)
#     
#     method_df$method <- method_name
#     return(method_df)
#     
#   })
#   
#   total_df <- rbindlist(per_method_tables)
#   return(total_df)
# }



## Will calculate RMSE distances for each method and each run. H matrix.
get_D_difference_values <- function(res_methods, data_used) {
  
  
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <- data_used[[current_run]]

      curr_D <- run_result$D_inv
      true_D <- curr_data$D_inv
      distances_table <- get_D_distances_on_diagonal(true_D, curr_D)
      distances_table$run <- current_run
      return(distances_table)
    })
    method_df <- rbindlist(per_run_tables)
    
    method_df$method <- method_name
    return(method_df)
    
  })
  
  total_df <- rbindlist(per_method_tables)
  return(total_df)
}



get_significance_values <-  function(distance_matrix, comparisons, vertical_gap = 0.03, metric='rmse') {
  stat_test_results <-  lapply(c(1:length(comparisons)), function(i){
    comparison <- comparisons[[i]]
    data <- distance_matrix[distance_matrix$method %in% comparison]
    p.values <- sapply(split(data, data$component), 
                       function(x){
                         wilcox.test(x[[metric]]~method, data=x)$p.value})
    labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","n.s."))
    
    y.values <- sapply(split(distance_matrix, distance_matrix$component), function(x){max(x[[metric]])})+ i*vertical_gap
    #y.values <- sapply(split(distance_matrix, distance_matrix$component), function(x){max(x$rmse)})+ i*0.2 # log
    return (list(labels=labels, y.values=y.values))
  })
  return (stat_test_results)
}





# Move energy from H to W to make everything comparable
correct_solution_to_normalize_H_and_W <- function(solution_W, solution_H) {
  D_w_inv <- diag(colSums(solution_W))
  D_w <- diag(1/colSums(solution_W))
  D_h_inv <- diag(rowSums(solution_H))
  D_h <- diag(1/rowSums(solution_H))
  
  new_W <-  solution_W %*% D_w  # should get col normalized W
  new_H <-  D_h %*% solution_H # should get row normalized H
  D_inv <-  D_w_inv %*% D_h_inv
  return (list(new_W=new_W, new_H=new_H, D_inv=D_inv))
}


correct_method_matrices_to_make_them_comparable <- function(res_methods) {
  
  corrected_methods_results <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_result <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_H <- run_result@fit@H
      curr_W <- run_result@fit@W
      curr_H[curr_H < 0] <-  0
      curr_W[curr_W < 0] <-  0
      correction_result <- correct_solution_to_normalize_H_and_W(curr_W, curr_H)
      curr_W <- correction_result$new_W
      curr_H <- correction_result$new_H
      run_result@fit@H <-  curr_H
      run_result@fit@H <-  curr_W
      return(list(H=curr_H, W=curr_W, D_inv=correction_result$D_inv))
    })
    return(per_run_result)
  })
  names(corrected_methods_results) <-  names(res_methods)
  return(corrected_methods_results)
  
}

correct_data_to_make_data_comparable <-  function(data_used) {
  per_run_data <- lapply(c(1:length(data_used)), function(current_run) {
    curr_data <-  data_used[[current_run]]
    true_H <- curr_data$true_H
    true_W <- curr_data$true_W
    correction_result <- correct_solution_to_normalize_H_and_W(true_W, true_H)
    return (list(true_W=correction_result$new_W,  true_H=correction_result$new_H, V=curr_data$V , D_inv=correction_result$D_inv))})
  return(per_run_data)
}



plot_metrics_box <- function(res_methods, data_used, folder_to_save_results, title="Matrices", reference_method="DualSimplex", 
                             W_metric="rmse", H_metric="rmse", V_metric="rmse",
                             W_signif_step=0.1, H_signif_step=0.1, V_signif_step=0.1,
                             W_signif_margin_top=0.01, H_signif_margin_top=0.01,  V_signif_margin_top=0.01,
                             H_limit=NULL, W_limit=NULL, V_limit=NULL, log_scale_W = F, log_scale_V=T, log_scale_H=F) {
  #' Get all boxplots
  #'
  #' Creates a bunch of boxplots in specified folder
  #' @param res_methods fit objects obtained from nmf experiment method ($fit_objects). this if per-method-per-run list of results of nmf() function
  #' @param data_used true data used for each run. Expected list of num_runs elements which contains. $true_H, $true_W
  #' @param  title just common title base string for each plot. 
  #' @param reference_method reference method name to calculate significance values. 
  #' @param W_metric  which metric to use for W plots. can be 'rmse', 'cosine' , 'pearson', 'spearman' will be calulated per column
  #' @param H_metric  which metric to use for H plots. can be 'rmse', 'cosine' , 'pearson', 'spearman' will be calculated per row
  #' @param V_metric  which metric to use for V plots. can be 'rmse' , 'pearson', 'spearman'  values will be the mean for columns
  corrected_res_methods <- correct_method_matrices_to_make_them_comparable(res_methods)
  corrected_data <-  correct_data_to_make_data_comparable(data_used)
  
  
  W_matrix_distances <-  get_W_metric_values(res_methods = corrected_res_methods, data_used = corrected_data, metric=W_metric)
  #W_matrix_distances[method == 'Linseedv2 NMF']$method <- 'DualSimplex'
  H_matrix_distances <-  get_H_metric_values(res_methods = corrected_res_methods, data_used = corrected_data, metric=H_metric)
  #H_matrix_distances[method == 'Linseedv2 NMF']$method <- 'DualSimplex'
  V_matrix_distances <-  get_V_metric_values(res_methods =corrected_res_methods, data_used = corrected_data, metric= V_metric )
  D_matrix_distances <-  get_D_difference_values(res_methods = corrected_res_methods, data_used = corrected_data)
  
  per_metric_titles <- c()
  per_metric_titles["rmse"] <-  expression("RMSE")
  per_metric_titles["pearson"] <-  expression(1 - pearson^2 )
  per_metric_titles["spearman"] <-   expression(1 - spearman^2 )
  
  

  
  all_methods <- sort(names(res_methods))
  comparisons  <- lapply(all_methods[all_methods != reference_method], function(x) {c(reference_method, x)})
  
  W_matrix_distances$method <- relevel(factor(W_matrix_distances$method), ref=reference_method)
  H_matrix_distances$method <- relevel(factor(H_matrix_distances$method), ref=reference_method)
  V_matrix_distances$method <- relevel(factor(V_matrix_distances$method), ref=reference_method)
  D_matrix_distances$method <- relevel(factor(D_matrix_distances$method), ref=reference_method)
  
  median_df_W <- W_matrix_distances %>%  
    group_by(method, component) %>%  
    summarize(average = median(.data[[W_metric]])) %>% 
    ungroup() 
  median_df_H <- H_matrix_distances %>%  
    group_by(method, component) %>%  
    summarize(average = median(.data[[H_metric]])) %>% 
    ungroup() 
  median_df_D <- D_matrix_distances %>%  
    group_by(method, component) %>%  
    summarize(average = median(diff)) %>% 
    ungroup()
  
  
  gw <- ggplot(W_matrix_distances, aes(x=factor(component), y=.data[[W_metric]], fill=method)) +
    
    geom_boxplot(position = position_dodge2(width = 0.1)) +  
    geom_line(data=median_df_W, aes(x=component, y=average, color=method), 
              position = position_jitterdodge(jitter.width = 0.1),  size=1.2, alpha=0.4) +  
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size=0.4) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25,base_family = 'sans') +
    theme(plot.title = element_text(face='bold', size=21)) +
    theme(axis.title.x  = element_text(size=20)) +
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=14)) +
    theme(axis.text.y  = element_text(size=14)) +
    theme(legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 5, unit = "pt")) +
    labs(title = "Distance to W",
         x="Component", 
         y= per_metric_titles[W_metric] #subtitle = "LOL", x="LOLX", y="LOLY"
    )
  stat_test_results <- get_significance_values(W_matrix_distances, comparisons, metric=W_metric, vertical_gap = 0.01 )

  # for (comp_index in 1:(length(stat_test_results) -1 )) {
  #   print(stat_test_results[[comp_index]]$y.values)
  # gw <- gw + geom_signif(y_position = stat_test_results[[comp_index]]$y.values, 
  #                        xmin = unique(W_matrix_distances$component) - (3*0.5/6),
  #                        xmax = unique(W_matrix_distances$component) -(0.5/6)+ 0.04*comp_index ,
  #                        annotations = stat_test_results[[comp_index]]$labels,
  #                        tip_length = 0)
  # }
  gw <-  gw 
  
  summarized_W <- ggplot(W_matrix_distances, aes(x=method, y=.data[[W_metric]], fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1)) +
    stat_compare_means(comparisons=comparisons, label = "p.signif", tip.length = 0,  ref.group = "DualSimplex", size=6) +
    
    #geom_signif(comparisons=comparisons, 
               # map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
               # margin_top = W_signif_margin_top,  textsize = 6, #tip_length = 0.1,# vjust=0.1,
              #  step_increase = W_signif_step)+
    geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25,base_family = 'sans') +
    labs(title = "Distance to original columns", y=per_metric_titles[W_metric], x=NULL) + 
    theme(plot.title = element_text(face='bold', size=21)) + 
    theme(axis.title.x  = element_text(size=20)) +  
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=18)) +  
    theme(axis.text.y  = element_text(size=14)) + 
    theme(legend.position = "none", legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1)) + 
    theme(plot.margin = margin(t = 3, r = 0, b = 3, l = 5, unit = "pt"))
  if (is.null(W_limit)) {
    #W_limit <- 1.2 * max(W_matrix_distances[[W_metric]])
  }

  if (log_scale_W) {
    print("scale W")
    summarized_W <-  summarized_W + scale_y_continuous(trans='log10', limits = c(NA, NA))
  }
  else {
   summarized_W <-  summarized_W + ylim(0, NA)
  }

  

  
  per_component_w_plot_list <- lapply(unique(W_matrix_distances$component), function(component_num) {
    curr_W_distances <-  W_matrix_distances[component == component_num]
    curr_plot <- ggplot(curr_W_distances, aes(x=method, y=.data[[W_metric]], fill=method)) +
      geom_boxplot(position = position_dodge2(width = 0.1)) +
      stat_compare_means(label = "p.signif",  ref.group = "DualSimplex", size=6) +
      # geom_signif(comparisons=comparisons, 
      #             map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
      #             margin_top = W_signif_margin_top, tip_length = 0, textsize = 6,# vjust=0.1,
      #             step_increase = W_signif_step)+
      geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
      scale_fill_brewer(palette = "Set1")+
      scale_color_brewer(palette = "Set1") +
      theme_classic(base_size=25,base_family = 'sans') +
      labs(title = paste0("Column ", component_num), y=per_metric_titles[W_metric], x=NULL #subtitle = "LOL", x="LOLX", y="LOLY"
      ) + 
      theme(plot.title = element_text(size=20)) +  
      theme(axis.title.x  = element_text(size=20)) +  
      theme(axis.title.y  = element_text(size=20))+
      theme(axis.text.x  = element_text(size=18)) +  
      theme(axis.text.y  = element_text(size=14)) + 
      theme(legend.position = "none", legend.text=element_text(size=16))+
      theme(legend.title=element_text(size=16)) +
      theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1)) + 
      theme(plot.margin = margin(t = 3, r = 0, b = 3, l = 5, unit = "pt"))
      

    if (is.null(W_limit)) {
      W_limit <- 1.7 * max(curr_W_distances[[W_metric]])
    }
    curr_plot <-  curr_plot + ylim(0, W_limit)
    return (curr_plot)
    
  })
  
  # H matrix
  gh <- ggplot(H_matrix_distances, aes(x=factor(component), y=.data[[H_metric]], fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1)) +  
    geom_line(data=median_df_H, aes(x=component, y=average, color=method), 
              position = position_jitterdodge(jitter.width = 0.1),  size=1.2, alpha=0.4) +  
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size=0.4) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25,base_family = 'sans') +
    theme(plot.title = element_text(face='bold', size=21)) +
    theme(axis.title.x  = element_text(size=20)) +
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=14)) +
    theme(axis.text.y  = element_text(size=14)) +
    theme( legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 5, unit = "pt")) +
    labs(title =  "Distance to H" ,x="Component", y=per_metric_titles[H_metric]  #subtitle = "LOL", x="LOLX", y="LOLY"
    )  
  stat_test_results <- get_significance_values(H_matrix_distances, comparisons, vertical_gap = 0.001, metric = H_metric)
  for (comp_index in 1:(length(stat_test_results) -1 )) {
  gh <- gh + geom_signif(y_position = stat_test_results[[comp_index]]$y.values, 
                         xmin = unique(H_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(H_matrix_distances$component)-(0.5/4)+ 0.04 ,
                         annotations = stat_test_results[[comp_index]]$labels,
                         tip_length = 0)
  }
  # gh <- gh + geom_signif(y_position = stat_test_results[[2]]$y.values, 
  #                        xmin = unique(H_matrix_distances$component)- (2*0.5/4),
  #                        xmax = unique(H_matrix_distances$component)+(0.5/4) - 0.02 ,
  #                        annotations = stat_test_results[[2]]$labels,
  #                        tip_length = 0)
  # 
  # gh <- gh + geom_signif(y_position = stat_test_results[[3]]$y.values, 
  #                        xmin = unique(H_matrix_distances$component)- (2*0.5/4),
  #                        xmax = unique(H_matrix_distances$component)+ (2*0.5/4)+ 0.04 ,
  #                        annotations = stat_test_results[[3]]$labels,
  #                        tip_length = 0)
  
  per_component_h_plot_list <- lapply(unique(H_matrix_distances$component), function(component_num) {
    curr_H_distances <-  H_matrix_distances[component == component_num]
    curr_plot <- ggplot(curr_H_distances, aes(x=method, y=.data[[H_metric]], fill=method)) +
      geom_boxplot(position = position_dodge2(width = 0.1)) +
      stat_compare_means(label = "p.signif",  ref.group = "DualSimplex", size=6) +
      # geom_signif(comparisons=comparisons, 
      #             map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
      #             margin_top = H_signif_margin_top, tip_length = 0, textsize = 6,# vjust=0.1,
      #             step_increase = H_signif_step)+
      geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
      scale_fill_brewer(palette = "Set1")+
      scale_color_brewer(palette = "Set1") +
      theme_classic(base_size=25,base_family = 'sans') +
      labs(title = paste0("Row ", component_num),x=NULL, y=per_metric_titles[H_metric] #subtitle = "LOL", x="LOLX", y="LOLY"
      ) + 
      theme(plot.title = element_text(size=20)) +  
      theme(axis.title.x  = element_text(size=20)) +  
      theme(axis.title.y  = element_text(size=20))+
      theme(axis.text.x  = element_text(size=18)) +  
      theme(axis.text.y  = element_text(size=14)) + 
      theme(legend.position = "none", legend.text=element_text(size=16))+
      theme(legend.title=element_text(size=16)) +
      theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1)) + 
      theme(plot.margin = margin(t = 3, r = 0, b = 3, l = 5, unit = "pt"))
    if (is.null(H_limit)) {
      H_limit <- 1.5 * max(curr_H_distances[[H_metric]])
    }
    curr_plot <-  curr_plot + ylim(0, H_limit)
    return (curr_plot)
    
  })
  
  
  summarized_H <- ggplot(H_matrix_distances, aes(x=method, y=.data[[H_metric]], fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1)) +
    stat_compare_means(comparisons=comparisons,label = "p.signif", tip.length = 0, ref.group = "DualSimplex", size=6) +
    # geom_signif(comparisons=comparisons, 
    #             map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
    #             margin_top = H_signif_margin_top, tip_length = 0, textsize = 6,# vjust=0.1,
    #             step_increase = H_signif_step)+
    geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25,base_family = 'sans') +
    labs(title = "Distance to original rows",x=NULL, y=per_metric_titles[H_metric] #subtitle = "LOL", x="LOLX", y="LOLY"
    ) + 
    theme(plot.title = element_text(size=20)) +  
    theme(plot.title = element_text(face='bold', size=21)) + 
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=18)) +  
    theme(axis.text.y  = element_text(size=14)) + 
    theme(legend.position = "none", legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1)) + 
    theme(plot.margin = margin(t = 3, r = 0, b = 3, l = 5, unit = "pt"))
  if (is.null(H_limit)) {
  #  H_limit <- 1.5 * max(H_matrix_distances[[H_metric]])
  }
  if (log_scale_H) {
    print("scale H")
    summarized_H <-  summarized_H + scale_y_continuous(trans='log10', limits = c(NA, NA))
  }
  else {
   summarized_H <-  summarized_H + ylim(0, NA)
  }
  
  
  
  
  gv <- ggplot(V_matrix_distances, aes(x=method, y=.data[[V_metric]], fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1))+
    stat_compare_means(comparisons=comparisons,label = "p.signif", tip.length = 0, ref.group = "DualSimplex", size=6) +
    # geom_signif(comparisons=comparisons, 
    #             map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
    #             textsize=6,
    #             margin_top = V_signif_margin_top, tip_length = 0,
    #             step_increase = V_signif_step,)+
    scale_fill_brewer(palette = "Set1")+
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
    labs(title = "Distance to original columns",x=NULL, y=per_metric_titles[V_metric]  #subtitle = "LOL", x="LOLX", y="LOLY"
    )+ 
    theme_classic(base_size=25, base_family = 'sans') +
    theme(plot.title = element_text(face='bold', size=21)) +  
    theme(axis.title.x  = element_text(size=20)) +  
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=18)) +  
    theme(axis.text.y  = element_text(size=14)) + 
    theme(legend.position = "none", legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16))+
    theme(axis.text.x  = element_text(angle=45, vjust=1, hjust=1)) + 
    theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 5, unit = "pt")) 
  
    if (is.null(V_limit)) {
    #  V_limit <- 1.5 * max(V_matrix_distances[[V_metric]])
    }
  if (log_scale_V) {
    print("scale V")
   # print(paste("V limit", V_limit))
    gv <-  gv + scale_y_continuous(trans='log10', limits = c(NA, NA))
  } else {
   gv <-  gv + ylim(0, NA)
  }
  

  
  # Plot for W matrix
  gd <- ggplot(D_matrix_distances, aes(x=factor(component), y=diff, fill=method)) +
    
    geom_boxplot(position = position_dodge2(width = 0.1)) +  
    geom_line(data=median_df_D, aes(x=component, y=average, color=method), 
              position = position_jitterdodge(jitter.width = 0.1),  size=1.2, alpha=0.4) +  
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size=0.4) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25,base_family = 'sans') +
    theme(plot.title = element_text(face='bold', size=21)) +
    theme(axis.title.x  = element_text(size=20)) +
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=14)) +
    theme(axis.text.y  = element_text(size=14)) +
    theme(legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(plot.margin = margin(t = 2, r = 0, b = 2, l = 5, unit = "pt")) +
    
    labs(title = "Distance to D",
         x="Component", 
         y="percentage_change" #subtitle = "LOL", x="LOLX", y="LOLY"
    )
  
  stat_test_results <- get_significance_values(D_matrix_distances, comparisons, metric='diff')
  for (comp_index in 1:(length(stat_test_results) -1 )) {
    gd <- gd + geom_signif(y_position = stat_test_results[[comp_index]]$y.values, 
                           xmin = unique(W_matrix_distances$component)- (3*0.5/6),
                           xmax = unique(W_matrix_distances$component)-(0.5/6)+ 0.04*comp_index ,
                           annotations = stat_test_results[[comp_index]]$labels,
                           tip_length = 0)
  }
  
  summarized_plot <- ggarrange(gv, summarized_W, summarized_H, nrow = 1, common.legend = TRUE, legend='right')
  
  

  
  t <- ggarrange(gw, gh, gd, gv, nrow = 1, common.legend = TRUE, legend='right')
  filename <- paste0(folder_to_save_results,  title, ".svg")
  ggsave(file=filename, plot=t, width=18, height=5, device=svglite)
  
  filename <- paste0(folder_to_save_results, "V_", title, ".svg")
  ggsave(filename, gv, width=3.5, height=5, device=svglite)
  filename <- paste0(folder_to_save_results, "H_", title, ".svg")
  ggsave(filename, gh, width=8, height=5, device=svglite)
  filename <- paste0(folder_to_save_results, "W_", title, ".svg")
  ggsave(filename, gw, width=8, height=5, device=svglite)
  filename <- paste0(folder_to_save_results, "W_summary_", title, ".svg")
  ggsave(filename, summarized_W, width=3.5, height=5, device=svglite)
  filename <- paste0(folder_to_save_results, "H_summary_", title, ".svg")
  ggsave(filename, summarized_H, width=3.5, height=5, device=svglite)
  
  # for (comp_ind in 1:length(per_component_w_plot_list)) {
  #   filename <- paste0(folder_to_save_results, "W_", comp_ind, title, ".svg")
  #   ggsave(filename, per_component_w_plot_list[[comp_ind]], width=7.5, height=5, device=svglite)
  #   filename <- paste0(folder_to_save_results,"H_", comp_ind, title, ".svg")
  #   ggsave(filename, per_component_h_plot_list[[comp_ind]], width=7.5, height=5, device=svglite)
  # }
  per_comp_common_w <- ggarrange(plotlist=per_component_w_plot_list, nrow = 1, common.legend = TRUE, legend='none')
  per_comp_common_w <- annotate_figure(per_comp_common_w, top=text_grob(
    "Distance to original matrix W", face="bold", size= 21))
  
  filename <- paste0(folder_to_save_results,"W_common_", title, ".svg")
  ggsave(filename, per_comp_common_w, width=13, height=5, device=svglite)
  
  per_comp_common_h <- ggarrange(plotlist=per_component_h_plot_list, nrow = 1, common.legend = TRUE, legend='none')
  per_comp_common_h <- annotate_figure(per_comp_common_h, top=text_grob(
    "Distance to original matrix H", face="bold", size= 21))
  filename <- paste0(folder_to_save_results,"H_common_", title, ".svg")
  ggsave(filename, per_comp_common_h, width=14, height=5, device=svglite)
  
  
  return (list(W_matrix_distances=W_matrix_distances, 
               H_matrix_distances=H_matrix_distances,
               V_matrix_distances=V_matrix_distances, 
               D_matrix_distances=D_matrix_distances,
               plot=t,
               V_plot=gv,
               H_plot=gh,
               W_plot=gw,
               per_comp_plot_w=per_comp_common_w,
               per_comp_plot_h=per_comp_common_h,
               summarized_W = summarized_W,
               summarized_H=summarized_H,
               summarized_plot=summarized_plot
  ))
}






#Plot  W H V matrices
plot_w_h_v <- function(W, H, V, title="", Vrow_name ="words", Vcol_name="documents", component_name="topics (K)") {
  grid.newpage()
  
  pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2,heights = unit(c(0.1, 0.9),"null"))))
  grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  ht_W = Heatmap(W,name="W",column_title=component_name,row_title=Vrow_name, cluster_rows=FALSE, cluster_columns = FALSE, show_column_names=FALSE, show_row_names=FALSE)
  draw(ht_W, newpage = FALSE, )
  
  upViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  ht_H = Heatmap(H, name="H",column_title=Vcol_name, row_title=component_name, cluster_rows=FALSE,cluster_columns = FALSE,show_column_names=FALSE, show_row_names=FALSE)
  draw(ht_H, newpage = FALSE, )
  upViewport()
  
  Heatmap(as.matrix(V), name="V", column_title=Vcol_name, row_title=Vrow_name, cluster_rows=FALSE,cluster_columns = FALSE, show_column_names=FALSE, show_row_names=FALSE, )
  
}

# Plot triangles(X and Omega matrices) for linseed object.
check_triangles <- function(res_methods) {
  linseed_results <- res_methods[["DualSimplex"]]
  sapply(linseed_results, function(single_result) {
    lo_obj <- single_result$obj
    lo_obj$plot_projected("zero_distance", "zero_distance",with_solution=TRUE, use_dims = list(2:3))
  })
  
}
tf <- function(m) t(m)[, nrow(m):1]
# Plot all result pictures for all methods 
plot_pictures <- function(W_matrix, picture_title="title", folder_to_save_results) {
  filename <- paste0(folder_to_save_results, picture_title, ".jpeg")
  jpeg(filename, width = 1850, height = 400)
  
  number_of_pictures <- dim(W_matrix)[[2]]
  size <- as.integer(sqrt(dim(W_matrix)[[1]]))
  par(mfrow = c(1,number_of_pictures),mar=c(0.5, 0.5, 0.5,0.5), oma=c(0.5, 0.5, 0.5,0.5)) # all sides have 3 lines of space
  res <- lapply(c(1:number_of_pictures), function(picture_index){
    matrix_to_draw <- tf(matrix(W_matrix[, picture_index], nrow = size, ncol = size, byrow = TRUE))
    image(matrix_to_draw,
          col=grey(seq(0, 1, length = 64)), axes = FALSE)
    return(matrix_to_draw)
    
  })
  #title(main = picture_title, font.main = 4)
  #title(main= picture_title, outer = TRUE)
  dev.off()
  return(res)
  
}


plot_pictures_for_result <- function(method_results, model_index, folder_to_save_results, individual_pictures=F) {
  res <- lapply(method_results, function (method_multiple_results){
    method_result <- method_multiple_results[[model_index]]
    curr_W <- method_result@fit@W
    if (individual_pictures) {
      plot_single_pictures(curr_W, method_result@method, folder_to_save_results)
    } else {

    plot_pictures(curr_W, method_result@method, folder_to_save_results)
    }
    #par(oma=c(0,0,0,0)) # all sides have 3 lines of space
  })
  return(res)
  
}



plot_single_pictures <- function(W_matrix, picture_title="title", folder_to_save_results) {
  number_of_pictures <- dim(W_matrix)[[2]]
  size <- as.integer(sqrt(dim(W_matrix)[[1]]))
 
  res <- lapply(c(1:number_of_pictures), function(picture_index){
    filename <- paste0(folder_to_save_results, picture_title,picture_index,  ".jpeg")
    jpeg(filename, width = 400, height = 400)
    par(mar=c(0.5, 0.5, 0.5,0.5), oma=c(0.5, 0.5, 0.5,0.5)) # all sides have 3 lines of space
    matrix_to_draw <- tf(matrix(W_matrix[, picture_index], nrow = size, ncol = size, byrow = TRUE))
    image(matrix_to_draw,
          col=grey(seq(0, 1, length = 64)), axes = FALSE)
    dev.off()
    return(matrix_to_draw)
    
  })
  #title(main = picture_title, font.main = 4)
  #title(main= picture_title, outer = TRUE)

  return(res)
  
}



#### End section




#### Just left this for history (two metric which are commonly used for NMF) ####

Euclidean.objective_fun <- function(target, w, h, eta, beta){
  
  1/2 * ( sum( (target - (w %*% h))^2 ) 
          # + eta * sum(w^2) 
          # + beta * sum( colSums( h )^2 )
  )
}

KL.objective_fun <- function(y, w, h){
  
  a <-  w %*% h
  return (sum(y * log(y/a) - a + y))
}

#### End


### Utils
## Extended sinkhorn procedure (to get true coordinates)
sinkhorn_transform <- function(data_V, data_W, data_H, n_iters = 20) {
  D_vs_row_inv <- list()
  D_vs_column_inv <- list()
  D_hs_row_inv <-  list()
  D_ws_col_inv <-  list()
  
  
  V_row <- NULL
  V_column <- data_V
  
  H_row  <- NULL
  H_column <-  data_H
  W_row  <- NULL
  W_column <-  data_W
  
  for (i in seq_len(n_iters)) {
    D_vs_row_inv[[i]] <- rowSums(V_column)
    D_hs_row_inv[[i]] <- rowSums(H_column)
    #print(paste("Dh", toString(dim( diag(D_hs_row_inv[[i]])))))
    
    V_row <- diag(1 / D_vs_row_inv[[i]]) %*% V_column
    H_row <- diag(1 / D_hs_row_inv[[i]]) %*% H_column
    W_row <- diag(1 / D_vs_row_inv[[i]]) %*% W_column %*% diag(D_hs_row_inv[[i]])
    #print(paste("V_row", toString(dim(V_row))))
    #print(paste("W_row", toString(dim(W_row))))
    #print(paste("H_row", toString(dim(H_row))))
    D_vs_column_inv[[i]] <- colSums(V_row)
    D_ws_col_inv[[i]] <- colSums(W_row)
    #print(paste("Dw", toString(dim(diag(D_ws_col_inv[[i]])))))
    #print(paste("Dv_col", toString(dim(diag(D_vs_column_inv[[i]])))))
    
    V_column <- V_row %*% diag(1 / D_vs_column_inv[[i]])
    W_column <- W_row %*% diag(1 / D_ws_col_inv[[i]])
    
    H_column <- diag(D_ws_col_inv[[i]]) %*% H_row %*% diag(1 / D_vs_column_inv[[i]])
    #print(paste("V_col", toString(dim(V_column))))
    #print(paste("W_col", toString(dim(W_column))))
    #print(paste("H_col", toString(dim(H_column))))
  }
  list(
    V_row = V_row,
    V_column = V_column,
    W_row = W_row,
    W_column = W_column,
    H_row = H_row,
    H_column = H_column,
    D_vs_row_inv = D_vs_row_inv,
    D_vs_column_inv = D_vs_column_inv,
    D_hs_row_inv = D_hs_row_inv,
    D_ws_column_inv = D_ws_col_inv
  )
}



