library('NMF')
library(ggplot2)
library(gtable)
library(data.table)

library(ggsignif)
library('ComplexHeatmap')
library('Rcpp')
library('GEOquery')
library(linseed)
library(dplyr)
library(raster)
library(grid)
library(gridExtra)
library(ggrastr)
library(ggpubr)

require(linseed2)

#### Data Generation methods ####
make_syntetic_data <- function(M, N, K ) {
  V <- syntheticNMF(M, K, N, noise=TRUE)
  rownames(V) <- c(paste0("feature", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  V$coefficients
  true_W <- V$basis
  
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
make_picture_data <- function(M, N, K, pictures_W ) {
  
  H <- linseed:::sampleFromSimplexUniformly(N, K, 10000)
  V <- pictures_W %*% H
  
  rownames(V) <- c(paste0("pixels", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  H
  true_W <- pictures_W
  
  return (list(V=V, true_H=true_H, true_W=true_W))
  
}

## Generate all matrices using provides pictures_W matrix. H is random.
make_picture_data_random_mixtures <- function(M, N, K, pictures_W ) {
  
  H <- matrix(runif(K*N),nrow=K)
  V <- pictures_W %*% H
  
  rownames(V) <- c(paste0("pixels", 1:nrow(V)))
  colnames(V) <- c(paste0("mixture", 1:ncol(V)))
  true_H <-  H
  true_W <- pictures_W
  
  return (list(V=V, true_H=true_H, true_W=true_W))
  
}


#### End section


#### Logic methods ####

# Main realization of LinseedNMF method. Should be passed to nmf() package. Will train linseedv2 object
linseedv2_nmf_algorithm <- function(x, seed){
  print("New run")
  coef_hinge_H <-  1
  coef_hinge_W <-  1
  coef_pos_D_h <-  0
  coef_pos_D_w <-  0
  
  
  factorization_rank <- nbasis(seed)
  lo2 <- Linseed2Solver$new()
  lo2$set_data(x)
  lo2$project(factorization_rank)
  lo2$init_solution("random")
  
  iterations <-  300
  
  #print(paste("number of components", factorization_rank))
  
  
  
  NUM_RUNS <-  20
  LR_DECAY_STEPS <-  3
  lr_decay <- 0.5
  
  lr_x <-  0.1
  lr_omega <- 0.1
  
  for  (lr_step in 1:LR_DECAY_STEPS) {
    pb <- progress_bar$new(total = NUM_RUNS)
    for (i in 1:NUM_RUNS) {
      pb$tick()
      #print(paste("\titeration",x ))
      # starts C++ backend for optimization with gradient descent
      # optimization will run in 2 stages
      # number of steps for each stage specified with global_iterations parameter
      lo2$optim_solution(iterations, 
                         optim_config(coef_hinge_H = coef_hinge_H,
                                      coef_hinge_W = coef_hinge_W,
                                      coef_der_X = lr_x, 
                                      coef_der_Omega = lr_omega
                         ))
      
    }
    lr_x <- lr_x * lr_decay
    lr_omega <-  lr_omega * lr_decay
    
  }
  solution <- lo2$finalize_solution()
  
  # W matrix
  # W matrix
  W <- solution$W
  print(dim(W))
  W[W < 0] <-  0
  basis(seed) <- W
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  
  H <- solution$H
  H[H < 0] <-  0
  coef(seed) <- H
  # use the rotated matrix to get the mixture coefficient
  # use a scaling factor (just to illustrate the use of extra parameters)
  #coef(seed) <- t(abs(pca$x[,1:factorization.rank])) / scale.factor
  # return updated data
  seed$obj <- lo2
  
  
  distance <- nmfDistance(objective(seed))(seed, x)
  print(paste("objective is ", distance) )
  if (is.na(distance)) {
    print("Seems like objective is NA. Could happen randomly. Reruning the training")
    seed <- linseedv2_nmf_algorithm(x, seed)
  }
  
  return(seed)
}


#### End section







#### Plot methods ####

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
  
  sqrt((1/(length(w)) * sum( (target - w )^2 )) 
       # + eta * sum(w^2) 
       # + beta * sum( colSums( h )^2 )
  )
}

# Will calculate RMSE distances for each component (column)
get_rmse_distances_by_column <- function(original_matrix, curr_matrix){
  distances <-  sapply(c(1:current_K), function(original_component_index){
    candidates <- lapply(c(1:current_K), function(current_component_index){
      # print(paste(original_component_index,current_component_index))
      return (normalized_RMSE_function(original_matrix[,original_component_index], curr_matrix[,current_component_index]))
    })
    # print(unlist(candidates))
    return (min(unlist(candidates)))
  })
  return (data.table(component=c(1:current_K), rmse=unlist(distances)))
}

# Will calculate RMSE distances for each component (row)
get_rmse_distances_by_row <- function(original_matrix, curr_matrix){
  current_K <-  dim(original_matrix)[[1]]
  distances <-  sapply(c(1:current_K), function(original_component_index){
    
    candidates <- lapply(c(1:current_K), function(current_component_index){
      # print(paste(original_component_index,current_component_index))
      
      return (normalized_RMSE_function(original_matrix[original_component_index,], curr_matrix[current_component_index,]))
    })
    # print(unlist(candidates))
    return (min(unlist(candidates)))
  })
  return (data.table(component=c(1:current_K), rmse=unlist(distances)))
}

## Will calculate RMSE distances for each method and each run. W matrix.
get_W_rmse_values <- function(res_methods, data_used) {
  
  
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    
    
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <-  data_used[[current_run]]
      
      curr_W <- run_result@fit@W
      true_W <-curr_data$true_W
      
      distances_table <- get_rmse_distances_by_column(curr_data$true_W, curr_W)
      
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
get_H_rmse_values <- function(res_methods, data_used) {
  
  
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <- data_used[[current_run]]
      curr_H <- run_result@fit@H
      true_H <- curr_data$true_H
      
      distances_table <- get_rmse_distances_by_row(curr_data$true_H, curr_H)
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


get_V_rmse_values <- function(res_methods, data_used) {
  
  
  per_method_tables <- lapply(names(res_methods), function (method_name){
    curr_method_result <- res_methods[[method_name]]
    per_run_tables <- lapply(c(1:length(curr_method_result)), function(current_run) {
      run_result <- curr_method_result[[current_run]]
      curr_data <-  data_used[[current_run]]
      curr_H <- run_result@fit@H
      curr_W <- run_result@fit@W
      curr_H[curr_H < 0] <-  0
      curr_W[curr_W < 0] <-  0
      true_H <- curr_data$true_H
      true_W <- curr_data$true_W
      
      distance <- normalized_RMSE_function(true_W %*% true_H, curr_W %*% curr_H)
      distances_table <-  data.table(rmse=distance)
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




get_significance_values <-  function(distance_matrix, comparisons, vertical_gap = 0.03) {
  stat_test_results <-  lapply(c(1:length(comparisons)), function(i){
    comparison <- comparisons[[i]]
    data <- distance_matrix[distance_matrix$method %in% comparison]
    p.values <- sapply(split(data, data$component), 
                       function(x){
                         wilcox.test(rmse~method, data=x)$p.value})
    labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","n.s."))
    
    y.values <- sapply(split(distance_matrix, distance_matrix$component), function(x){max(x$rmse)})+ i*vertical_gap
    #y.values <- sapply(split(distance_matrix, distance_matrix$component), function(x){max(x$rmse)})+ i*0.2 # log
    return (list(labels=labels, y.values=y.values))
  })
  return (stat_test_results)
}


## Plot main boxplot with RMSE statistics for multiple runs, multiple methods for each component. (All 3 matrices)
plot_metrics_box <- function(res_methods, data_used, folder_to_save_results, title="Matrices") {
  W_matrix_distances <-  get_W_rmse_values(res_methods = res_methods, data_used = data_used)
  W_matrix_distances[method == 'Linseedv2 NMF']$method <- 'DualSimplex'
  H_matrix_distances <-  get_H_rmse_values(res_methods = res_methods, data_used = data_used)
  H_matrix_distances[method == 'Linseedv2 NMF']$method <- 'DualSimplex'
  V_matrix_distances <-  get_V_rmse_values(res_methods = res_methods, data_used = data_used)
  V_matrix_distances[method == 'Linseedv2 NMF']$method <- 'DualSimplex'
  W.summary <- W_matrix_distances %>% 
    group_by(method, component) %>% 
    summarise(
      sd = sd(rmse, na.rm=TRUE),
      value = mean(rmse, na.rm=TRUE)
    )
  
  H.summary <- H_matrix_distances %>% 
    group_by(method, component) %>% 
    summarise(
      sd = sd(rmse,na.rm=TRUE),
      value = mean(rmse, na.rm=TRUE)
    )
  
  V.summary <- V_matrix_distances %>% 
    group_by(method) %>% 
    summarise(
      sd = sd(rmse, na.rm=TRUE),
      value = mean(rmse, na.rm=TRUE)
    )
  
  
  comparisons = list(c("nsNMF", "DualSimplex"), c("lee", "DualSimplex"),c("brunet", "DualSimplex"))
  W_matrix_distances$method <- relevel(factor(W_matrix_distances$method), ref="DualSimplex")
  H_matrix_distances$method <- relevel(factor(H_matrix_distances$method), ref="DualSimplex")
  V_matrix_distances$method <- relevel(factor(V_matrix_distances$method), ref="DualSimplex")
  
  median_df_W <- W_matrix_distances %>%  
    group_by(method, component) %>%  
    summarize(average = median(rmse)) %>% 
    ungroup() 
  median_df_H <- H_matrix_distances %>%  
    group_by(method, component) %>%  
    summarize(average = median(rmse)) %>% 
    ungroup() 
  
  
  # Plot for W matrix
  gw <- ggplot(W_matrix_distances, aes(x=factor(component), y=rmse, fill=method)) +
    
    geom_boxplot(position = position_dodge2(width = 0.1)) +  
    geom_line(data=median_df_W, aes(x=component, y=average, color=method), 
              position = position_jitterdodge(jitter.width = 0.1),  size=1.2, alpha=0.4) +  
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size=0.4) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25) +
    theme(plot.title = element_text(face='bold', size=21)) +
    theme(axis.title.x  = element_text(size=20)) +
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=14)) +
    theme(axis.text.y  = element_text(size=14)) +
    theme(legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(plot.margin = margin(t = 2, r = 10, b = 2, l = 10, unit = "pt")) +
    
    labs(title = "RMSE error to original matrix W",
         x="Component", 
         y="RMSE" #subtitle = "LOL", x="LOLX", y="LOLY"
    )
  
  stat_test_results <- get_significance_values(W_matrix_distances, comparisons)
  gw <- gw + geom_signif(y_position = stat_test_results[[1]]$y.values, 
                         xmin = unique(W_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(W_matrix_distances$component)-(0.5/4)+ 0.04 ,
                         annotations = stat_test_results[[1]]$labels,
                         tip_length = 0)
  
  gw <- gw + geom_signif(y_position = stat_test_results[[2]]$y.values, 
                         xmin = unique(W_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(W_matrix_distances$component)+(0.5/4) - 0.02 ,
                         annotations = stat_test_results[[2]]$labels,
                         tip_length = 0)
  
  gw <- gw + geom_signif(y_position = stat_test_results[[3]]$y.values, 
                         xmin = unique(W_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(W_matrix_distances$component)+ (2*0.5/4)+ 0.04 ,
                         annotations = stat_test_results[[3]]$labels,
                         tip_length = 0)
  
  per_component_w_plot_list <- lapply(unique(W_matrix_distances$component), function(component_num) {
    curr_W_distances <-  W_matrix_distances[component == component_num]
    curr_plot <- ggplot(curr_W_distances, aes(x=method, y=rmse, fill=method)) +
      geom_boxplot(position = position_dodge2(width = 0.1)) +
      geom_signif(comparisons=comparisons, 
                  map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
                  margin_top = 0.1, tip_length = 0, textsize = 6,# vjust=0.1,
                  step_increase = 0.15)+
      geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
      scale_fill_brewer(palette = "Set1")+
      scale_color_brewer(palette = "Set1") +
      theme_classic(base_size=25) +
      labs(title = paste0("RMSE to W column ", component_num),x="Method", y="RMSE" #subtitle = "LOL", x="LOLX", y="LOLY"
      ) + 
      theme(plot.title = element_text(face='bold', size=21)) +  
      theme(axis.title.x  = element_text(size=20)) +  
      theme(axis.title.y  = element_text(size=20))+
      theme(axis.text.x  = element_text(size=14)) +  
      theme(axis.text.y  = element_text(size=14)) + 
      theme(legend.text=element_text(size=16))+
      theme(legend.title=element_text(size=16)) +
      theme(plot.margin = margin(t = 3, r = 10, b = 3, l = 10, unit = "pt")) +
      ylim(0, 1.5 * max(curr_W_distances$rmse))
    return (curr_plot)
    
  })
  
  # H matrix
  gh <- ggplot(H_matrix_distances, aes(x=factor(component), y=rmse, fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1)) +  
    geom_line(data=median_df_H, aes(x=component, y=average, color=method), 
              position = position_jitterdodge(jitter.width = 0.1),  size=1.2, alpha=0.4) +  
    geom_point(position = position_jitterdodge(jitter.width = 0.1), size=0.4) +
    scale_fill_brewer(palette = "Set1")+
    scale_color_brewer(palette = "Set1") +
    theme_classic(base_size=25) +
    theme(plot.title = element_text(face='bold', size=21)) +
    theme(axis.title.x  = element_text(size=20)) +
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=14)) +
    theme(axis.text.y  = element_text(size=14)) +
    theme(legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16)) +
    theme(plot.margin = margin(t = 2, r = 10, b = 2, l = 10, unit = "pt")) +
    labs(title = "RMSE to original matrix H",x="Component", y="RMSE"  #subtitle = "LOL", x="LOLX", y="LOLY"
    )  
  stat_test_results <- get_significance_values(H_matrix_distances, comparisons, vertical_gap = 0.015)
  gh <- gh + geom_signif(y_position = stat_test_results[[1]]$y.values, 
                         xmin = unique(H_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(H_matrix_distances$component)-(0.5/4)+ 0.04 ,
                         annotations = stat_test_results[[1]]$labels,
                         tip_length = 0)
  
  gh <- gh + geom_signif(y_position = stat_test_results[[2]]$y.values, 
                         xmin = unique(H_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(H_matrix_distances$component)+(0.5/4) - 0.02 ,
                         annotations = stat_test_results[[2]]$labels,
                         tip_length = 0)
  
  gh <- gh + geom_signif(y_position = stat_test_results[[3]]$y.values, 
                         xmin = unique(H_matrix_distances$component)- (2*0.5/4),
                         xmax = unique(H_matrix_distances$component)+ (2*0.5/4)+ 0.04 ,
                         annotations = stat_test_results[[3]]$labels,
                         tip_length = 0)
  
  per_component_h_plot_list <- lapply(unique(H_matrix_distances$component), function(component_num) {
    curr_H_distances <-  H_matrix_distances[component == component_num]
    curr_plot <- ggplot(curr_H_distances, aes(x=method, y=rmse, fill=method)) +
      geom_boxplot(position = position_dodge2(width = 0.1)) +
      geom_signif(comparisons=comparisons, 
                  map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
                  margin_top = 0.1, tip_length = 0, textsize = 6,# vjust=0.1,
                  step_increase = 0.2)+
      geom_point(position =position_jitterdodge(jitter.width = 0.1)) +
      scale_fill_brewer(palette = "Set1")+
      scale_color_brewer(palette = "Set1") +
      theme_classic(base_size=25) +
      labs(title = paste0("RMSE to H row ", component_num),x="Method", y="RMSE" #subtitle = "LOL", x="LOLX", y="LOLY"
      ) + 
      theme(plot.title = element_text(face='bold', size=21)) +  
      theme(axis.title.x  = element_text(size=20)) +  
      theme(axis.title.y  = element_text(size=20))+
      theme(axis.text.x  = element_text(size=14)) +  
      theme(axis.text.y  = element_text(size=14)) + 
      theme(legend.text=element_text(size=16))+
      theme(legend.title=element_text(size=16)) +
      theme(plot.margin = margin(t = 3, r = 10, b = 3, l = 10, unit = "pt")) +
      ylim(0, 1.7 * max(curr_H_distances$rmse))
    return (curr_plot)
    
  })
  
  
  
  gv <- ggplot(V_matrix_distances, aes(x=method, y=rmse, fill=method)) +
    geom_boxplot(position = position_dodge2(width = 0.1))+
    geom_signif(comparisons=comparisons, 
                map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "n.s."=2),
                margin_top = 0.1, tip_length = 0,
                step_increase = 0.1,)+
    scale_fill_brewer(palette = "Set1")+
    geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
    scale_y_continuous(trans='log10') +
    labs(title = "RMSE to original matrix V",x="Method", y="RMSE"  #subtitle = "LOL", x="LOLX", y="LOLY"
    )+ 
    theme_classic(base_size=25) +
    theme(plot.title = element_text(face='bold', size=21)) +  
    theme(axis.title.x  = element_text(size=20)) +  
    theme(axis.title.y  = element_text(size=20))+
    theme(axis.text.x  = element_text(size=14)) +  
    theme(axis.text.y  = element_text(size=14)) + 
    theme(legend.text=element_text(size=16))+
    theme(legend.title=element_text(size=16))+
    theme(plot.margin = margin(t = 2, r = 10, b = 2, l = 10, unit = "pt"))
  
  
  t <- ggarrange(gw, gh, gv, nrow = 1, common.legend = TRUE, legend='right')
  #  top = title
  #  )%>%
  #  ggexport(filename = filename, width = 18, height = 5)
  filename <- paste0(folder_to_save_results,  title, ".svg")
  ggsave(file=filename, plot=t, width=18, height=5, device=svg)
  filename <- paste0(folder_to_save_results, "V_", title, ".svg")
  ggsave(filename, gv, width=7.5, height=5, device=svg)
  filename <- paste0(folder_to_save_results, "H_", title, ".svg")
  ggsave(filename, gh, width=7.5, height=5, device=svg)
  filename <- paste0(folder_to_save_results, "W_", title, ".svg")
  ggsave(filename, gw, width=7.5, height=5, device=svg)
  for (comp_ind in 1:length(per_component_w_plot_list)) {
    filename <- paste0(folder_to_save_results, "W_", comp_ind, title, ".svg")
    ggsave(filename, per_component_w_plot_list[[comp_ind]], width=7.5, height=5, device=svg)
    filename <- paste0(folder_to_save_results,"H_", comp_ind, title, ".svg")
    ggsave(filename, per_component_h_plot_list[[comp_ind]], width=7.5, height=5, device=svg)
  }
  per_comp_common_w <- ggarrange(plotlist=per_component_w_plot_list, nrow = 1, common.legend = TRUE, legend='right')
  filename <- paste0(folder_to_save_results,"W_common_", title, ".svg")
  ggsave(filename, per_comp_common_w, width=23, height=5, device=svg)
  
  per_comp_common_h <- ggarrange(plotlist=per_component_h_plot_list, nrow = 1, common.legend = TRUE, legend='right')
  filename <- paste0(folder_to_save_results,"H_common_", title, ".svg")
  ggsave(filename, per_comp_common_h, width=23, height=5, device=svg)
  
  
  return (list(W_matrix_distances=W_matrix_distances, 
               H_matrix_distances=H_matrix_distances,
               V_matrix_distances=V_matrix_distances, 
               plot=t,
               V_plot=gv,
               H_plot=gh,
               W_plot=gw,
               per_comp_plot_w=per_comp_common_w,
               per_comp_plot_h=per_comp_common_h
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
  linseed_results <- res_methods[["Linseedv2 NMF"]]
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
    # plot(raster(matrix(W_matrix[, picture_index], nrow = size, ncol = size, byrow = TRUE)),
    #      col = gray.colors(100, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), 
    #      axes=FALSE, box=FALSE, legend=FALSE, ann=FALSE)
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


plot_pictures_for_result <- function(method_results, model_index, folder_to_save_results) {
  res <- lapply(method_results, function (method_multiple_results){
    method_result <- method_multiple_results[[model_index]]
    curr_W <- method_result@fit@W
    plot_pictures(curr_W, method_result@method, folder_to_save_results)
    #par(oma=c(0,0,0,0)) # all sides have 3 lines of space
  })
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


