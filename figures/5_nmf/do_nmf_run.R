## This is same script as 2_nmf_noizy....   just to run this outside R studio

rm(list = ls())

# path to save everything
dir_to_save_fig <-  "../../out/NMF_simulated/"
dir.create(file.path(".", dir_to_save_fig), showWarnings = F, recursive = T)
unloadNamespace("DualSimplex")
devtools::load_all(path='../../../linseed2/')
source('nmf_methods.R')
library(linseed)
library(raster)

run_experiment_picture_random_noizy <- function(current_K, current_M, current_N, pictures_W, n_data_generations = 3, 
                                                start_with_hinge_W=1, start_with_hinge_H=1,
                                                LR_DECAY_STEPS = 7, RUNS_EACH_STEP = 500,    PARAMETERS_INCREASE_STEPS = 2, 
                                                lr_x =  0.1,
                                                lr_omega = 0.1,
                                                noise_deviation=3.5, noise_nature="proportional" ) {
  
  #data_generated <- lapply(c(1:n_data_generations), function(run_index) make_simulated_data(current_M, current_N, current_K, noise_deviation))
  
    data_generated <- lapply(c(1:n_data_generations), 
			     function(run_index) make_picture_data_random_mixtures(current_M, current_N, current_K, 
										   pictures_W,noise_deviation)) 
    total_result <- lapply(c(1:n_data_generations), function(run_index) {
    print(paste("data gen ", run_index))
    syntetic_data <- data_generated[[run_index]]
    current_res_methods <-  list()
    current_res_methods <- nmf(syntetic_data$noizy_V, current_K, list('snmf/l','ls-nmf','lee', 'brunet', 'nsNMF'), .options='vp10', nrun=10)
    
    print("Compute NMF method 'HALSacc'")
    halsacc_res <- nmf(syntetic_data$noizy_V, current_K, hNMF::HALSacc, nrun=10, .options='vp10', name="HALSacc")
    #halsacc_res <- nmf(syntetic_data$V, current_K, halsacc_nmf_algorithm, nrun=10, .options='tp1', name="HALSacc")
    current_res_methods[["HALSacc"]] <- halsacc_res
    #print("Compute NMF method 'PGNMF'")
    pgnmf_res <- nmf(syntetic_data$noizy_V, current_K, hNMF::PGNMF, nrun=10, .options='vp10', name="PGNMF")
    # # pgnmf_res <- nmf(syntetic_data$V, current_K, pgnmf_nmf_algorithm, nrun=10, .options='tp1', name="PGNMF")
    current_res_methods[["PGNMF"]] <- pgnmf_res
    #print("Compute NMF method 'ALS'")
    als_nmf_res <- nmf(syntetic_data$noizy_V, current_K, als_nmf_algorithm, nrun=10, .options='p10', name="ALS")
    current_res_methods[["ALS"]] <- als_nmf_res
    print("Compute NMF method 'DualSimplex'")
    linseed_res <- nmf(syntetic_data$noizy_V, 
                       current_K, 
                       dualsimplex_nmf_algorithm, 
                       nrun=10, .options='vp10', name="DualSimplex", 
                       start_with_hinge_W=start_with_hinge_W, 
                       start_with_hinge_H=start_with_hinge_H,
                       LR_DECAY_STEPS = LR_DECAY_STEPS, 
                       RUNS_EACH_STEP = RUNS_EACH_STEP, 
                       PARAMETERS_INCREASE_STEPS=PARAMETERS_INCREASE_STEPS, 
                       lr_x=lr_x, lr_omega=lr_omega)
    current_res_methods[["DualSimplex"]] <- linseed_res
    return(current_res_methods)
  })
  
  
  all_results_flat <- unlist(total_result, recursive = F)
  unique_names <- unique(names(all_results_flat))
  per_method_per_run_result <- lapply(unique_names, function(cur_nm) all_results_flat[names(all_results_flat) == cur_nm])
  names(per_method_per_run_result) <- unique_names
  # distances_to_check <- plot_metrics(res_methods = per_method_per_run_result, data_used=data_generated, title=paste("K =", as.character(current_K), "M =", as.character(current_M),"N =",as.character(current_N)))
  
  
return (list(data_used=data_generated, fit_objects=per_method_per_run_result))
}




components_count <-  c(3, 4, 5)
mixtures_count <-  c(25, 50, 100, 500)
features_count <- c(200, 500, 1000, 2000)

size <- 128
number_of_pictures_to_select <- 4

pictures_W <- rhdf5::h5read("img/merged_imgs.h5", "matrix")
pictures_W <- t(pictures_W)[,1:number_of_pictures_to_select]

#current_M <- 1000
#current_K <- 4
#current_N  <- 800
# For N = 500 use   start_with_hinge_W=0.1 start_with_hinge_H=0.1
# For N = 100 use   start_with_hinge_W=10 start_with_hinge_H=1
noise_nature <-  "proportional"
current_M <- nrow(pictures_W)
current_K <- ncol(pictures_W)
current_N <- 100
for (noise_deviation in c(0,1,2)) {
  print(paste("Current noise deviation", noise_deviation))
  result_random <- run_experiment_picture_random_noizy(current_K = current_K, current_M= current_M, current_N = current_N, 
                                                       n_data_generations = 10, 
                                                       pictures_W = pictures_W ,  
                                                       noise_deviation=noise_deviation,
                                                       noise_nature = noise_nature,
                                                       start_with_hinge_W=10, 
                                                       start_with_hinge_H=1,
                                                       LR_DECAY_STEPS = 11,
                                                       RUNS_EACH_STEP = 2000, 
                                                       PARAMETERS_INCREASE_STEPS = 3,
                                                       lr_x =  0.1,
                                                       lr_omega = 0.1)
  # distances_to_check <- plot_metrics_box(res_methods = result_random$fit_objects, 
  #                                      folder_to_save_results = dir_to_save_fig,
  #                                      data_used=result_random$data_used, 
  #                                       title=paste("random_noise_", as.character(noise_deviation), 
  #                                                   "_K =", as.character(current_K), "M =", as.character(current_M),"N =",as.character(current_N)))
  # print("Make plots")
  # result_random$distances_to_check <- distances_to_check
  print("Save to file")
  file_name <- paste0(dir_to_save_fig, noise_nature, "_random_noise_",  as.character(noise_deviation),  '_random_proportions_', 
                      current_K ,'K_', current_N, 'N_', current_M, 'M_random-H_10runs_10init_random_W10_H1.rds')
  saveRDS(result_random, file=file_name)
}

