# Non-negative matrix factorization and deconvolution as dual simplex problem

This repository is an official starting point to explore [Dual Simplex NMF/deconvolution method](https://www.biorxiv.org/content/10.1101/2024.04.09.588652v1)
It contains code to reproduce figures from the paper, and at the same time, provides examples on how to use the [DualSimplex package](https://github.com/artyomovlab/dualsimplex).

```
Non-negative matrix factorization and deconvolution as dual simplex problem
Denis Kleverov, Ekaterina Aladyeva, Alexey Serdyukov, Maxim Artyomov
bioRxiv 2024.04.09.588652; doi: https://doi.org/10.1101/2024.04.09.588652
```

### Project structure
```
- data — all the external data, used in figures
- figures — notebooks for figures reproduction
- out — generated svgs and dualsimplex checkpoints will be placed here
- R — supporting code, imported in figures
```

### Running
1. Select a figure to reproduce.
2. Script `setup.R` (executed at the beginnig of the each script) will install the DualSimplex package using the github
3. Some figures require additional data, copy [large.tar.gz](https://drive.google.com/file/d/1BP2jk5ug4-UNKkLew84miL8JtrdY17Ng/view?usp=sharing) into `data/large`. 
4. Go to the `figures` directory and open the corresponding notebook.
5. Run cells in the notebook one by one. Optionally, tweak some parameters to see alternative outcomes.
6. See resulting figures in the `out` directory.


## Figures in this repository
### 2. Sinkhorn procedure
Simple visualization of the Sinkhorn procedure applied to factorizable matrix 
([2_sinkhorn_visualization.Rmd](figures/2_sinkhorn_visualization.Rmd))

### 3. Main algorithm
Deconvolution of simulated bulk RNA-seq gene expression dataset with main approach 
([3c_simulated_gene_expression_main_algorithm.Rmd](figures/3c_simulated_gene_expression_main_algorithm.Rmd))
### 4. Picture unmixing with NMF
- Simple synthetic NMF case just for testing
([4_1_nmf_setup_ensure_working_with_simulated_data.Rmd](figures/4_1_nmf_setup_ensure_working_with_simulated_data.Rmd))
- Picture unmixing 
    - Generate the data, solve problem with different methods ([4_2_nmf_noisy_pictures.Rmd](figures/4_2_nmf_noisy_pictures.Rmd))
    - Visualize results  ([4_3_make_result_nmf_figure_plots.Rmd](figures/4_3_make_result_nmf_figure_plots.Rmd))

### 5. Complete deconvolution of bulk RNA-seq data
-  Brain/Liver/Lung mixtures GSE19830 ([5a_gse19830.Rmd](figures/5a_gse19830.Rmd))
-  4 immune cell types GSE11058 ([5b_gse11058.Rmd](figures/5b_gse11058.Rmd))
-  TCGA HNSC dataset ([5c_bulk_hnsc_preprocessing.Rmd](figures/5c_bulk_hnsc_preprocessing.Rmd)), ([5de_bulk_hnsc_deconvolution.Rmd](figures/5de_bulk_hnsc_deconvolution.Rmd)), ([5f_bulk_hnsc_single_cell_validation.Rmd](figures/5f_bulk_hnsc_single_cell_validation.Rmd)), ([5g_bulk_hnsc_clinical_correlations.Rmd](figures/5g_bulk_hnsc_clinical_correlations.Rmd))

### 6. DREAM challenge data analysis with Dual Simplex
- Get marker genes from pure samples
([6_1_dream_challenge_prepare_markers.Rmd](figures/6_1_dream_challenge_prepare_markers.Rmd))
- Reference free solution
([6_2_dream_challenge_reference_free.Rmd](figures/6_2_dream_challenge_reference_free.Rmd))
- Signature based solution
([6_3_dream_challenge_signature_based.Rmd](figures/6_3_dream_challenge_signature_based.Rmd))

### 7. Special cases when our method is good
- Deconvolution with missing component-specific features (corners) 
([7_1_deconvolition_with_missing_corners.Rmd](figures/7_1_deconvolition_with_missing_corners.Rmd))
- Positivity optimization. just alternative approach to optimization in projected space
([7_2_simulated_gene_expression_alternative_approach.Rmd](figures/7_2_simulated_gene_expression_alternative_approach.Rmd))
### S3. NMF with simulated data matrices
-  Factorize simulated dataset with multiple methods ([s3c_1_nmf_simulated_data.Rmd](figures/s3c_1_nmf_simulated_data.Rmd))
-  Visualize results ([s3c_1_nmf_simulated_data.Rmd](figures/s3c_2_make_result_figure_plots.Rmd))

### S4. Further analysis for TCGA HNSC bulk RNA-seq dataset
Pathway analysis, signature genes expression heatmap, multple initializations  ([s4_hnsc_further.Rmd](figures/s4_hnsc_further.Rmd))

### S5. Signature base deconvolution with DualSimpelx approach
-  Generate data ([s5_1_signature_based_deconvolution_make_mixed_data.Rmd](figures/s5_1_signature_based_deconvolution_make_mixed_data.Rmd))
-  Brain/Liver/Lung mixtures GSE19830  ([s5_2_signature_based_deconvolution_GSE11058.Rmd](figures/s5_2_signature_based_deconvolution_GSE11058.Rmd))
-  Simulated brain dataset ([s5_3_signature_based_deconvolution_brain.Rmd](figures/s5_3_signature_based_deconvolution_brain.Rmd))
-  
### S7. Single cell analysis (clustering)
-  Single cell analysis script([s7_single_cell_hnsc.Rmd](figures/s7_single_cell_hnsc.Rmd))
### S8. Multiple solution NMF. How our method behaves
-  Solution of multiple solution NMF ([s8_multiple_solution_deconvolution.Rmd](figures/s8_multiple_solution_deconvolution.Rmd))
### Supplementary notes scripts
-  Script used to generate supplementary note about optimization ([supplementary_notes_scripts.Rmd](figures/supplementary_notes_scripts.Rmd))


## Authors
Contributors names and contact info 
- Denis Kleverov ([@denis_kleverov](https://twitter.com/denis_kleverov)) ([linkedIn](https://linkedin.com/in/denklewer)) ([email](mailto:denklewer@gmail.com))
-  Ekaterina Aladyeva ([AladyevaE](https://twitter.com/AladyevaE)) 
-  Alexey Serdyukov  ([email](mailto:leshaserdyukov@gmail.com))
-  prof. Maxim Artyomov ([@maxim_artyomov](https://twitter.com/maxim_artyomov)) ([email](mailto:martyomov@wustl.edu))


## Troubleshooting
> Dependency: package 'xxx' is not available (for R version x.y.z)

Install package directly from source link from CRAN. For example:

```install.packages(https://cran.r-project.org/src/contrib/RcppML_0.3.7.tar.gz, repos = NULL)```



> Can't plot UMAP with plot_projected on Mac

Unfortunately, umap library has a bug (only on MacOS) that doesn't allow to
add new points to umap after it's calculated, which is crucial for DualSimplex.
If that is the case for you, call `plot_projected(use_dims = 2:3)`,
or other dimensions, to see simplexes without dimensionality reduction.
