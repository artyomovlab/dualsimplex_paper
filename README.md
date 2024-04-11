# Non-negative matrix factorization and deconvolution as dual simplex problem

This repository is an official starting point to explore [Dual Simplex NMF/deconvolution method](https://arxiv.org/abs/2030.12345)
It contains code to reproduce figures from the paper, and at the same time, provides examples on how to use the [DualSimplex package](https://github.com/artyomovlab/dualsimplex).

```
Paper full reference here
```

### Project structure
```
- data — all the external data, used in figures
- figures — notebooks for figures reproduction
- out — generated svgs and dualsimplex checkpoints will be placed here
- R — supporting code, imported in figures
```

### Running
1. Clone the DualSimplex repository, and put it in the same directory as this repository, for setup.R to work.
2. Select a figure to reproduce.
3. If you chose Figure 6 or Figure 7, download and unpack the contents of [large.tar.gz](https://drive.google.com/drive/folders/1SgL1sCW4ItfY1wH5kqh__5FJPO5W0ltS) into `data/large`. 
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
### 4. Minimal formulation

Deconvolution of simulated bulk RNA-seq gene expression dataset with alternative aproach 
([4d_simulated_gene_expression_alternative_approach.Rmd](figures/4d_simulated_gene_expression_alternative_approach.Rmd))

`you will need to switch the branch`
### 5. Picture unmixing with NMF
- Simple synthetic NMF case just for testing
([5_1_nmf_setup_ensure_working_with_simulated_data.Rmd](figures/5_1_nmf_setup_ensure_working_with_simulated_data.Rmd))
- Picture unmixing 
    - Generate the data, solve problem with different methods ([5_2_nmf_noizy_pictures.Rmd](figures/5_2_nmf_noisy_pictures.Rmd))
    - Visualize results  ([5_3_make_result_nmf_figure_plots.Rmd](figures/5_3_make_result_nmf_figure_plots.Rmd)
)
### 6. Single cell data
-  Method aplied to single cell clustering ([6_single_cell_hnsc.Rmd](figures/6_single_cell_hnsc.Rmd))

### 7. Complete deconvolution of bulk RNA-seq data
-  Brain/Liver/Lung mixtures GSE19830 ([7a_gse19830.Rmd](figures/7a_gse19830.Rmd))
-  4 immune cell types GSE11058 ([7b_gse11058_7b.Rmd](figures/7b_gse11058_7b.Rmd))
-  TCGA HNSC dataset ([7c_bulk_hnsc_preprocessing.Rmd](figures/7c_bulk_hnsc_preprocessing.Rmd)), ([7de_bulk_hnsc_deconvolution.Rmd](figures/7de_bulk_hnsc_deconvolution.Rmd)), ([7f_bulk_hnsc_single_cell_validation.Rmd](figures/7f_bulk_hnsc_single_cell_validation.Rmd)), ([7g_bulk_hnsc_clinical_correlations.Rmd](figures/7g_bulk_hnsc_clinical_correlations.Rmd))

### S3. NMF with simulated data matrices
-  Factorize simulated dataset with multiple methods
-  Visualize results


### S4. Different number of clusters for single cell data
Comparison of the clustering solutions for different methods


### S5. Further analysis for TCGA HNSC bulk RNA-seq dataset
Pathway analysis, signature genes expression heatmap, multple initializations

### S6. Signature base deconvolution with DualSimpelx approach
-  Brain/Liver/Lung mixtures GSE19830
-  Simulated brain dataset


## Authors
Contributors names and contact info 
- Denis Kleverov ([@denis_kleverov](https://twitter.com/denis_kleverov)) ([linkedIn](https://linkedin.com/in/denklewer) )
-  Ekaterina Aladyeva ([AladyevaE](https://twitter.com/AladyevaE)) 
-  Alexey Serdyukov  ([email](mailto:leshaserdyukov@gmail.com))
-  prof. Maxim Artyomov ([@maxim_artyomov](https://twitter.com/maxim_artyomov)) ([email](mailto:martyomov@wustl.edu))



## Figures in this repository
### 2. Sinkhorn procedure
Simple visualization of the Sinkhorn procedure applied to factorizable matrix 
([2_sinkhorn_visualization.Rmd](figures/2_sinkhorn_visualization.Rmd))

### 3. Main algorithm
Deconvolution of simulated bulk RNA-seq gene expression dataset with main approach 
([3c_simulated_gene_expression_main_algorithm.Rmd](figures/3c_simulated_gene_expression_main_algorithm.Rmd))
### 4. Minimal formulation

Deconvolution of simulated bulk RNA-seq gene expression dataset with alternative aproach 
([4d_simulated_gene_expression_alternative_approach.Rmd](figures/4d_simulated_gene_expression_alternative_approach.Rmd))

`you will need to switch the branch`
### 5. Picture unmixing with NMF
- Simple synthetic NMF case just for testing
([5_1_nmf_setup_ensure_working_with_simulated_data.Rmd](figures/5_1_nmf_setup_ensure_working_with_simulated_data.Rmd))
- Picture unmixing 
    - Generate the data, solve problem with different methods ([5_2_nmf_noizy_pictures.Rmd](figures/5_2_nmf_noisy_pictures.Rmd))
    - Visualize results  ([5_3_make_result_nmf_figure_plots.Rmd](figures/5_3_make_result_nmf_figure_plots.Rmd)
)
### 6. Single cell data
-  Method aplied to single cell clustering ([6_single_cell_hnsc.Rmd](figures/6_single_cell_hnsc.Rmd))

### 7. Complete deconvolution of bulk RNA-seq data
-  Brain/Liver/Lung mixtures GSE19830 ([7a_gse19830.Rmd](figures/7a_gse19830.Rmd))
-  4 immune cell types GSE11058 ([7b_gse11058_7b.Rmd](figures/7b_gse11058_7b.Rmd))
-  TCGA HNSC dataset ([7c_bulk_hnsc_preprocessing.Rmd](figures/7c_bulk_hnsc_preprocessing.Rmd)), ([7de_bulk_hnsc_deconvolution.Rmd](figures/7de_bulk_hnsc_deconvolution.Rmd)), ([7f_bulk_hnsc_single_cell_validation.Rmd](figures/7f_bulk_hnsc_single_cell_validation.Rmd)), ([7g_bulk_hnsc_clinical_correlations.Rmd](figures/7g_bulk_hnsc_clinical_correlations.Rmd))

### S3. NMF with simulated data matrices
-  Factorize simulated dataset with multiple methods
-  Visualize results


### S4. Different number of clusters for single cell data
Comparison of the clustering solutions for different methods


### S5. Further analysis for TCGA HNSC bulk RNA-seq dataset
Pathway analysis, signature genes expression heatmap, multple initializations

### S6. Signature base deconvolution with DualSimpelx approach
-  Brain/Liver/Lung mixtures GSE19830
-  Simulated brain dataset


## Authors
Contributors names and contact info 
- Denis Kleverov ([@denis_kleverov](https://twitter.com/denis_kleverov)) ([linkedIn](https://linkedin.com/in/denklewer) )
-  Ekaterina Aladyeva ([AladyevaE](https://twitter.com/AladyevaE)) 
-  Alexey Serdyukov  ([email](mailto:leshaserdyukov@gmail.com))
-  prof. Maxim Artyomov ([@maxim_artyomov](https://twitter.com/maxim_artyomov)) ([email](mailto:martyomov@wustl.edu))


### Troubleshooting
> Dependency: package 'xxx' is not available (for R version x.y.z)

Install package directly from source link from CRAN. For example:

```install.packages(https://cran.r-project.org/src/contrib/RcppML_0.3.7.tar.gz, repos = NULL)```



> Can't plot UMAP with plot_projected on Mac

Unfortunately, umap library has a bug (only on MacOS) that doesn't allow to
add new points to umap after it's calculated, which is crucial for DualSimplex.
If that is the case for you, call `plot_projected(use_dims = 2:3)`,
or other dimensions, to see simplexes without dimensionality reduction.