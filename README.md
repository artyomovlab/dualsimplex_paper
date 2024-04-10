# Non-negative matrix factorization and deconvolution as dual simplex problem
This repository is an official starting point to explore [Dual Simplex NMF/deconvolution method](https://arxiv.org/abs/2030.12345)


```
Paper full reference here
```

Here we posted all the scripts/data necessary to reproduce figures from our paper
## Requirements
You will need to install R

## Figures in this repository
### 2. Sinkhorn procedure
Simple visualization of the Sinkhorn procedure applied to factorizable matrix 
([2_sinkhorn_visualization.Rmd](figures/2_sinkhorn_visualization.Rmd))

### 3. Main algorithm
Deconvolution of simulated bulk RNA-seq gene expression dataset with main approach 
([3_simulated_gene_expression_main_algorithm.Rmd](figures/3_simulated_gene_expression_main_algorithm.Rmd))
### 4. Minimal formulation

Deconvolution of simulated bulk RNA-seq gene expression dataset with alternative aproach 
([4_simulated_gene_expression_alternative_approach.Rmd](figures/4_simulated_gene_expression_alternative_approach.Rmd))

`you will need to switch the branch`
### 5. Picture unmixing with NMF
- Simple synthetic NMF case just for testing
([5_1_nmf_setup_ensure_working_with_simulated_data.Rmd](figures/5_1_nmf_setup_ensure_working_with_simulated_data.Rmd))
- Picture unmixing 
    - Generate the data, solve problem with different methods ([5_2_nmf_noizy_pictures.Rmd](figures/5_2_nmf_noizy_pictures.Rmd))
    - Visualize results  ([5_3_make_result_nmf_figure_plots.Rmd](figures/5_3_make_result_nmf_figure_plots.Rmd)
)
### 6. Single cell data
-  Method aplied to single cell clustering ([example](google.com))

### 7. Complete deconvolution of bulk RNA-seq data
-  Brain/Liver/Lung mixtures GSE19830 ([7a_gse19830.Rmd](figures/7a_gse19830.Rmd))
-  4 immune cell types GSE11058 ([7b_gse11058_7b](figures/7b_gse11058_7b.Rmd))
-  TCGA HNSC dataset ([example](google.com))

### S3. NMF with simulated data matrices
-  Factorize simulated dataset with multiple methods ([s3c_1_nmf_simulated_data.Rmd](figures/s3c_1_nmf_simulated_data.Rmd))
-  Visualize results  ([s3c_2_make_result_figure_plots.Rmd](figures/s3c_2_make_result_figure_plots.Rmd))


### S4. Different number of clusters for single cell data
Comparison of the clustering solutions for different methods ([example](google.com))


### S5. Further analysis for TCGA HNSC bulk RNA-seq dataset
Pathway analysis, signature genes expression heatmap, multple initializations ([example](google.com))

### S6. Signature base deconvolution with DualSimpelx approach
-  Brain/Liver/Lung mixtures GSE19830 ([example](google.com))
-  Simulated brain dataset ([example](google.com))

## Authors
Contributors names and contact info 
- Denis Kleverov ([@denis_kleverov](https://twitter.com/denis_kleverov)) ([linkedIn](https://linkedin.com/in/denklewer) )
-  Ekaterina Aladyeva ([AladyevaE](https://twitter.com/AladyevaE)) 
-  Alexey Serdyukov  ([linkedIn](https://www.linkedin.com/in/alexey-serdyukov-52624b213/))
-  prof. Maxim Artyomov ([@maxim_artyomov](https://twitter.com/maxim_artyomov)) ([email](martyomov@wustl.edu))



