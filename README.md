
# About CeL-Gen
The CeL-Gen method aims to map the cell-type specificity of eQTL from bulk gene expression data of heterogeneous tissues.
 
The method takes as input a large cohort of individuals, where the input for each individual includes:
- Genotyping. 
- Bulk expression of genes in a certain tissue. 
- The relative abundance (proportions) of the various cell types in the tissue (it is possible to use computational deconvolution methods to predict cell-type proportions from bulk genomics data (Newman et al., 2015)). 
- The known cell lineage tree that includes the cell types in the heterogeneous tissue under study. 

The output is the cell-type specificity of eQTLs: an inferred eQTL for each gene together with its effect size in each cell type, assuming that each eQTL has only two levels of genetic effect sizes, each of which exists in a subset of the cell types. 

In searching for the cell-type specificity, CeL-Gen limits the analysis based on the known architecture of the cell lineage. 
Particularly, rather than testing all cell-type subsets, CeL-Gen tests all differentiation steps (called ‘branches’) in the lineage tree – each branch partitions the cell types into two subsets (descendants/non-descendants) of distinct effect size. 
This strategy relies on the understanding that variation in effects typically coincides with the cell lineage and follows the fundamental principle of a parsimonious sequence of alterations. 


# Requirements:
* numpy~=1.16.3
* pandas~=0.24.2
* statsmodels~=0.11.1
* python 3.7

# Using CeL-Gen


        python3 eQTL_analysis_cel_gen.py <output_path> 
                                        [--genotype GENOTYPE_PATH] 
                                        [--gene_expression GENE_EXPRESSION_PATH]
                                        [--cell_fractions CELL_FRACTIONS_PATH]
                                        [--tree_structure TREE_PATH]
                                        [--gene_positions GENE_POSITIONS]
                                        [--snp_positions SNP_POSITIONS]

The files are expected to be in the following formats (example files are available in the example_data folder):

**GENE_EXPRESSION_PATH**

The location of the gene expression file. The script expects a csv file with the following format:
* The first column should contain the genes' names, as they appear in the GENE_POSITIONS files.
* Each of the next columns should describe the gene expression measured for a different sample.
* The header row should include the samples' names exactly as they appear in the genotype file.

Please note that any preprocessing should be done before using the script, as CelGen does not apply any normalization.

**GENOTYPE_PATH**

The location of the genotype file. The script expects a csv file with the following format:

* The first column should contain the SNPs' names, as they appear in the SNP_POSITIONS files.
* Each of the next columns should describe a different sample's genotype, converted to numbers (0-2).
* The header row should include the samples' names exactly as they appear in the gene expression file.

**TREE_PATH**

A file describing the known linage tree's structure. 

Each row on the file (except the header row) should describe a different cell type.
The file should contain the following columns:
* name - the name of the cell type, exactly as it appears in the cell fractions file.
* parent - the name of the cell type's parent cell type (or empty if the row represents the lineage tree root).
* children - the names of the children cell types, separated by ";".

**GENE_POSITIONS**

The path of a csv file describing the positions of each gene. The expected format:
* The first column should include the gene's name, exactly as it appears in the gene expression file.
* The next three columns should be:
    - chromosome_name: the name of the chromosome on which the gene is located. Should match the names in SNP_POSITIONS.
    - start_position: the start position of the gene (in bp).
    - end_position: the end position of the gene (in bp).

**SNP_POSITIONS**
The path of a csv file describing the positions of each SNP. The expected format:
* The first column should include the SNP's name, exactly as it appears in the genotype file.
* The next three columns should be:
    - chromosome: the name of the chromosome on which the gene is located. Should match the names in GENE_POSITIONS.
    - bp: the position of the SNP (in bp).
    
## Output file
### How to read the file
The output file includes the following columns:
#### General columns
* gene_id
* snp_id
#### Null model columns (describing a model with a one effect size for all the cell types)
* generic_log_likelihood - the log likelihood of the one-effect model.
* generic_ssr - sum of squared (whitened) residuals.
* generic_ess - the explained sum of squares.
* generic_centered_tss - the total (weighted) sum of squares centered about the mean.
* generic_r2 - R-squared of the one-effect model.
* generic_r2_adj - adjusted R-squared.
* generic_p_value - goodness of fit p-value.
* generic_aic - Akaike’s information criteria.
* generic_bic - Bayes’ information criteria.
* generic_genotype_coeff - the genotype term coefficient.
* generic_intercept_coeff.

#### Branch model columns
For cell type x, a model where x is the root of the affected tree (the branch between x and x's parent is the alteration switch) is fitted, and the following values are reported in the output file:
* <cell_type>_restricted_full_f_test - the p-value of an F-test comparing between the null model (a single effect model) 
  and the <cell_type>-subtree model. The p-value is -log10-transformed.  
* <cell_type>_full_log_likelihood - the log likelihood of the <cell_type>-subtree model. 
* <cell_type>_ssr - sum of squared (whitened) residuals.
* <cell_type>_ess - the explained sum of squares.
* <cell_type>_centered_tss - the total (weighted) sum of squares centered about the mean.
* <cell_type>_r2 - R-squared of the <cell_type>-subtree model.
* <cell_type>_r2_adj - adjusted R-squared.
* <cell_type>_p_value - goodness of fit p-value.
* <cell_type>_full_aic - Akaike’s information criteria.
* <cell_type>_full_bic - Bayes’ information criteria.
* <cell_type>_genotype_coeff - the genotype term coefficient.
* <cell_type>_genotype_cell_coeff - the coefficient of the genotype-cell frequencies interaction term. 
* <cell_type>_intercept_coeff.
