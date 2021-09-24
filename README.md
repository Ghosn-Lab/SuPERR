# SuPERR-Seq
These files are intended to support our manuscript titled "Comprehensive multi-omics single-cell data integration reveals greater heterogeneity in the human immune system".
This work can be found freely available at https://www.biorxiv.org/content/10.1101/2021.07.25.453651v2

# Manual Gating of single-cell sequencing data
In the manuscript "Comprehensive multi-omics single-cell data integration reveals greater heterogeneity in the human immune system" we refer to Manual Gating using a customized strategy of biaxal plots, which was implemented in MATLAB.
Following the code contained in 'Anti_seq_manual_gating_BM.m' and 'Anti_seq_manual_gating_PBMC.m', we can easily reproduce the polygons used to produce the gates, and therefore reproduce the results of the manuscript.


In the manuscript "Comprehensive multi-omics single-cell data integration reveals greater heterogeneity in the human immune system" we refer to a Cluster Mismatch Statistic (CMS). CMS is applied exactly as in our recent publication "Data Matrix Normalization and Merging Strategies Minimize Batch-specific Systemic Variation in scRNA-Seq Data", which is freely available at https://www.biorxiv.org/content/10.1101/2021.08.18.456898v1

This document is intended to guide the user through the steps needed to reproduce the CMS calculation.

# CMS

The Cell Misclassification Statistic (CMS) is a scoring metric we apply to evaluate whether cells in a single-cell RNA-seq experiment retain the same classification variable (e.g. cell type, function, state) between replicate dataset analyses.
CMS is presented here exactly as in our complementary work "Data Matrix Normalization and Merging Strategies Minimize Batch-specific Systemic Variation in scRNA-Seq Data" which is freely available at https://www.biorxiv.org/content/10.1101/2021.08.18.456898v1.

## Software Requirements

To calculate a CMS score we need only an up-to-date version of base R (tested on R 4.0.5 running on MacOS 11.5)
No special installation is required.

## Example

To calculate a CMS we will need a dataframe "x" containing the cell barcodes and the classifications which we wish to test (one row for each cell).
We provide a small example using demo data here:
The expected runtime for this demo is less than one second.


> head(x)

#>                     	ID_1       		ID_2    
#> AAACCTGCAAGCGAGT		"B_Cells"		"B_Cells"  
#> AAACCTGCACACAGAG		"NK_CD56Hi"		"NK_CD56Hi"
#> AAACCTGGTAAACACA		"NK"       		"NK"       
#> AAACCTGGTCGGATCC		"T_CD8"    		"NK"    
#> AAACCTGGTCTCTTTA		"T_CD4"    		"T_CD4"    
#> AAACCTGGTTTAAGCC		"NK"       		"NK"
#> AAACCTGTCAACACCA		"NK"       		"NK"                
#> AAACCTGTCTATCCCG		"Monocyte" 		"Monocyte"

# Find which IDs from workflow 1 match IDs from workflow 2
> match <- x[, 'ID_1'] == x[, 'ID_2']

# Divide the number of matches by the total number of cells
> CMS <- sum(match, na.rm = T) / length(match)

# Invert the score so that 1 is 100% mismatch and 0 is 100% match
> CMS <- 1 - CMS
> print(CMS)
#> 0.125

In this example, 7 of 8 cells (87.5%) change classification between analyses (represented by ID_1 and ID_2).
The calculated CMS is 0.125, or to say that a 12.5% rate of misclassification is occurring between analyses.


CMS can be compared between analyses to evaluate the cell classification fidelity of a given dataset analysis workflow.

