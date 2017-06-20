### Compute enrichment/depletion  profiles based on a pangenome matrix, metadata and a phylogenetic tree

These scripts are the core of the enrichment/depletion calls in the 9 taxons analyzed. 

These scripts were developed/tests using r/3.3.1 , Scoary version 1.6.10 , and python/3.5.1
The next R packages with their corresponding version are required

phylolm_2.5  , bugwas_0.0.0.9000 , phytools_0.5-64  , maps_3.1.1   ,ape_4.1 

 
#### Hypergeometric test

oh.hypergeometric.R - R script that takes as input a pangenome matrix and returns two tables. One of enrichment betweens PA and NPA and the other between RA and Soil genomes. The tables include two versions of the Hypergeometric test, one in which the values in the matrix are binarized and one in which they are not. The script requires the metadata file located in  . It assumes the metadata file is in the . folder.


##### PhyloGLM

Because the vast number of orthogroups we parallelized  running Phyloglm. 

oh.split_matrix_for_pagelphyloglm.R - R script that takes a the metadata file contained in and rearrange a provided pangenome matrix. This script creates the matrix format required  by PhyloGLM. The script automatically process PA/NPA and RA/Soil comparisons. 

oh.split_matrix_for_pagelphyloglm_binarize.R - R script that is equivalent to  oh.split_matrix_for_pagelphyloglm_binarize.R but instead of using raw values, it binarizes the provided pangenome matrix.


oh.phyloglm.R - R script that runs PhyloGLM over a prepared matrix output either by oh.split_matrix_for_pagelphyloglm.R or oh.split_matrix_for_pagelphyloglm_binarize.R. The script uses the phylogenetic tree proided in. 


oh.merge_phyloglm.R - Rscript to concatenate multiple runs of Phyloglm. 

#### Scoary

oh.create_scoary_input.R - R Script that takes a pangenome matrix, the metadata file contained in and the phylogenetic tree contained in and creates three files per comparison (PA/NPA and RA/Soil). The files created are traits_scoary_pa_npa.tsv matrix_scoary_pa_npa.tsv, tree_scoary_pa_npa.newick, traits_scoary_ra_soil.tsv, matrix_scoary_ra_soil.tsv and tree_scoary_ra_soil.newick. These files are the input to scoary.


oh.scoary.sh - Shell script that runs scoary taking the three files created per comparison. E.g oh.scoary.sh traits_scoary_pa_npa.tsv matrix_scoary_pa_npa.tsv tree_scoary_pa_npa.newick

#### Bugwas

oh.bugwas.R - RScript that takes a pangenome matrix and run bugwas. The analysis is automatically run for PA/NPa and RA/Soil comparison.

### Wrapper

oh.wrapper_enrichment.sh - This wrapper exemplifies the calling of all the scripts. Edit the script to include the path where all the other enrichment files are located.  And example of the functionality of all the scripts is presented in the subfolder test_enrichments.
### Example

In the test_enrichments file we have the pangenome matrix matrix_Xanthomonadaceae_cog.tsv . As well, in the same folder we have the metadata file available at and the phylogenetic tree available at.

To compute all the enrichments we simply call the wrapper:

oh.wrapper_enrichment.sh matrix_Xanthomonadaceae_cog.tsv

