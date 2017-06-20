### Scanner new genomes with the HMM profiles derived from the Orthogroups

The scripts in this folder help to annotane new sequenced genomes with th HMM profiles computed.

These scripts were developed/tested using perl/5.18.2 

##### 1)

profiler_hmm.pl - This script takes an aminoacid fasta file (new genome)  and a taxon number(see below) to map all the HMM corresponding to that taxon to the query genome.
The core function to scan the new genomes with an HMM are taken from Phyla_Amphora https://github.com/martinwu/Phyla_AMPHORA. At this point you need to a priori know to which of the 9 taxonomic categories your genome belong to.  
**Please edit this variable in the script my $general_folder, you should  copy the path where you downloaded the HMM profiles**. **The HMM profiles can be obtained in the following link**
Options:
	-taxon: 	1. Acinetobacter
			2. Actinobacteria_one
			3. Actinobacteria_two
			4. Alphaproteobacteria
			5. Bacillales
			6. Bacteroidetes
			7. Burkholderiales
			8. Pseudomonas
			9. Xanthomonadaceae
	-Evalue: HMMER evalue cutoff. Default: 1e-7 

##### 2)

create_matrix_from_hmm_tables.pl - This script will take the output of profiler_hmm.pl and create a pangenome-matrix of all the HMM that mapped to all the genomes  (faa files) present in the Working Directory.

##### 3) 

wrapper_profiler_create_matrix.sh - This bash script is a wrapper for scripts 1 and 2. If you have a vast number of faa files to scan to, it is recommended to parallelize profiler_hmm.pl. This wrapper script is just intended to show the functionality of the two previous commands.
**Please edit this variable in the script  dpath , you should copy the path where scripts 1 and 2 are located.**

#### Example

An example of usage can be found in test_scanner_orthogroups_alphaproteobacteria. In this folder 5 Alphaproteobacteria faa files are included.

The pipeline was called using the wrapper script. It was called this way

./wrapper_profiler_create_matrix.sh 4

The pipeline creates 3 folders: dir_candidates ,dir_hmmsearch and  dir_tables

dir_candidates contains the HMM hits found in the query genomes.

dir_hmmsearch contains lof files of the hmmsearch command.

dir_tables contain a datame per faa file (genome) depicting the gene_id -> HMM profile relationship. These dataframes are utilized to construct the pangenome matrix. The file matrix_hmm.tsv is the pangenome matrix created.
