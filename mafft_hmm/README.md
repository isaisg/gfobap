###MAFFT alignments and HMM Profiles #######
The scripts contained in this folder are two perl scripts that parallelize the creation and submmission of jobs in a SLURM based cluster.

oh.parallelize_mafft_orthofinder.pl - This Script takes a folder with fasta files corresponding to Orthogroups and creates alignments of each file. This corresponds to an alignment per Orthogroup.

oh.parallelize_mafft_orthofinder.pl - This Script takes a folder with alignment files corresponding to Orthogroups and create a HMM profile using hmmer-build 

