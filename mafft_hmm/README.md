### MAFFT alignments and HMM Profiles

The scripts contained in this folder are two perl scripts that parallelize the creation and submmission of jobs in a SLURM based cluster.
The scripts were tested using perl/5.18.2 , MAFFT v7.305b (2016/Aug/16) and HMMER 3.1b2 (February 2015).

oh.parallelize_mafft_orthofinder.pl - This Script takes a folder with fasta files corresponding to Orthogroups and creates alignments of each file. This corresponds to an alignment per Orthogroup.

oh.parallelize_mafft_orthofinder.pl - This Script takes a folder with alignment files corresponding to Orthogroups and create a HMM profile using hmmer-build.

The main intention of the two scripts described above  was to minimize the time needed to create all the alignments and HMM profiles for all the orthogroups reported in the article. 

The MAFFT command to compute the alignments was:

**mafft --anysymbol --thread 48  infile > outfile**

The HMMER command to create the HMM profile was:

**hmmbuild --amino --cpu 48 outname infile**
