### Script to benchmark the purity/integrity of a given set or orthogroups based on a comparison to PhylaAmphora markers

The data folder contains a set of dataframes, one per taxon and method utlized (Orthofinder/Ulcust) that depict the distribution of a particular PhylaAmphora Marker across a given set of orthogroups.
Asssuming a PhylaAmphora marker consists of true homologues, then if we scan that Marker across N given genomes and then we cluster all the proteins of those N genomes, we would expect that the Marker scanned should be contained in 1 single cluster of proteins with the exact number of proteins as found by the HMM scanning of the marker across the genomes.
Using this logic we devised a purity index and fragmentation index. 

The dataframes were constructed using the perl script -compareclusters2phylaamphora.pl and compareuclust2phylaamphora.pl . These scripts takes as input a directoey where PhylaAmphora was called for a set of N genomes. 
The results of these comparisons can be visualized using the R script phyla_amphora_of_uclust_comparison.R
