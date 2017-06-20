### Orthoinder Scripts using Diamond

Wrapper Scripts to run Orthofinder using Diamond a the search engine.These scripts were optimized to be used in the slurm based cluster. 
The scripts assume that OrthoFinder and Diamond are in your path.

These scripts were developed and tested using perl/5.18.2 , diamond v0.8.36.98 and Orthofinder-1.1.4


##### 1)

Orthofinder outputs a file named Orthogroups.txt containing the proteins contained in any existen orthogroup. The next script appends a given prefix to each Orthogroup so we can merge multiple Orthogroup tables afterwards.

oh.add_orthogroups_prefix.pl prefix -The script receive as input the given prefix we wanna append to the Orthogorups. E.g oh.add_orthogroups_prefix.pl Xanthomonadaceae run in the WorkingDir of the Xanthomonadaceae will create an Orthogroups file with the Structure Xanthomondaceae_OG000000 ... etc.

##### 2)

oh.orthofinder_orthogroups2matrix.pl - This  script  creates a pangenome matrix from an Orthogroup.txt file output by orthofinder. This is achieved by mapping each given protein to a given genome. The script outputs two version of the pangenome matrix, one with the raw counts of orthogroups per genome and the other one just the binarize version of the raw matrix.
The script takes 3 arguments: 

-Folder where the .faa files analyzed are container
-Orthogroup.txt file that we wanna convert to matrix
-Prefix to attach to the name of the matrices output 

perl oh.orthofinder_orthogroups2matrix.pl faa_folder Orthogroup.txt prefix


##### 3)

oh.create_dataframe_from_orthogroups.pl - This scripts creates a dataframe with three columns (Orthogroup_Id,Genome,Gene_Id). These dataframes can be utlized to find which genes correspond to a given Orthogroups in a particular genome of interest. The script trakes 2 arguments

-Orthogroups.txt file that contains the orthogroups.
-Folder where the .faa files analyzed are container

perl oh.create_dataframe_from_orthogroups.pl Orthogroup.txt faa_folder


