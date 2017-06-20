### Orthoinder Scripts using Diamond

Wrapper Scripts to run Orthofinder using Diamond a the search engine.These scripts were optimized to be used in the slurm based cluster. 
The scripts assume that OrthoFinder and Diamond are in your path.

These scripts were developed and tested using perl/5.18.2 , diamond v0.8.36.98 and Orthofinder-1.1.4


###### 1)

Given a folder with .faa files whose names is faa_dir. We run the first script of the pipeline

oh.orthofinderprepare.sh faa_dir - This shell script is calling orthofinder and stopping it after preparing the needed input files to run its algorithm.This script will create a folder with the name Results_Date, where Date is the Given month and day when the script was called. Inside this folder, there will be another folder named WorkingDirectory. WorkingDirectory folder contains the input .fa files plus some text files needed for Orthofinder to run properly.

###### 2)

The second script should be run inside the WorkingDirectory created by the first script.
oh.makediamond.sh - This script creates a diamond database per each .fa file present in the folder.

###### 3)
#Because orthofinder needs n^2 blast-like searches to run properly and we have taxons with aroun 600 genomes we implemented  a script that parallelizes the .faa diamond serches in the slurm cluster utlized.

oh.parallelize_diamond_orthofinder.pl - This script should be called inside the WorkingDirectory.

###### 4)

oh.orthofinderfromblast.sh WorkingDirectory - This script runs Orthofinder utlizing the diamond searches. The directory where the .fa files and the diamond searches should be passed as a parameter.
