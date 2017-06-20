#!/bin/bash

# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gdobap.
#
#    gdobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gdobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gdobap.  If not, see <http://www.gnu.org/licenses/>.

#This wrapper should be run in the folder where the input matrix is located
#Mofidy the spath  to wherever the enrichment scripts are located
spath=/home/isai/Documents/results/3837_pa_npa_soil_june_2016_resubmit/github_repo/gdobap/enrichment_tests;
file=$1;
#Run the hypergeometric test
echo "Running hypergeometric";
Rscript --vanilla "$spath"/oh.hypergeometric.R $file
echo "Running bugwas";
#Run bugwas
Rscript --vanilla "$spath"/oh.bugwas.R $file
#Create scoary input
echo "Creating Scoary input";
Rscript --vanilla "$spath"/oh.create_scoary_input.R $file;
echo "Running Scoary";
#Run scoary
scoary -t traits_scoary_pa_npa.tsv -g matrix_scoary_pa_npa.tsv -n tree_scoary_pa_npa.newick -s 2 -e 1000 -p 1.0  --threads 8
#Create phyloglm input binary and not binary
echo "Creating phyloglm input"
Rscript --vanilla "$spath"/oh.split_matrix_for_pagelphyloglm_binarize.R $file 1000000000000 binary
Rscript --vanilla "$spath"/oh.split_matrix_for_pagelphyloglm.R $file 10000000000000 raw  
#Copy the tree to the phyloglm folder
cp 3837_genomes_31scg_june2016.newick split_matrix_pa_npa_binary
cp 3837_genomes_31scg_june2016.newick split_matrix_pa_npa_raw
cp 3837_genomes_31scg_june2016.newick split_matrix_ra_soil_binary
cp 3837_genomes_31scg_june2016.newick split_matrix_ra_soil_raw
#Run phyloglm
echo "Running phyloglm";
#Change directory 
cd split_matrix_pa_npa_binary;
Rscript --vanilla "$spath"/oh.phyloglm.R splitted_file_1.tsv
cd ..;
cd split_matrix_pa_npa_raw;
Rscript --vanilla "$spath"/oh.phyloglm.R splitted_file_1.tsv
cd ..;
cd split_matrix_ra_soil_binary;
Rscript --vanilla "$spath"/oh.phyloglm.R splitted_file_1.tsv
cd ..;
cd split_matrix_ra_soil_raw;
Rscript --vanilla "$spath"/oh.phyloglm.R splitted_file_1.tsv
cd ..;
