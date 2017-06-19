# (C) Copyright 2016 Isai Salas Gonzalez
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

#!/usr/bin/bash

dpath="/pine/scr/i/s/isai/3837_march_2017/archive_publication_files/scanner_orthogroups";
option=$1;
eval=1e-7;

if [ $# -eq 0 ]
  then
    echo "One argument is mandatory Taxon to choose";
    echo "If second argument is provided then the eval threshold is modified. Default 1e-7";
    echo "Example chossing Bacillales but keeping the eval 1e-7: ./wrapper_profiler_create_matrix.pl 5";
    echo "Example modifying evalue: ./wrapper_profiler_create_matrix.pl 5 1e-10";
    echo "Options:
	-taxon: 	1. Acinetobacter
			2. Actinobacteria_one
			3. Actinobacteria_two
			4. Alphaproteobacteria
			5. Bacillales
			6. Bacteroidetes
			7. Burkholderiales
			8. Pseudomonas
			9. Xanthomonadaceae
	"
	exit 1;
fi
if [ $# -eq 2 ]
	then
		eval=$2;
fi

#echo $option;
#echo $eval;

##Check if there are faa files in the folder
tot=`find -maxdepth 1 -name "*.faa" | wc -l`;
#echo $tot;
if [ $tot == 0 ]
	then
		echo "Folder without .faa files present. Please check the input";
		exit 1;
fi;

ls *.faa | while read file; do 
echo "Working on $file";
perl $dpath/profiler_hmm.pl -taxon $option -Evalue $eval $file;
done
echo "Creating matrix";
perl $dpath/create_matrix_from_hmm_tables.pl 
