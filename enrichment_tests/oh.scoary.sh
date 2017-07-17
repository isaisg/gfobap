#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 2-1
#SBATCH --mem 80000

# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gfobap.
#
#    gfobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gfobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gfobap.  If not, see <http://www.gnu.org/licenses/>.


trait=$1;
gene=$2;
tree=$3;
scoary -t $trait -g $gene -n $tree -s 2 -e 1000 -p 1.0  --threads 8
