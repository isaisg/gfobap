#!/bin/bash
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -t 10-10
#SBATCH --mem 200g
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


dir=$1;
orthofinder -b $dir -a 24 -og
