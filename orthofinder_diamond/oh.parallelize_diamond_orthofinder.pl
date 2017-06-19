# (C) Copyright 2016 Isai Salas Gonzalez
# #
# #    This file is part of gdobap.
# #
# #    gdobap is free software: you can redistribute it and/or modify
# #    it under the terms of the GNU General Public License as published by
# #    the Free Software Foundation, either version 3 of the License, or
# #    (at your option) any later version.
# #
# #    gdobap is distributed in the hope that it will be useful,
# #    but WITHOUT ANY WARRANTY; without even the implied warranty of
# #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# #    GNU General Public License for more details.
# #
# #    You should have received a copy of the GNU General Public License
# #    along with gdobap.  If not, see <http://www.gnu.org/licenses/>.
#
#!/usr/bin/perl
use strict;
use warnings;

##Unlink previous files 
my $del="rm batch*.txt";
system $del;

my $folder=".";

my $minijobs=600;
#my $batches=61343;
if (defined $ARGV[0])
{
	 $minijobs=$ARGV[0];
}

my $send="--very-sensitive";
if (defined $ARGV[1])
{
	$send=$ARGV[1];
}
#my $minijobs=$ARGV[0];
#my $minijobs=1200;
opendir (my $dh,$folder) || die "$folder could not be opened\n";
my @fa_files=grep {/\.fa$/} readdir $dh;
my @udb_files=grep {/\.udb$/} readdir $dh;
my @blast_files=grep {/Blast.*\.txt/} readdir $dh;
closedir $dh;

my %blasts=();
##Put the blast files in a hash
foreach my $blast (@blast_files)
{
	$blasts{$blast}++;
}


my %combinatorics=();
my $contador=0;
my $batch_num=1;

my %all_jobs=();

my $jobo=join "","batch_job.",$batch_num,".sh";
$all_jobs{$jobo}++;
open(OUT,">",$jobo) || die;
print OUT "#!/bin/bash\n";
print OUT "#SBATCH -n 4\n";
print OUT "#SBATCH -N 1\n";
print OUT "#SBATCH -t 2-23\n";
print OUT "#SBATCH --mem 10000\n";
#print OUT "#SBATCH -n 48\n";
#print OUT "#SBATCH -N 1\n";
#print OUT "#SBATCH --qos bigmem_access\n";
#print OUT "#SBATCH -p bigmem\n";
#print OUT "#SBATCH --mem=40000\n";
#print OUT "#SBATCH -t 2-23\n";
print OUT"\n";
print STDERR $batch_num,"\n";
foreach my $file (@fa_files)
{
	my $base=$file;
	$base=~s/\.fa//;
	foreach my $file2 (@fa_files)
	{
		if ($contador==$minijobs)
		{
			close OUT;
			$contador=0;
			$batch_num++;
			print STDERR $batch_num,"\n";
			$jobo=join "","batch_job.",$batch_num,".sh";
			$all_jobs{$jobo}++;
			open(OUT,">",$jobo) || die;
			print OUT "#!/bin/bash\n";
			print OUT "#SBATCH -n 4\n";
			print OUT "#SBATCH -N 1\n";
			print OUT "#SBATCH -t 2-23\n";
			print OUT "#SBATCH --mem 10000\n";
			#print OUT "#SBATCH -n 48\n";
			#print OUT "#SBATCH -N 1\n";
			#print OUT "#SBATCH --qos bigmem_access\n";
			#print OUT "#SBATCH -p bigmem\n";
			#print OUT "#SBATCH --mem=40000\n";
			#print OUT "#SBATCH -t 10-23\n";
			print OUT "\n";
		}
		$contador++;
		my $inbase=$file2;
		$inbase=~s/\.fa//;
		my $x=$base;
        	my $y=$inbase;
		$x=~s/Species//;
        	$y=~s/Species//;
		my $upath=join "",$inbase,".diamonddb";
		my $outfile=join "","Blast",$x,"_",$y,".txt";
		next if  exists $blasts{$outfile};
		next if exists $combinatorics{$base}{$inbase};
		print OUT "diamond blastp --threads 4 --db $upath --outfmt 6 --query $file --max-target-seqs 500 --evalue 0.001 $send --quiet --out $outfile\n";
		$combinatorics{$base}{$inbase}++;
			
	}
}



foreach my $j(keys %all_jobs)
{
	print STDERR "Running $j\n";
	my $command="sbatch $j";
	system $command;

}
