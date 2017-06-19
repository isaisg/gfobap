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
#!/usr/bin/perl
use strict;
use warnings;

my $dir=".";

my $minijobs=600;
if (defined $ARGV[0])
{
         $minijobs=$ARGV[0];
}

#Read the dir

opendir (my $dh,$dir) || die "$dir could not be opened\n";
my @files=grep {/\.fasta$/} readdir $dh;
closedir $dh;



my $contador=0;
my $batch_num=1;

my %all_jobs=();

my $jobo=join "","batch_job.",$batch_num,".sh";
$all_jobs{$jobo}++;
open(OUT,">",$jobo) || die;
print OUT "#!/bin/bash\n";
print OUT "#SBATCH -n 48\n";
print OUT "#SBATCH -N 1\n";
print OUT "#SBATCH -t 10-10\n";
print OUT "#SBATCH --mem 80000\n";
print OUT"\n";
print STDERR $batch_num,"\n";


foreach my $file (@files)
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
		print OUT "#SBATCH -n 48\n";
		print OUT "#SBATCH -N 1\n";
		print OUT "#SBATCH -t 10-10\n";
		print OUT "#SBATCH --mem 80000\n";
		print OUT"\n";

	}
	my $num=`/usr/bin/grep -c ">" $file`;
	chomp($num);
	next if $num ==1;
	$contador++;
	my $outname=$file;
	$outname=~s/\.fasta/\.hmm/;
	$file=~s/\.fasta/\.mafft/;
	#print $file,"\t",$num,"\n";
	my $command="hmmbuild --amino --cpu 48 $outname $file";
	print OUT $command,"\n";
}



foreach my $j(keys %all_jobs)
{
	print STDERR "Running $j\n";
	#my $command="sbatch $j";
	#system $command;

}

