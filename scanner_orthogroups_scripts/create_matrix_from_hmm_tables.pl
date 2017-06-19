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

my $table_folder="dir_tables";

opendir(my $dh,$table_folder) || die "$table_folder could not be opened\n";
my @files=grep {/^table_/} readdir $dh;
closedir $dh;

my %universe_hmm=();
my %matrix=();
#Read the table files
foreach my $file (@files)
{
	my $name=$file;
	$name=~s/table_query_to_hmm_//;
	$name=~s/\.faa//;
	my $path_file=join "/",$table_folder,$file;
	open (IN,"<",$path_file) || die "$path_file could not be opened\n";
	while (<IN>)
	{
		chomp();
		next if /^Gene_id/;
		my @ele=split "\t",$_;
		$matrix{$name}{$ele[1]}++;
		$universe_hmm{$ele[1]}++;
	}
	close IN;
}


#Print the matrix
my $outfile=join "/",$table_folder,"matrix_hmm.tsv";
open (OUT,">",$outfile) || die "$outfile could not be opened\n";
print OUT "Gene";
foreach my $genome (sort keys %matrix)
{
	print OUT "\t",$genome;
}
print OUT "\n";

foreach my $og (sort keys %universe_hmm)
{
	print OUT $og;
	foreach my $genome (sort keys %matrix)
	{
		if (exists $matrix{$genome}{$og})
		{
			print OUT "\t",$matrix{$genome}{$og};
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close OUT;
