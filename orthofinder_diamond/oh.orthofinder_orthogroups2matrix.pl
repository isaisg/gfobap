# (C) Copyright 2017 Isai Salas Gonzalez
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

#!/nas/longleaf/apps/perl/5.18.2/bin/perl
use strict;
use warnings;

my $faa_dir=$ARGV[0];
my $ortho_file=$ARGV[1];
my $prefix=$ARGV[2];
my $outfile=join "","matrix_raw_orthogroups_",$prefix,".tsv";
#Assumes that the faa file have a .fa termination (Such as the one requested by orthofinder)
opendir(my $dh,$faa_dir) || die "$faa_dir could not be opened\n";
my @fa_files=grep {/\.faa$/} readdir $dh;
closedir $dh;


#Read each fa file and keep the identifiers in a hash
my %genes2genome=();
#Create another hash to keep all the genomes
my %genomes=();
foreach my $file (@fa_files)
{
	my $genome_name=$file;
	$genome_name=~s/\.faa//;
	my $fpath=join "",$faa_dir,"/",$file;
	open (IN,"<",$fpath) || die "$fpath could not be opened\n";
	while (<IN>)
	{
		chomp();
		next unless /^>/;
		my $line=$_;
		$line=~s/\>//;
		$genes2genome{$line}=$genome_name;
		$genomes{$genome_name}++;
	}
	close IN;
}



###Now read the orthofiles
my %ortho2genome=();
open (OR,"<",$ortho_file) || die "$ortho_file could not be opened\n";
while (<OR>)
{
	chomp();
	my @ele=split " ",$_;
	$ele[0]=~s/\://;
	for (my $i=1;$i < scalar @ele;$i++)
	{
		my $gene=$ele[$i];
		my $genome=$genes2genome{$gene};
		$ortho2genome{$ele[0]}{$genome}++;
	}
}
close OR;

###PRint the matrix
open (OUT,">",$outfile) || die "$outfile could not be opened\n";
print OUT "Gene";
foreach my $genome (sort {$a<=>$b} keys %genomes)
{
	print OUT "\t",$genome;
}
print OUT "\n";

foreach my $og (sort  keys %ortho2genome)
{
	print OUT $og;
	foreach my $genome (sort {$a <=>$b} keys %genomes)
	{
		if (exists $ortho2genome{$og}{$genome})
		{
			print OUT "\t",$ortho2genome{$og}{$genome};
		}
		else
		{
			print OUT "\t0";
		}
	}
	print OUT "\n";
}
close OUT;

###Input the binarize version of the same matrix
my $outfile_binary=join "","matrix_binary_orthogroups_",$prefix,".tsv";
open (OUT,">",$outfile_binary) || die "$outfile_binary could not be opened\n";
print OUT "Gene";
foreach my $genome (sort {$a<=>$b} keys %genomes)
{
        print OUT "\t",$genome;
}
print OUT "\n";

foreach my $og (sort  keys %ortho2genome)
{
        print OUT $og;
        foreach my $genome (sort {$a <=>$b} keys %genomes)
        {
                if (exists $ortho2genome{$og}{$genome})
                {
                        print OUT "\t1";
                }
                else
                {
                        print OUT "\t0";
                }
        }
        print OUT "\n";
}
close OUT;

