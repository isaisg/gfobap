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
#!/usr/bin/perl
use strict;
use warnings;

my $orthogroups_file=$ARGV[0];
my $faa_folder=$ARGV[1];

my %orthogroups=();
open(IN,"<",$orthogroups_file) || die "$orthogroups_file could not be opened\n";
while (<IN>)
{
	chomp();
	my @ele=split " ",$_;
	$ele[0]=~s/\://g;
	for (my $i=1;$i < scalar @ele;$i++)
	{
		$orthogroups{$ele[0]}{$ele[$i]}++;
	}
}
close IN;

#Open the folder
opendir(my $dh,$faa_folder) || die "$faa_folder could not be opened\n";
my @files=grep {/\.faa/} readdir $dh;
closedir $dh;
my %gene2genome=();
foreach my $file (@files)
{
	my $name=$file;
	$name=~s/\.faa//;
	my $ffile=join "",$faa_folder,"/",$file;
	open (FILE,"<",$ffile) || die "$ffile could not be opened";
	while (<FILE>)
	{
		chomp();
		next unless /^>/;
		my $line=$_;
		$line=~s/\>//;
		$gene2genome{$line}=$name;
	}
	close FILE;
}
print "Orthogroup_Id\tGenome\tGene_Id\n";
foreach my $orthogroup (sort keys %orthogroups)
{
	foreach my $gene (sort keys $orthogroups{$orthogroup})
	{
		print $orthogroup,"\t",$gene2genome{$gene},"\t",$gene,"\n";
	}
}
