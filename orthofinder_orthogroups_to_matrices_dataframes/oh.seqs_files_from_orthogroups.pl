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
use Bio::SeqIO;

my $faa_dir=$ARGV[0];
my $orthogroups_file=$ARGV[1];

#Create a folder to contain the sequences
mkdir "orthogroups_sequences" unless -d "orthogroups_sequences";

opendir(my $dh,$faa_dir) || die "$faa_dir could not be opened\n";
my @fa_files=grep {/\.faa$/} readdir $dh;
closedir $dh;

#Read all the faa files and create a hash with the id and sequence
my %genes2seq=();
foreach my $file (@fa_files)
{
	my $genome_name=$file;
	$genome_name=~s/\.faa//;
	my $fpath=join "",$faa_dir,"/",$file;
	my $in  = Bio::SeqIO->new(-file => $fpath ,-format => 'Fasta');
	while ( my $Seq = $in->next_seq() )
	{
		my $id=$Seq -> id();
		my $seq= $Seq -> seq();
		$genes2seq{$id}=$seq;
	}
}



#Read the orthogroups file and output a fasta sequence fo each orthogroup with all the members of the orthogroup
open (IN,"<",$orthogroups_file) || die "$orthogroups_file could not be opened\n";
while (<IN>)
{
	chomp();
	my @ele=split " ",$_;
	my $id=$ele[0];
	$id=~s/\://;
	my $outfile=join "","orthogroups_sequences/",$id,".fasta";
	open (OUT,">",$outfile ) || die "$outfile could not be opened\n";
	for (my $i=1;$i < scalar @ele;$i++)
	{
		my $hit=$ele[$i];
		my $seq=$genes2seq{$hit};	
		print OUT ">$hit\n$seq\n";
	}
	close OUT;
}
close IN;
