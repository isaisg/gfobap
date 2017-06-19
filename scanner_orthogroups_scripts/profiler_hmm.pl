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
use Bio::SeqIO;
use Getopt::Long;

#Script to annotate a folder of .faa files using the a folder of hmm profiles
#$general_folder is the only variable that needs to be modified in the script.
my $general_folder="/pine/scr/i/s/isai/3837_march_2017/archive_publication_files/scanner_orthogroups";

my $db_folder=join "",$general_folder,"/hmm_databases";

my $evalue = 1e-7;
my $option_taxon=undef;
my $help = undef;
my %seq = ();


my $usage = qq~
Usage: $0 -taxon [1-9] sequence-file

Options:
	-taxon: 	1. Acinetobacter
			2. Actinobacteria_one
			3. Actinobacteria_two
			4. Alphaproteobacteria
			5. Bacillales
			6. Bacteroidetes
			7. Burkholderiales
			8. Pseudomonas
			9. Xanthomonadaceae
	-Evalue: HMMER evalue cutoff. Default: 1e-7 
	-Help: print help;
~;



GetOptions (	
		'Evalue=f'=>\$evalue,
		'taxon=i'=>\$option_taxon,
		'Help'=>\$help) || die "Invalid command line options\n";

die $usage if $help;
die $usage unless $ARGV[0];
die $usage unless $option_taxon;

my %number2taxon=(
 	"1" => "acinetobacter",
	"2" => "actinobacteria_one",
	"3" => "actinobacteria_two",
	"4" => "alphaproteobacteria",
	"5" => "bacillales",
	"6" => "bacteroidetes",
	"7" => "burkholderiales",
	"8" => "pseudomonas",
	"9" => "xanthomonadaceae"	
);

#File to work over 
my $input_seq= $ARGV[0];
#Select database according to the option provided
my $chosen_taxon=$number2taxon{$option_taxon};

my $hmm_database=join "",$db_folder,"/orthogroups_database_",$chosen_taxon,".hmm";

#Create two folders
my $dirt="dir_tables";
my $dirc="dir_candidates";
my $dirh="dir_hmmsearch";

mkdir $dirh unless -d $dirh;
mkdir $dirt unless -d $dirt;
mkdir $dirc unless -d $dirc;

#Get candidates
get_candidates();


#####################################################################################################################
#Subroutine taken from Phyla amphora        
sub get_candidates {
	#Add and inverted hits hash to contain the query-hmm relationtship
        my (%score, %candidates, %query_length, %hmm_length, %query_match, %hmm_match, %hits,%inverted_hits) =();
        #Remove the -Z option
	system ("hmmsearch --cpu 48 -E $evalue --domE $evalue --domtblout $dirh/$$.hmmsearch -o /dev/null $hmm_database $input_seq");         # fix the number of sequences in the database for E-value calculation

        open (IN, "$dirh/$$.hmmsearch") || die "Can't open $dirh/$$.hmmsearch";
        while (<IN>) {
                chop;
                next if /^#/;
                my ($query, $query_accession, $qlength, $hmm, $hmm_accession, $hmm_length, $evalue, $score, $bias, $domain, $domain_number, $dom_evalue, $ievalue, $dom_score, $dom_bias, $hmm_start, $hmm_stop,$query_start, $query_stop, $rest) = split /\s+/;
                if ((! exists $score{$query}) or ($score{$query} < $score)) {
                        $score{$query} = $score;
                        $candidates{$query} = $hmm;
                }
                $query_length{$query} = $qlength;
                $hmm_length{$hmm} = $hmm_length;
                $query_match{$hmm}{$query} += ($query_stop - $query_start);
                $hmm_match{$hmm}{$query} += ($hmm_stop - $hmm_start);
        }
        close IN;
        
        while ( my ($query, $hmm) = each (%candidates) ) {
                next if ( ($query_match{$hmm}{$query}/$query_length{$query} < 0.7) and ($hmm_match{$hmm}
{$query}/$hmm_length{$hmm} < 0.7) );    #ignore the hit if the match is partial
                $hits{$hmm}{$query} = 1;
		$inverted_hits{$query}=$hmm;
        }

        my $seqin = new Bio::SeqIO('-file'=>$input_seq);
        while (my $seq = $seqin->next_seq) {
                $seq{$seq->id} = $seq;
        }
        
        for my $marker (keys %hits) {
                my $seqout = new Bio::SeqIO('-file'=>">>$dirc/$$.$marker.candidate.pep",'-format'=>'fasta');
                for my $seqid (keys %{$hits{$marker}}) {
                        $seqout->write_seq($seq{$seqid});
                }
        }
	#Print a file containing the query -hmm relationship
	my $outtable=join "",$dirt,"/table_query_to_hmm_",$input_seq;
	open(OUT,">",$outtable) || die "$outtable could not be opened\n";
	print OUT "Gene_id\tHMM_hit\n";
	for my $q (keys %inverted_hits)
	{
		print OUT $q,"\t",$inverted_hits{$q},"\n";
	}
	close OUT;
        #return \%hits;
}


