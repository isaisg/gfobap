#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $og_file=$ARGV[1];
my $faa_dir=$ARGV[0];
##Read all the genomes
#The idea is that we will have the complete repository of genes clustered for the particular analysis
opendir (my $dh,$faa_dir) || die "$faa_dir could not be opened\n";
my @faa_files=grep {/\.faa$/} readdir $dh;
closedir $dh;
my %gene2genome=();
my %genome2gene=();
foreach my $file (@faa_files)
{
	my $name=$file;
	$name=~s/\.faa//;
	my $fpath=join "/",$faa_dir,$file;
	open (IN,"<",$fpath) || die "$fpath could not be opened\n";
	while (<IN>)
	{
		chomp();
		next unless /^>/;
		my $line=$_;
		$line=~s/\>//;
		$gene2genome{$line}=$name;
		$genome2gene{$name}{$line}++;
	}
	close IN;
}


#Now read the phyla amphora markers data frame
my $pa_df="df_phyla_amphora_markers_3837_genomes.tsv";
my %scg=();
my %invscg=();
my %tax=();
open (IN,"<",$pa_df) || die "$pa_df could not be opened\n";
while (<IN>)
{
	chomp();
	my @ele=split "\t",$_;
	my $scgind=$ele[1];
	my $id=$ele[2];
	$scg{$scgind}{$id}++;
	$invscg{$id}=$scgind;
}
close IN;


#Read the Orthogroups file output by the Orthofinder pipeline
my %gene2cluster=();
my %cluster2gene=();
my $mcl_file=$og_file;
open (MCL,"<",$mcl_file) || die "$mcl_file could not be opened\n";
my $count=0;
while (<MCL>)
{
	$count++;
	my $cluster_id=join "","cluster_",$count;
	chomp();
	#Optimized to take the OrthoFinderOutput
	my @ele=split " ",$_;
	shift @ele;
	foreach my $gene (@ele)
	{
		$cluster2gene{$cluster_id}{$gene}++;
		$gene2cluster{$gene}=$cluster_id;
	}
}
close MCL;

#
foreach my $sc (keys %scg)
{
	my %total=();
	my %done_cluster=();
	my %found_cluster=();
	my %total_query=();
	my %missing=();
	foreach my $gene (keys $scg{$sc})
	{
		if (exists $gene2cluster{$gene})
		{
			$total{$gene}++;
			my $cluster_id=$gene2cluster{$gene};
			next if exists $done_cluster{$cluster_id};
			$done_cluster{$cluster_id}++;
			foreach my $qgene (keys $cluster2gene{$cluster_id})
			{
				$found_cluster{$cluster_id}{$qgene}++;
				$total_query{$qgene}++;
			}
		}
		else
		{
			if(exists $gene2genome{$gene})
			{
				#Here keep record of the genes that were not put into any Orthogroup
				$missing{$gene}++;
				#Sum another value to the total
				$total{$gene}++;	
			}
		}
	}
	
	my $total_amphora=scalar keys %total;
	my $total_found=scalar keys %total_query;
	my $total_missing=scalar keys %missing;
	next if $total_amphora == 0;
	#Print the general values
	if($total_found ==0)
	{
		print $sc,"\t",$total_amphora,"\t",0,"\t",$total_missing,"\t",0,"\t",0,"\n";
	}
	else
	{
	    	my $prop=$total_amphora/$total_found;
		my $part=scalar keys %found_cluster;
		print $sc,"\t",$total_amphora,"\t",$total_found,"\t",$total_missing,"\t",$prop,"\t",$part,"\n";
	}

}

