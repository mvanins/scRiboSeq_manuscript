#!/usr/bin/perl
# Michael VanInsberghe
# 2020-03-24
# Extracts unique gene name, type from gtf
# outputs as tsv
# extractGeneInfo.pl <gtf> 
#

use warnings;
use strict;

die("No GTF supplied!\n$0 input.gtf\n") unless $#ARGV == 0;

open my $fh, "$ARGV[0]" or die "Can't open input GTF $ARGV[0]: $!\n";

my %genes;
while(<$fh>){
	chomp;
	next if /^#/;
	next if /^\s*$/;

	my @line = split "\t";
	next unless $line[2] eq "gene";
	my $geneinfo = $line[8];

	# my $havana_gene = match_tag($line[8], qr/havana_gene "(.*?)";/);
	my $gene_id = match_tag($line[8], qr/gene_id "(.*?)";/);
	my $gene_name = match_tag($line[8], qr/gene_name "(.*?)";/);
	my $gene_type = match_tag($line[8], qr/gene_type "(.*?)";/);
	my $source = $line[1];
	my $level = match_tag($line[8], qr/level (\d*?);/);
	my $havana_gene = match_tag($line[8], qr/havana_gene "(.*?)";/);
	my $hgnc_id = match_tag($line[8], qr/hgnc_id "(.*?)";/);
	my $tag = match_tag($line[8], qr/tag "(.*?)";/);

	$genes{$gene_id} = join "\t", ($gene_id,$gene_name,$gene_type,$source,$level,$havana_gene,$hgnc_id,$tag);

	# print "$genes{$gene_id}\n";
}


print join("\t", ("gene_id","gene_name","gene_type","source","level","havana_gene","hgnc_id","tag")) . "\n";
foreach my $gi (keys %genes){
	print "$genes{$gi}\n";
}

sub match_tag {
	my $line = $_[0];
	my $pattern = $_[1];

	if ( $line =~ m/$pattern/ ){
		return $1;
	}else{
		return "";
	}
}
