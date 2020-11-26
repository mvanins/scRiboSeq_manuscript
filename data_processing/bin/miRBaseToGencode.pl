#!/usr/bin/perl
# Michael VanInsberghe 2020-02-20
# Converts miRBase v22 GFF3 to GTF format similar to Gencodev33 annotations
# outputs mature miRNAs
# cat miRBase.gff | miRBaseToGencode.pl > miRBase.gtf
#

use warnings;
use strict;

# my $mirbase = shift; # miRBase GFF3

#open my $input, "$mirbase" or die "Can't open input miRBase GFF3 $mirbase: $!\n";
while(<STDIN>){
	#chomp;
	
	# skip gff version declaration
	next if /^##gff/;
	
	if (/^##date/){
		print "$_";
		next;
	}
	
	# if header, add another hash
	if (/^#\s/){
		print "#$_";
		next;
	}
	
	chomp;	
	my @line = split "\t";
	
	# skip miRNA_primary_transcript
	next if $line[2] eq "miRNA_primary_transcript";
	
	# add mirbase as source
	$line[1] = "miRBase";
	
	(my $id, my $alias, my $name, my $derives) = $line[8] =~ m/ID=(\w+);Alias=(\w+);Name=(.+);Derives_from=(\w+)/g;
	#print "$id\t$alias\t$name\t$derives\n";
	
	my $geneDescription = "gene_id \"$id\"; gene_type \"mature_miRNA\"; gene_name \"$name\"; level 2; derives_from \"$derives\";";
	my $transcriptDescription = "gene_id \"$id\"; transcript_id \"$id\"; gene_type \"mature_miRNA\"; gene_name \"$name\"; level 2; derives_from \"$derives\"; transcript_type \"mature_miRNA\"; transcript_name \"$name\";";
	my $exonDescription = $transcriptDescription . " exon_number 1;";
	
		
	print join("\t",@line[0..1]) . "\tgene\t" . join("\t",@line[3..7]) . "\t$geneDescription\n";
	print join("\t",@line[0..1]) . "\ttranscript\t" . join("\t",@line[3..7]) . "\t$transcriptDescription\n";
	print join("\t",@line[0..1]) . "\texon\t" . join("\t",@line[3..7]) . "\t$exonDescription\n";

}

