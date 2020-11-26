#!/usr/bin/perl
# Michael VanInsberghe
# 2020-03-21
# Collapses identical sequences in input. Cluster names contain all input sequences
# collapseSequences.pl input1.fasta (input2.fasta input3.fasta ... )

use strict;
use warnings;
use File::Basename;

die("No fasta files supplied!\n $0 input_fasta_files.fa\n") unless $#ARGV >= 0;

open my $outgtf, ">>cluster_info.gtf" or die "Can't open output cluster_info.gtf: $!\n";

my %unique;

foreach my $fasta (@ARGV){
	open my $input, "$fasta" or die "Can't open input fasta $fasta: $!\n";

	my @sequence = ();
	
	while(<$input>){
		chomp;
		push @sequence, $_;
		
		if ($#sequence == 1){
			my ($name,$seq) = @sequence;
			@sequence = ();
			
			$name =~ s/^>//;
			$seq = uc $seq;
	
			if (exists $unique{$seq}){
				$unique{$seq} .= "_$name";
			}else{
				$unique{$seq} = "$name";
			}
		}
	}
	
	close $input;
}


my $cluster = 1;

foreach my $k (keys %unique){
	# print  ">$unique{$k}\n$k\n";
	print  ">tRNACluster\_$cluster\n$k\n";
	
	my $genename = $unique{$k};
	$genename =~ s/\(-\)/@/g;
	$genename =~ s/\(+\)/\$/g;
	$genename =~ s/-/./g;
	$genename =~ s/\(//g;
	$genename =~ s/\)//g;
	$genename =~ s/@/-/g;
	$genename =~ s/\$/+/g;

	my @line = ("tRNACluster\_$cluster",
				"tRNAScan",
				"gene",
				"1",
				length($k),
				".",
				"+",
				".");

	my $geneid = "tRNACluster\_$cluster";
	print $outgtf join("\t",@line[0..1]) . "\tgene\t" . join("\t",@line[3..7]) . "\tgene_id \"$geneid\"; gene_type \"mature_tRNA\"; gene_name \"$genename\"; level 2;\n";
	print $outgtf join("\t",@line[0..1]) . "\ttranscript\t" . join("\t",@line[3..7]) . "\tgene_id \"$geneid\"; transcript_id \"$geneid\_t\"; gene_type \"mature_tRNA\"; gene_name \"$genename\"; transcript_type \"mature_tRNA\"; transcript_name \"$genename\"; level 2;\n";
	print $outgtf join("\t",@line[0..1]) . "\texon\t" . join("\t",@line[3..7]) . "\tgene_id \"$geneid\"; transcript_id \"$geneid\_\"; gene_type \"mature_tRNA\"; gene_name \"$genename\"; transcript_type \"mature_tRNA\"; transcript_name \"$genename\"; level 2; exon_number 1;\n";
	
	# print $outgtf "tRNACluster\_$cluster\ttRNAScan\t1\t" . length($k) . "\t.\t+\t.\tgene_id \"$genename\"; gene_type \"mature_tRNA\"; gene_name \"$genename\"\n";

	$cluster++;
}
