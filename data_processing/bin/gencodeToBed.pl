#!/usr/bin/perl
# Michael VanInsberghe
# 2020-03-20
# Convertts gencode GTF to a bed with useful names for bedtools
#
use warnings;
use strict;

die("No GTF supplied!\n$0 input.gtf\n") unless $#ARGV == 0;

open my $fh, "$ARGV[0]" or die "Can't open input gtf $ARGV[0]: $!\n";

while(<$fh>){
	chomp;
	next if /^#/;
	next if /^\s*$/;

	my @line = split "\t";

	my $gene_name = match_tag($line[8],qr/gene_name "(.*?)";/);
	my $gene_id = match_tag($line[8],qr/gene_id "(.*?)";/);
	# (my $gene_name) = $line[8] =~ m/gene_name "(.*?)";/;
	# (my $gene_id) = $line[8] =~ m/gene_id "(.*?)";/;

	my $name = $gene_id . "__" . $gene_name;
	my @out = ($line[0], $line[3]-1, $line[4], $name, 0, $line[6]);

	print join("\t", @out) . "\n";
}
close $fh;

sub match_tag {
	# to get around potential misassignment in direct match
    my $line = $_[0];
    my $pattern = $_[1];

    if ( $line =~ m/$pattern/ ){
        return $1;
    }else{
        return "";
    }
}

