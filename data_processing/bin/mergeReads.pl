#!/usr/bin/perl
# Michael VanInsberghe 2020-03-31
# merges reads from two fastq files
# uses "read1" read names

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($read1, $read2);
my ($out);

GetOptions(	"read1=s" => \$read1,
			"read2=s" => \$read2) or die("Error in command line arguments");

die("No fastq files supplied") unless ($read1 && $read2);

my ($r1, $r2);
my ($R1, $R2);

if ($read1 =~ /.gz$/){
    open $r1, "zcat $read1 |" or die "Can't open gzipped read1 $read1: $!\n";
}else{
	open $r1, "$read1" or die "Can't open read1 $read1: $!\n";
}

if ($read2 =~ /.gz$/){
	open $r2, "zcat $read2 |" or die "Can't open gzipped read2 $read2: $!\n";
}else{
	open $r2, "$read2" or die "Can't open read2 $read2: $!\n";
}

my $endOfFiles = 1;

my ($r1h, $r1s, $r1p, $r1q);
my ($r2h, $r2s, $r2p, $r2q);
while($endOfFiles){
	chomp($r1h = readline($r1));
	chomp($r1s = readline($r1));
	chomp($r1p = readline($r1));
	chomp($r1q = readline($r1));

	chomp($r2h = readline($r2));
	chomp($r2s = readline($r2));
	chomp($r2p = readline($r2));
	chomp($r2q = readline($r2));

	my $sequence = $r1s . $r2s;
	my $quality = $r1q . $r2q;

	print "$r1h\n$sequence\n+\n$quality\n";

	$endOfFiles = ( ($endOfFiles & !eof($r1)) & !eof($r2) );	
}

