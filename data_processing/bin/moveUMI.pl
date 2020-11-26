#!/usr/bin/perl
# Michael VanInsberghe
# move UMI sequence for solo input
# 2020-03-22

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($read1, $read2);
my ($out1, $out2);
my $length = 10;
my $offset = 0;

GetOptions(	"read1=s" => \$read1,
			"read2=s" => \$read2,
			"out1=s" => \$out1,
			"out2=s" => \$out2,
			"length=i" => \$length,
			"offset=i" => \$offset) or die ("Error in command line arguments");

die("No Fastq files supplied") unless ($read1 && $read2);
# my $length = $UMILength || 10;
# my $offset = $UMIOffset || 0;

my ($r1, $r2);
my ($R1, $R2);

if ($read1 =~ /.gz$/){
	open $r1, "zcat $read1 |" or die "Can't open gzipped read1 $read1: $!\n";
#	my $filebase = basename($read1, (".fastq.gz", "fq.gz"));
#	open $R1, " | gzip -c > $filebase\_moved.fastq.gz" or die "Can't open output gzip pipe $filebase: $!\n";
}else{
	open $r1, "$read1" or die "Can't open read1 $read1: $!\n";
}

if ($read2 =~ /.gz$/){
	open $r2, "zcat $read2 |" or die "Can't open gzipped read2 $read2: $!\n";
#	my $filebase = basename($read2, (".fastq.gz", "fq.gz"));
#	open $R2, " | gzip -c > $filebase\_moved.fastq.gz" or die "Can't open output gzip pipe $filebase: $!\n";
}else{
	open $r2, "$read2" or die "Can't open read2 $read2: $!\n";
}


if ($out1 =~ /.gz$/){
	open $R1, " | gzip -c > $out1" or die "Can't open output gzip pipe $out1: $!\n";
}else{
	open $R1, "$out1" or die "Can't open output read1 $out1: $!\n";
}

if ($out2 =~ /.gz$/){
	open $R2, " | gzip -c > $out2" or die "Can't open output gzip pipe $out2: $!\n";
}else{
	open $R2, "$out2" or die "Can't open output read1 $out2: $!\n";
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
	
	my $umis = substr $r1s, $offset, $length;
	my $umiq = substr $r1q, $offset, $length;

	my @r1h = split " ", $r1h;
	$r1h[0] .= "_$r2s\_$umis"; # cbc_UMI
	my @r2h = split " ", $r2h;
	$r2h[0] .= "_$r2s\_$umis"; # cbc_UMI

	my $o1s = substr $r1s, $length;
	my $o1q = substr $r1q, $length;

	my $o2s = $r2s . $umis; # cbc_UMI
	my $o2q = $r2q . $umiq; # cbc_UMI

	
	print $R1 join(" ", @r1h) . "\n$o1s\n+\n$o1q\n";
	print $R2 join(" ", @r2h) . "\n$o2s\n+\n$o2q\n";
	
	$endOfFiles = ( ($endOfFiles & !eof($r1)) & !eof($r2) );
}

close $r1;
close $r2;
close $R1;
close $R2;	
