#!/usr/bin/perl
# Michael VanInsberghe
# 2020-03-20
# merge gtf files, keeping both headers, sort output
#

use strict;
use warnings;

die("No gtf files supplied") unless $#ARGV >= 0;

my $header;
my $body;

foreach my $f (@ARGV){
	open my $fh, "$f" or die "Can't open gtf $f: $!\n";
	while(<$fh>){
		chomp;
		if (/^#/){
			$header .= "$_\n";
		}else{
			$body .= "$_\n";
		}
	}
	close $fh;
}

print STDOUT "$header";

open SP, " | sort -k1,1V -k4,4n -k5,5rn -k3,3r " or die "Can't open sort output pipe: $!\n";
print SP "$body";
close SP;
