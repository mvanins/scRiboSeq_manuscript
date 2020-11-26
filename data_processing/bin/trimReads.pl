#!/usr/bin/env perl
# Michael VanInsberghe 2020-10-09
# Trims input read to a maximum length

use strict;
use warnings;
use Getopt::Long;

my $maxlength = (undef);

GetOptions(
  "maxlength=i" => \$maxlength
);

unless(defined $maxlength){
  die "Error: Please provide a length to trim to: $!\n";
}

my $count = 1;

while(<STDIN>){
  chomp;
  if($count == 1){
    print "$_\n";
    $count++;
    next;
  }
  if($count == 2){
    my $trimmed = substr $_, 0, $maxlength;
    print "$trimmed\n";
    $count++;
    next;
  }
  if($count == 3){
    print "$_\n";
    $count++;
    next;
  }
  if($count == 4){
    my $trimmed = substr $_, 0, $maxlength;
    print "$trimmed\n";
    $count = 1;
    next;
  }
}

