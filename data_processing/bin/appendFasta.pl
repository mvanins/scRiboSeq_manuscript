#!/usr/bin/perl
# Michael VanInsberghe
# 2021-03-01
# append or prepend given sequence to input fasta
#

use strict;
use warnings;
use Getopt::Long;

my $prepend = '';
my $append = '';
my $input;

GetOptions("prepend=s" => \$prepend,
           "append=s" => \$append);

die("No prepend/append sequence provided in one of --prepend, --append") unless( $prepend || $append);

my $name;
my $sequence;

while(<STDIN>){
  chomp;
  next if /^\s*$/;

  if(/^>/){
    print(">$name\n$prepend$sequence$append\n") if($name);
    ($name = $_) =~ s/^>//;
    $sequence = '';
    next;
  }else{
    $sequence .= $_;
    next;
  }
}

