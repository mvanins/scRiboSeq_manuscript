#!/usr/bin/perl
# Michael VanInsberghe
# 2021-02-28
# symmetrically expands BED12 file similar to bedtools slop, but also expanding block regions
#

use strict;
use warnings;
use Getopt::Long;

my $slop;
my $index;

GetOptions("slop=i" => \$slop,
           "index=s" => \$index);

die("No expansion provided in --slop") unless($slop);
die("No genome index provided in --index") unless($index);

# read genome
my %genomeSize;

open my $idx, "$index" or die "Can't open genome index file $index: $!\n";
while(<$idx>){
  chomp;
  next if /^\s*$/;

  my ($chr, $size, @temp) = split '\s';
  $genomeSize{$chr} = $size;
}
close $idx;

while(<STDIN>){
  chomp;
  next if /^\s*$/;

  my ($chr, $start, $stop, $name, $score, $strand, $thickStart, $thickEnd, $rgb, $blockCount, $blockSizes, $blockStarts) = split '\t';
  $blockSizes =~ s/,$//;  # remove trailing , if exists
  $blockStarts =~ s/,$//;
  
  my $lSlop = correctSlop(\%genomeSize, $chr, $start, -1 * $slop);
  my $rSlop = correctSlop(\%genomeSize, $chr, $stop, $slop);

    print join("\t",
      $chr, 
      $start + $lSlop, 
      $stop + $rSlop, 
      $name,
      $score, 
      $strand, 
      $thickStart+$lSlop, 
      $thickEnd+$rSlop, 
      $rgb, 
      $blockCount) . "\t";

  if($blockCount == 1){
    # if no exon
    print join("\t",
      $blockSizes - $lSlop + $rSlop, 
      $blockStarts) . "\n";
  }else{
    my @sizes = split(',', $blockSizes);
    my @starts = split(',', $blockStarts);

    $sizes[0] = $sizes[0] - $lSlop;
    $sizes[-1] = $sizes[-1] + $rSlop;

    $starts[-1] = $starts[-1] - $lSlop;

    print join("\t",
      join(',', @sizes),
      join(',', @starts) ) . "\n";
  }

}

sub correctSlop {
  # checks that the slop doesn't spill past chromosome boundaries
  my $genomeSize_ref = $_[0];
  my $chr = $_[1];
  my $position = $_[2];
  my $change = $_[3];

  if($change < 0){
    # change is subtract, check if result is less than 1
    if($position + $change > 1){
      return $change;
    }else{
      return 1 - $position;
    }
  }

  if($change > 0){
    # change is addition, check if result is greater than chromosome boundary
    if($position + $change <= $genomeSize_ref->{$chr}){
      return $change;
    }else{
      return $genomeSize_ref->{$chr} - $position;
    }
  }
}

