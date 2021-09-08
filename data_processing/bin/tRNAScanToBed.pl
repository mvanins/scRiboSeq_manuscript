#!/usr/bin/perl
# Michael VanInsberghe
# 2021-02-28
# Converts tRNAScan-SE output to bed12 for downstream processing
# Default bed from tRNAScan-SE doesn't contain information on which hits are pseudogenes 
# Based on tRNAscan2bed12.pl from tRNA sequencing best practice workflow

use strict;
use warnings;

while(<STDIN>){
  chomp;
  next if /^\s*$/;

  my @hit = split /\s+/;

  my $chr = $hit[0];
  my $chrHitCount = $hit[1];
  my $left = $hit[2];
  my $right = $hit[3];
  my $three = $hit[4];
  my $codon = $hit[5];
  my $intronStart = $hit[6];
  my $intronEnd = $hit[7];
  my $score = int($hit[8]*10);
  $score = ($score > 1000)?(1000):($score);
  
  my $name = "$chr.tRNA$chrHitCount-$three$codon";
  $name = $name . "-$hit[11]" if($hit[11]);

  my $strand;
  my $start;
  my $stop;
  if($right > $left){
    $strand = "+";
    $start = $left - 1;
    $stop = $right;
  }else{
    $strand = "-";
    $start = $right - 1;
    $stop = $left;
  }

  my $blockCount = 1;
  my $blockSize = abs($stop - $start);
  my $blockStart = 0;
  if($intronStart && $intronEnd){
    $blockCount = 2;

    if($strand eq "+"){
      $blockSize = join(',', $intronStart - $left, $right - $intronEnd);
      $blockStart = join(',', 0, $intronEnd - $left + 1);
    }elsif($strand eq "-"){
      $blockSize = join(',', $intronEnd - $right, $left - $intronStart);
      $blockStart = join(',', 0, $intronStart - $right);
    }
  }

  print join("\t", $chr, $start, $stop, $name, $score, $strand, $start, $stop, 0, $blockCount, $blockSize, $blockStart) . "\n";

}

