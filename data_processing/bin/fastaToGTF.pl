#!/usr/bin/env perl
# Michael VanInsberghe 2020-04-05
# Takes an input fasta files and creates a GTF where each entry covers the fasta
#

use strict;
use warnings;
use File::Basename;

die("No fasta files supplied!\n $0 input_fasta_file1.fa input_fasta_file2.fa\n") unless $#ARGV >=0;

#open my $outgtf, ">fasta_info.gtf" or die "Can't open output fasta_info.gtf: $!\n";

my $first = -1;
my ($name, $sequence);
foreach my $fasta (@ARGV){
  open my $input, "$fasta" or die "Can't open input fasta $fasta: $!\n";

  
  while(<$input>){
    chomp;
    next if /^\s*$/;
    
    if (/^>/){
      my $newname = $_;
      $newname =~ s/^>//;

      if($first > 0){
        print_gtf($name,$sequence);
      
      }else{
        $first = 1;
      }
    
      $name = $newname;
      $sequence = "";
      next;
    
    }else{
      $sequence .= $_;
    }

  }
  close $input;
}
print_gtf($name,$sequence);

sub print_gtf {
  my $name = $_[0];
  my $sequence = $_[1];

  my @seqinfo = split(' ' , $name, 2);
  my $geneid = $seqinfo[0];
  my $desc = $seqinfo[1] || $seqinfo[0];
  
  print STDOUT join("\t",($geneid,".","gene",1,length($sequence),".","+",".","gene_id \"$geneid\"; gene_name \"$desc\";")) . "\n";
  print STDOUT join("\t",($geneid,".","transcript",1,length($sequence),".","+",".","gene_id \"$geneid\"; transcript_id \"$geneid\_t\"; gene_name \"$desc\";")) . "\n";
  print STDOUT join("\t",($geneid,".","exon",1,length($sequence),".","+",".","gene_id \"$geneid\"; transcript_id \"$geneid\_t\"; gene_name \"$desc\"; exon_number 1;")) . "\n";

}

