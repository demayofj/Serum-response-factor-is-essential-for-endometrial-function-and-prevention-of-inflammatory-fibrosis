#!/usr/bin/perl

$infile = "doublets_scDblFinder/scDblFinder_output.txt";
@samples = split/\n/, `cut -f1 $infile | grep -v sample | uniq`;

open(ALL, ">doublets_scDblFinder/singlets.all.txt");
foreach $id (@samples) {
  open(OUT, ">doublets_scDblFinder/singlets.$id.txt");
  open(IN, "$infile");
  while (<IN>) {
    chomp $_; ($s, $bc, $status, $score) = split/\t/, $_;
    next unless ($s eq $id);
    next unless ($status eq "singlet");
    print OUT "$bc\n";
    print ALL "$bc\n";
  }
  close(IN); close(OUT);
}

