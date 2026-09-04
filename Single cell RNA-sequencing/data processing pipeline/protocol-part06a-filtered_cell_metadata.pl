#!/usr/bin/perl

open(IN, "doublets_scDblFinder/singlets.all.txt");
open(OUT, ">initial_QC/cell_ids-filtered.WithMetadata.txt");
while (<IN>) {
  chomp $_; $id = $_;
  @ar = split/\_/, $id;  # get sample name
  $gt = $ar[0]; chop $gt;  # remove replicate number
  $batch = substr($ar[0], -1);  # Ryan wants to check for batch effect (integration, I suspect he means)
  print OUT "$id\t$gt\t$batch\n";
}
close(IN); close(OUT);

