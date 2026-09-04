#!/usr/bin/perl

### Make table of assigned type and subtype (if relevant) for each cell.

%assiSub = ();

open(IN, "subcluster/epithelial/AssignedSubtypes.epithelial.txt");
while (<IN>) { chomp $_; ($cell, $x) = split/\t/, $_; $assiSub{$cell} = $x; } close(IN);

open(IN, "subcluster/fibroblast/AssignedSubtypes.fibroblast.txt");
while (<IN>) { chomp $_; ($cell, $x) = split/\t/, $_; $assiSub{$cell} =	$x; } close(IN);

open(IN, "subcluster/immune/AssignedSubtypes.immune.txt");
while (<IN>) { chomp $_; ($cell, $x) = split/\t/, $_; $assiSub{$cell} = $x; } close(IN);

%keyBroad = (0 => "Fibroblast", 1 => "Fibroblast", 6 => "Fibroblast", 7 => "Fibroblast", 8 => "Fibroblast", 18 => "Mesothelial",
             2 => "Endothelial", 19 => "LymphaticEndothelial", 4 => "Perivascular", 11 => "SmoothMuscle", 14 => "SmoothMuscle",
             3 => "Epithelial", 9 => "Epithelial",
             5 => "Immune", 13 => "Immune", 17 => "Immune", 16 => "Immune", 12 => "Immune", 10 => "Immune", 15 => "Immune");
open(OUT, ">assigned_cell_types.04feb2025.txt");
open(IN, "broad_clusters/AssignedClusters_res0.4.RPCA-BySample.txt");
while (<IN>) {
  chomp $_; ($cell, $C) = split/\t/, $_;
  $subtype = $assiSub{$cell};
  if ($subtype eq "") { $subtype = "none"; $mergedtype = $keyBroad{$C}; }
  else { $mergedtype = $subtype; }
  print OUT "$cell\t$C\t$keyBroad{$C}\t$subtype\t$mergedtype\n";
}
close(IN); close(OUT);

