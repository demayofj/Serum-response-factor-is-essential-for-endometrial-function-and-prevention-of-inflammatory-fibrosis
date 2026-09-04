#!/usr/bin/perl

### Summarize counts per sample per broad ctype and per subtype.

@listFib = ("FibroblastF2","FibroblastF3","FibroblastUnk","ProliferatingFibroblastSG2","ProliferatingFibroblastM");
@listEpi = ("GlandularEpithelial","LuminalEpithelial","ProliferatingGlandularEpithelial","LowQualEpithelial");
@listImm = ("cDC1","cDC2_mregDC","pDC","Neutrophil","Macrophage","Monocyte","ILC2","gammadeltaT","T_NKT","NK","ProliferatingNK","Plasma","Bcell");
@listBroad = ("Fibroblast","Epithelial","Perivascular","SmoothMuscle","Endothelial","LymphaticEndothelial","Mesothelial","Immune");
@samples = ("Ctrl1","Ctrl2","SrfKO1","SrfKO2");

$infile = "assigned_cell_types.04feb2025.txt";
$outfile = "assigned_cell_types-summary.04feb2025.txt";

%cts = ();  # can keep them all in one hash because there are no shared names between broad & subtype levels...
open(IN, "$infile");
while (<IN>) {
  chomp $_; ($cell, $C, $broad, $sub, $combined) = split/\t/, $_;
  ($s, $bc) = split/\_\_\_/, $cell;
  $cts{$s}{$broad}++;
  $cts{$s}{$sub}++;
}
close(IN);

open(OUT, ">$outfile");
print OUT "#broad cell types"; foreach $s (@samples) { print OUT "\t$s"; } print OUT "\n";
foreach $x (@listBroad) {
  print OUT "$x";
  foreach $s (@samples) { $n = 0+$cts{$s}{$x}; print OUT "\t$n"; }
  print OUT "\n";
}
print OUT "\n";

print OUT "#fibroblast subtypes"; foreach $s (@samples) { print OUT "\t$s"; } print OUT "\n";
foreach $x (@listFib) {
  print	OUT "$x";
  foreach $s (@samples)	{ $n = 0+$cts{$s}{$x}; print OUT "\t$n"; }
  print	OUT "\n";
}
print OUT "\n";

print OUT "#epithelial subtypes"; foreach $s (@samples) { print OUT "\t$s"; } print OUT "\n";
foreach	$x (@listEpi) {
  print OUT "$x";
  foreach $s (@samples) { $n = 0+$cts{$s}{$x}; print OUT "\t$n"; }
  print OUT "\n";
}
print OUT "\n";

print OUT "#immune subtypes"; foreach $s (@samples) { print OUT "\t$s"; } print OUT "\n";
foreach	$x (@listImm) {
  print OUT "$x";
  foreach $s (@samples) { $n = 0+$cts{$s}{$x}; print OUT "\t$n"; }
  print OUT "\n";
}
print OUT "\n";

close(OUT);
