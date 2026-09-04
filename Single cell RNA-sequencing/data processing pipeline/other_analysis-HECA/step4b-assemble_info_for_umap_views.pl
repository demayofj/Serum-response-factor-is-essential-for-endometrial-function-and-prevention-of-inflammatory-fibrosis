#!/usr/bin/perl

%human2color = (); @orderCT_H = ();
open(IN1, "colors_hecaNH.txt"); @in1 = <IN1>; close(IN1);
open(IN2, "celltypes_hecaNH.txt"); @in2 = <IN2>; close(IN2);
for ($i=0; $i<=$#in1; $i++) {
  chomp $in1[$i]; chomp $in2[$i];
  $human2color{$in2[$i]} = $in1[$i];
  push @orderCT_H, $in2[$i];
}

%assiH = (); @cellsH = (); @cellsM = ();
open(IN, "assigned_celltype.hecaNH.txt");
while (<IN>) {
  chomp $_; ($cell, $assi) = split/\t/, $_;
  $assiH{$cell} = $assi;
  push @cellsH, $cell;
}
close(IN);
open(IN, "tranferred_celltype_predictions.hecaNH_to_mouse.txt");
while (<IN>) {
  chomp	$_; ($cell, $assi, $score) = split/\t/, $_;
  next if ($cell eq "");
  $assiH{$cell} = $assi;
  push @cellsM, $cell;
}
close(IN);

%umap_coord = ();
open(IN, "UMAP_coordinates.hecaNH.txt");
while (<IN>) {
  chomp	$_; ($cell, $x,	$y) = split/\t/, $_;
  $umap_coord{$cell} = "$x\t$y";
}
close(IN);
open(IN, "UMAP_coordinates.projected_mouse.txt");
while (<IN>) {
  chomp $_; ($cell, $x, $y) = split/\t/, $_;
  $umap_coord{$cell} = "$x\t$y";
}
close(IN);

@orderCT_M = (); %assiM = (); %mouse2color = ();
open(IN, "../assigned_cell_types.04feb2025.txt");
while (<IN>) {
  chomp $_; @ar = split/\t/, $_;
  if ($ar[4] eq "ProliferatingFibroblastM") { $ar[4] = "ProliferatingFibroblastG2M"; }
  $assiM{$ar[0]} = $ar[4];
}
close(IN);
@orderCT_M = ("FibroblastF2","FibroblastF3","FibroblastUnk","ProliferatingFibroblastSG2","ProliferatingFibroblastG2M","GlandularEpithelial","ProliferatingGlandularEpithelial",
              "LuminalEpithelial","LowQualEpithelial","Perivascular","SmoothMuscle","Endothelial","LymphaticEndothelial","Mesothelial","Macrophage","Monocyte","Neutrophil",
              "cDC1","cDC2_mregDC","pDC","gammadeltaT","ILC2","T_NKT","NK","ProliferatingNK","Plasma","Bcell");
open(IN, "colors_27-for_mouse.txt"); @in = <IN>; close(IN);
for ($i=0; $i<=$#in; $i++) {
  chomp $in[$i];
  $mouse2color{$orderCT_M[$i]} = $in[$i];
}


open(OUT, ">info_for_umap_views.03mar2025.txt");
print OUT "cell\tspecies\tgenotype\tmouseCellType\thumanCellType\tcolorMouse\tcolorHuman\tumapX\tumapY\n";
foreach $cell (@cellsH) {
  print OUT "$cell\thuman\tNA\tnone\t$assiH{$cell}\tnone\t$human2color{$assiH{$cell}}\t$umap_coord{$cell}\n";
}
foreach $cell (@cellsM) {
  if ($cell =~ /Ctrl/) { $gt = "Ctrl"; } elsif ($cell =~ /SrfKO/) { $gt = "SrfKO"; } else { $gt = "fixme"; }
  print OUT "$cell\tmouse\t$gt\t$assiM{$cell}\t$assiH{$cell}\t$mouse2color{$assiM{$cell}}\t$human2color{$assiH{$cell}}\t$umap_coord{$cell}\n";
}
close(OUT);
