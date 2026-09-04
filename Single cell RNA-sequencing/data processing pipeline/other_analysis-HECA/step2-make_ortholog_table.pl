#!/usr/bin/perl

open(OUT1, ">genelist_ortholog_info.tmp");
open(OUT2, ">genelist_ortholog_pairs.txt");
open(OUT3, ">genelist_ortholog_mouse2human_translation.txt");  # human symbol if in 1:1 pair, otherwise keep mouse symbol

%mouse2human = (); %human2mouse = (); %orthGrp = ();
open(IN, "/ddn/gs1/home/grimmsa/references/HumanMouseOrthologs/27feb2025/HOM_MouseHuman_ortholog_groups.1-to-1.27feb2025.txt");
while (<IN>) {
  next if ($_ =~ /^\#/);
  chomp $_; ($ogrp, $xM, $listM, $xH, $listH) = split/\t/, $_;  # each "list" should only have one gene since this is 1-to-1 subset!
  ($mouse_sym, $mouse_gid) = split/\|/, $listM;
  ($human_sym, $human_gid) = split/\|/, $listH;
  $orthGrp{$ogrp} = "$mouse_sym\t$human_sym";
  $mouse2human{$mouse_sym} = $human_sym;
  $human2mouse{$human_sym} = $mouse_sym;
}
close(IN);

%ggH = (); open(IN, "genelist.hecaCells-human.txt"); while (<IN>) { chomp $_; $ggH{$_} = 1; } close(IN);
%ggM = (); open(IN, "genelist.ryan-mouse.txt"); while (<IN>) { chomp $_; $ggM{$_} = 1; } close(IN);

%found_pairs = (); %extraH = (); %extraM = ();
foreach $symM (keys %ggM) {
  if (exists $mouse2human{$symM}) {
    if (exists $ggH{$mouse2human{$symM}}) {
       $pair = "$mouse2human{$symM}\t$symM";
       $found_pairs{$pair} = 1;
       delete $ggH{$mouse2human{$symM}};
       print OUT3 "$symM\t$mouse2human{$symM}\n";
    }
    else { $extraM{$symM} = 1; print OUT3 "$symM\t$symM\n"; }
  }
  else { $extraM{$symM} = 1; print OUT3 "$symM\t$symM\n"; }
}
foreach $symH (keys %ggH) { $extraH{$symH} = 1; }

$ctPairs = keys %found_pairs;
$ctSoloH = keys %extraH;
$ctSoloM = keys %extraM;

print "paired = $ctPairs\nextraH = $ctSoloH\nextraM = $ctSoloM\n";
# paired = 13476
# extraH = 4260
# extraM = 13010

foreach $pair (sort (keys %found_pairs)) { print OUT2 "$pair\n"; print OUT1 "paired\t$pair\n"; }
foreach $symM (sort (keys %extraM)) { print OUT1 "onlyMouse\t\t$symM\n"; }
foreach $symH (sort (keys %extraH)) { print OUT1 "onlyHuman\t$symH\t\n"; }
close(OUT1); close(OUT2); close(OUT3);
system "sort -k2,2 -k3,3 genelist_ortholog_info.tmp > genelist_ortholog_info.txt";
system "rm genelist_ortholog_info.tmp";
