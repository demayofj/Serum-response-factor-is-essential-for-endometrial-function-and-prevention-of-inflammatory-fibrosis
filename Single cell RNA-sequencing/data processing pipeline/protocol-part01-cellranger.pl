#!/usr/bin/perl

$tool = "/ddn/gs1/home/grimmsa/tools/cellranger-6.0.0/cellranger";
$anno = "/ddn/gs1/home/papasbn/genome_index/10x_mm10_gencode.vM17.filtered";

open(IN, "id_sample_path.key.txt");
while (<IN>) {
  chomp $_; ($id, $s, $path) = split/\t/, $_;
  print "$tool count --include-introns $id $s $path --localcores=20 --localmem=120 --transcriptome=$anno\n";
}
close(IN);


print "\n\nRun from CellRanger folder...\n"; 
# swarm --file RunCR.txt --threads-per-process 20 --gb-per-process 128 --partition gpu
