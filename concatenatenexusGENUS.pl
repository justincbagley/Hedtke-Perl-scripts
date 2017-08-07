#! /usr/bin/perl

#########################################
# This program takes all nexus files    #
# in a directory and concatenates them. #
# One sequence is randomly selected for #
# each genus.							#
# Any genus sampled for only one gene   #
# is removed from the concatenated		#
# matrix.								#
# S. Hedtke Sept 2011					#
# v.2.0 Feb 2012 SMH					#
#########################################

# https://sites.google.com/site/shannonhedtke/Scripts
# This script was developed for Hedtke SM, Patiny S, Danforth BN. The Bee Tree of Life: a supermatrix approach to apoid phylogeny and biogeography.
#  Copyright (C) 2011  Shannon M Hedtke
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#   See <http://www.gnu.org/licenses/>.
# Portions of this code were modified from concatenatenexusfiles.pl version 1.01 , copyright Brian O'Meara, 2007
#	http://www.brianomeara.info ; licensed under the GNU Public License Version 2

## To use this script, place it in a folder that contains the nexus files you want to concatenate.
## Type : perl concatenatenexusGENUS.pl outputfileprefix
## where outputfileprefix is what you want to call your files.
## This program will build concatenated files for use in phylogenetic analysis:
## 	a nexus file (outputfileprefix.nex) with data and charset blocks (e.g., GARLI, PAUP)
## 	a phylip file (.phy) and a partition file (.part) (e.g., RAxML)
## Genera that only have one gene sampled are removed. Only genera with more than one gene
## sequenced will be in the concatenated files.

# You shouldn't need to change anything below, unless you want to.


my $outfileprefix=$ARGV[0];
my %species;
my %seqs;
my %concatseqs;
my @genenames;
my $ntax=0;
my $nchar=0;
my %seen;
my %end;
my %start;
my %codonset;
my %codonpart;
my %seenit;
my $lastend=0;
my $lastspec='';
my $lastgene='';
my $seq='';
my $spec='';
my $sel=0;
my %count;
my %lengthseqs;
my @genusnames;
my %numgenes;

my $nexusfile=$directory.$outfileprefix.'.nex';
my $phylipfile=$directory.$outfileprefix.'.phy';
my $partfile=$directory.$outfileprefix.'.part';

my @filearray=`ls *.nex`;

## This saves the gene names.
foreach $fname (@filearray) {
	$fname=~m/(.*)\.nex/i;
	$gene=$1;
print "\nGene: $gene";
	push @genenames,$gene;
	$start{$gene}=$lastend+1;
	$end{$gene}=$lastend;
	$lastend=0;
	open (FILE,$fname);
	$spec='';
	my $x=0;

## This gets the sequence information from the alignment files.
while (my $t=<FILE>) {
	$seq='';
	if ($t=~/Nexus/i || $t=~/Matrix/i || $t=~/\;/ || $t=~/\,/ || $t=~/format/i || $t eq /\s/) { next; }
	chomp $t;
# In file, each seq is: "Genus sp. subsp. id_GENBANKID    sequence"
# This first split separates the animal information from the genbank id
	my @poo=split (/\_/,$t);
# Then we need to split the animal information and just get the genus
	my @splitpoo=split(/\s/,$poo[0]);
	$genus=$splitpoo[0];
	$count{$genus.$gene}++;
# Then we need to split the genbank id & sequence 
	my @crap=split (/\s/,$poo[1]);
# and what we want is the sequence, which is after all the spaces.
	my $q=scalar@crap;
	$seq=$crap[$q-1];
	$length=length $seq;
	if ($genus ne /\s/ && $seen{$genus}==0) { $ntax++; $seen{$genus}++; push @genusnames,$genus; }
	if ($count{$genus.$gene}<2) {
		$lengthseqs{$genus.$gene}=$length;
		$seqs{$genus.$gene}=$seq;
		}
	elsif ($count{$genus.$gene}>1) {
			$flip=rand(10);
				if ($flip<4.5) { $lengthseqs{$genus.$gene}=$length; $seqs{$genus.$gene}=$seq; }
				if ($flip>4.5) { }
		}
	if ($seenit{$gene}==0 && $seq ne /\s/) { $seenit{$gene}++; $end{$gene}=$end{$gene}+length $seq; $lastend=$end{$gene}; $nchar=$nchar+length $seq; $sel=length $seq;}
}

# The following code was adapted from concatenatenexusfiles.pl written by Brian O'Meara (2007)
	my $rawcodonposset = `grep -i codonposset $fname`;
#This finds the part of the line that contains N:#-# but leaves off \3.
	if ($rawcodonposset=~m/N\:\s*([\d\-\s]+)/) {
		$codonset{$gene}="$codonset{$gene}\n"."charset $gene".".noncod = ";
		$codonpart{$gene}="$codonpart{$gene}\n"."DNA,$gene".".noncod=";
		$rawN=$1;
# This splits the N:#-#, #-# up by just the number parts.
# For every 1,3,5 #, want to follow with a -; every 2,4,6 would follow with , except for last
		@Narray=split(/\D+/,$rawN);
		my $k=0;
		foreach $Nelement (@Narray) {
			$newelement=$Nelement+($nchar-$sel);
			if ($k==scalar@Narray-1) {$codonset{$gene}="$codonset{$gene}"." $newelement"; $codonpart{$gene}="$codonpart{$gene}"." $newelement";}
			elsif ($k==0) { $codonset{$gene}="$codonset{$gene}"." $newelement -"; $codonpart{$gene}="$codonpart{$gene}"." $newelement-"; $k++;}
			else { $codonset{$gene}="$codonset{$gene}"." $newelement,"; $codonpart{$gene}="$codonpart{$gene}"." $newelement,"; $k--;}
}
		$codonset{$gene}="$codonset{$gene}".";";
		$codonpart{$gene}="$codonpart{$gene}".";";
	}
	if ($rawcodonposset=~m/1\:\s*([\d\-\s]+)/) {
		$codonset{$gene}="$codonset{$gene}\n"."charset $gene".".pos1 = ";
		$codonpart{$gene}="$codonpart{$gene}\n"."DNA,$gene".".pos1=";
		$raw1=$1;
		@Narray=split(/\D+/,$raw1);
		my $k=0;
		foreach $Nelement (@Narray) {
			$newelement=$Nelement+($nchar-$sel);
			if ($k==scalar@Narray-1) {$codonset{$gene}="$codonset{$gene}"." $newelement\\3"; $codonpart{$gene}="$codonpart{$gene}"." $newelement\\3";}
			elsif ($k==0) { $codonset{$gene}="$codonset{$gene}"." $newelement -"; $codonpart{$gene}="$codonpart{$gene}"." $newelement-"; $k++;}
			else { $codonset{$gene}="$codonset{$gene}"." $newelement\\3,"; $codonpart{$gene}="$codonpart{$gene}"." $newelement\\3,"; $k--;}
		}
		$codonset{$gene}="$codonset{$gene}".";";
		$codonpart{$gene}="$codonpart{$gene}".";";
	}
	if ($rawcodonposset=~m/2\:\s*([\d\-\s]+)/) {
		$codonset{$gene}="$codonset{$gene}\n"."charset $gene".".pos2 = ";
		$codonpart{$gene}="$codonpart{$gene}\n"."DNA,$gene".".pos2=";
		$raw2=$1;
		@Narray=split(/\D+/,$raw2);
		my $k=0;
		foreach $Nelement (@Narray) {
			$newelement=$Nelement+($nchar-$sel);
			if ($k==scalar@Narray-1) {$codonset{$gene}="$codonset{$gene}"." $newelement\\3"; $codonpart{$gene}="$codonpart{$gene}"." $newelement\\3";}
			elsif ($k==0) { $codonset{$gene}="$codonset{$gene}"." $newelement -"; $codonpart{$gene}="$codonpart{$gene}"." $newelement-"; $k++;}
			else { $codonset{$gene}="$codonset{$gene}"." $newelement\\3,"; $codonpart{$gene}="$codonpart{$gene}"." $newelement\\3,"; $k--;}
		}
		$codonset{$gene}="$codonset{$gene}".";";
		$codonpart{$gene}="$codonpart{$gene}".";";
	}
	if ($rawcodonposset=~m/3\:\s*([\d\-\s]+)/) {
		$codonset{$gene}="$codonset{$gene}\n"."charset $gene".".pos3 = ";
		$codonpart{$gene}="$codonpart{$gene}\n"."DNA,$gene".".pos3=";
		$raw3=$1;
		@Narray=split(/\D+/,$raw3);
		my $k=0;
		foreach $Nelement (@Narray) {
			$newelement=$Nelement+($nchar-$sel);
			if ($k==scalar@Narray-1) {$codonset{$gene}="$codonset{$gene}"." $newelement\\3"; $codonpart{$gene}="$codonpart{$gene}"." $newelement\\3";}
			elsif ($k==0) { $codonset{$gene}="$codonset{$gene}"." $newelement -"; $codonpart{$gene}="$codonpart{$gene}"." $newelement-"; $k++;}
			else { $codonset{$gene}="$codonset{$gene}"." $newelement\\3,"; $codonpart{$gene}="$codonpart{$gene}"." $newelement\\3,"; $k--;}
		}
		$codonset{$gene}="$codonset{$gene}".";";
		$codonpart{$gene}="$codonpart{$gene}".";";
	}
close FILE;
}

## Get ride of genera sampled only once.
for (my $i=0; $i<scalar@genusnames; $i++) {
	for (my $j=0; $j<scalar@genenames;$j++) {
		if ($count{$genusnames[$i].$genenames[$j]} ne '') {$numgenes{$genusnames[$i]}++;}
}
	if ($numgenes{$genusnames[$i]}<2) {
print "\n$genusnames[$i] was removed for having only $numgenes{$genusnames[$i]} genes.";
		$ntax=$ntax-1;
		$seen{$genusnames[$i]}='';
}
}

## Print the alignment.
open (NEXUS,">$nexusfile");
open (PHY, ">$phylipfile");
open (PART, ">$partfile");

print NEXUS "#NEXUS\nBegin data;\n dimensions ntax=$ntax nchar=$nchar;\nformat datatype=dna missing=? gap=-;\nmatrix";
print PHY "$ntax\t$nchar";

while ( my ($sp, $c) = each(%seen) ) {
	if ($seen{$sp} eq '') {next;}
	print NEXUS "\n$sp\t";
	print PHY "\n$sp\t";
	print "\n$sp";
	for (my $i=0; $i<scalar@genenames;$i++) {
		my $id = $sp.$genenames[$i];
		my $leng=$end{$genenames[$i]}-$start{$genenames[$i]};
		if ($seqs{$id} eq '') { for (my $j=0; $j<$leng+1; $j++) {print NEXUS "?"; print PHY "?";}}
		else { print NEXUS $seqs{$id}; print PHY $seqs{$id}; }
}
}
print NEXUS "\n;";
print NEXUS "\nend;";
print NEXUS "\nbegin sets;";
while ( my ($gene, $st) = each(%start) ) {
	my $range=$start{$gene}."-".$end{$gene};
	print NEXUS "\ncharset $gene = $range;";
	print NEXUS "$codonset{$gene}";
	print PART "\nDNA,$gene=$range";
	print PART "$codonpart{$gene}";
}

print NEXUS "\nend;";
close NEXUS;
close PHY;
close PART;
exit;
