#! /usr/bin/perl

#########################################
# This program takes all nexus files    #
# in a directory and concatenates them. #
# Any species sampled for only one gene #
# is removed from the concatenated		#
# matrix.								#
# S. Hedtke Sept 2011					#
# v.2.4 Feb 2012 SMH					#
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

## To use this script, place it in a folder that contains the nexus files you want to concatenate (and no other nexus files).
## Type : perl concatenatenexus.pl outputfileprefix
## where outputfileprefix is what you want to call your files.
## This program will build concatenated files for use in phylogenetic analysis:
## 	a nexus file (outputfileprefix.nex) with data and charset blocks (e.g., GARLI, PAUP)
## 	a phylip file (.phy) and a partition file (.part) (e.g., RAxML)
## The species names are re-formatted to avoid program errors (e.g., spaces, funky punctuation).
## Species that only have one gene sampled are removed. Only species with more than one gene
## sequenced will be in the concatenated files.

# You shouldn't need to change anything below, unless you want to.

my $outfileprefix=$ARGV[0];
my %species;
my %seqs;
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
my @speciesnames;
my %count;

my $nexusfile=$directory.$outfileprefix.'.nex';
my $phylipfile=$directory.$outfileprefix.'.phy';
my $partfile=$directory.$outfileprefix.'.part';

my @filearray=`ls *.nex`;

## Get the gene information from the list of nexus files.
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

## Get the sequence from each nexus file.
while (my $t=<FILE>) {
	$seq='';
	if ($t=~/Nexus/i || $t=~/Matrix/i || $t=~/\;/ || $t=~/\,/ || $t=~/format/i || $t eq /\s/) { next; }
	chomp $t;
	my @poo=split (/\_/,$t);
	$spec=$poo[0];
	$spec=~ s/(\s|\-|\/)/\_/g; $spec=~s/(\(|\)|\:|\'|\")//g;
	my @checksubsp=split(/\_/,$spec);
	if (scalar@checksubsp>1) {$spec=$checksubsp[0].'_'.$checksubsp[1];}
	$count{$spec.$gene}++;
	my @crap=split (/\s/,$poo[1]);
	my $q=scalar@crap;
	$seq=$crap[$q-1];
		if ($spec ne /\s/ && $seen{$spec}==0) { $ntax++; $seen{$spec}++; push @speciesnames,$spec;}
	$seqs{$spec.$gene}=$seqs{$spec.$gene}.$seq;
	if ($seenit{$gene}==0 && $seq ne /\s/) { $seenit{$gene}++; $end{$gene}=$end{$gene}+length $seq; $lastend=$end{$gene}; $nchar=$nchar+length $seq; $sel=length $seq;}
}

# The following code was adapted from concatenatenexusfiles.pl written by Brian O'Meara (2007)
	my $rawcodonposset = `grep -i codonposset $fname`;
# This finds the part of the line that contains N:#-# but leaves off \3.
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

my %numspec;

## This for loop checks to see how many genes are sampled for each species. If a species is undersampled, then it's marked so that it won't be printed to file.
for (my $i=0; $i<scalar@speciesnames; $i++) {
	for (my $j=0; $j<scalar@genenames;$j++) {
		if ($count{$speciesnames[$i].$genenames[$j]} ne '') {$numgenes{$speciesnames[$i]}++; $numspec{$genenames[$j]}++;}
}
	if ($numgenes{$speciesnames[$i]}<2) {
print "\n$speciesnames[$i] was removed for having only $numgenes{$speciesnames[$i]} genes.";
		$ntax=$ntax-1;
		$seen{$speciesnames[$i]}='';
}
	if ($numspec{$genenames[$j]}<3) {
print "\nwarning:$genenames[$j] has less than three species represented -- this isn't ideal for phylogenetic analyses.\n";
}
}

print "\n$ntax remain.\n";

## Print the files.
open (NEXUS,">$nexusfile");
open (PHY, ">$phylipfile");
open (PART, ">$partfile");

print NEXUS "#NEXUS\nBegin data;\n dimensions ntax=$ntax nchar=$nchar;\nformat datatype=dna missing=? gap=-;\nmatrix";
print PHY "$ntax\t$nchar";

while ( my ($sp, $c) = each(%seen) ) {
	if ($seen{$sp} eq '') {next;}
	print NEXUS "\n$sp\t";
	print PHY "\n$sp\t";
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
close PART;
close PHY;
exit;
