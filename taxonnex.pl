#! /usr/bin/perl

##############################################################
# This program takes files with sequences for different genes#
# and outputs the number of genes sampled for each taxa,     #
# organized by taxon rank as defined by the NCBI server.     #
#              S. Hedtke, September 2011		     #
# 	            rev. Oct 2011 SMH		    	     #
##############################################################

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


## To run this program, type : perl taxon.pl filename
## where "filename" is a file that contains a list of the names of your alignment files.
## For example:
#	 coi.nex
#        nak.nex
## The alignment files can be in nexus format or simple tab-format (species	sequence)
## The output will be a list of the taxon names organized by classification,
## with counts of the numbers of sequences sampled for each rank.
## For example, if 40 individuals in the taxon Apoidea were sequenced for coi and 3 for nak,
## And ten (of those 40) coi sequences were in the family Apidae :
#	Taxon	coi	nak
#	Apoidea	40	3
#	 Apidae 10	1
## Note that the classification will be according to the Entrez taxon database, which is imperfect!

# You can change the name of your output file here.
my $outfile='GenesSampledByTaxon.xls';

# You shouldn't need to change anything below here unless you want to.
# Note that while the variable names refer to bees (my study group),
# it should work for any group of organisms.

use Bio::Taxon;
use Bio::DB::GenBank;
use Bio::Tree::Tree;

my $filename=$ARGV[0];
my %superfamily;
my %family;
my %tribe;
my %subfamily;
my %genus;
my %seen;
my %counts;

my $errorfile=$outfile.".err";
$errorfile=~s/(\.xls|\.txt|\.xml)//;
open(ERROR,">$errorfile");

## Load the gene file names.
unless (open(NAMES,$filename)) { print "\n   I can't find the file $filename. Closing!\n"; exit;}

while (<NAMES>) { push (@names,$_); }

## Get the species information from the gene file and query the NCBI database for classification.
for (my $g=0;$g<scalar@names;$g++) {
	chomp $names[$g];
	my $filename=$names[$g];
	$names[$g]=~s/(\.nex|\.txt|\.aln)//g;
	my $gene=$names[$g];
print "\n\nGetting info from $gene alignment: ";
	unless (open(GENE,$filename)) { print "\n  I can't find the file $filename. Skipping.\n"; print ERROR "Skipped $filename.\n"; next;}
	while (my $t=<GENE>) {
	if ($t=~/\#/ || $t=~/\;/ || $t=~/matrix/i || $t=~/begin/i) {next;}
	my @info=split(/\_/,$t);
	if (scalar@info<2 || $info[1] eq '' || $info[0] eq '') {next;}
	my $species=$info[0];
print "...";
	if ($species eq '') { next; }
	my $db=Bio::DB::Taxonomy->new(-source => 'entrez');
	my $taxon = $db->get_taxon(-name => $species);

if ($taxon eq '') { print "*** $species wasn't found in Entrez taxonomy database.***"; print ERROR "$species in $gene alignment not found.\n"; next;}
	my $tree = Bio::Tree::Tree->new(-node => $taxon);
if ($tree eq '') { print "*** $species wouldn't load the stupid tree. ??? ***";next;}
	my $beesuper=$tree->find_node(-rank=>'superfamily');
		my $beesupername;
		if ($beesuper eq '') {$beesupername='UnknownSuper'.$species;}
		else {$beesupername=$beesuper->node_name;}
	my $beefam = $tree->find_node(-rank=>'family');
		my $beefamname;
		if ($beefam eq '') {$beefamname='UnknownFam'.$species;}
		else {$beefamname=$beefam->node_name;}
	my $beesub = $tree->find_node(-rank => 'subfamily');
		my $beesubname;
		if ($beesub eq '') {$beesubname='UnknownSubfam'.$species;}
		else {$beesubname=$beesub->node_name;}
	my $beetribe=$tree->find_node(-rank=>'tribe');
		my $beetribename;
		if ($beetribe eq '') {$beetribename='UnknownTribe'.$species;}
		else {$beetribename=$beetribe->node_name;}
	my $beegenus=$tree->find_node(-rank=>'genus');
		my $beegenusname;
		if ($beegenus eq '') {$beegenusname='UnknownGenus'.$species;}
		else { $beegenusname=$beegenus->node_name; }
## Count the number of sequences by taxonomic rank.
	$counts{$gene.$species}++;
	$counts{$gene.$beegenusname}++;
	$counts{$gene.$beefamname}++;
	$counts{$gene.$beesupername}++;
	$counts{$gene.$beetribename}++;
	$counts{$gene.$beesubname}++;
## Save the names of any new taxa.
	next if $seen{$species}++;
	push @{$genus{$beegenusname}},$species;
	next if $seen{$beegenusname}++;
	push @{$tribe{$beetribename}},$beegenusname;
	next if $seen{$beetribename}++;
	push @{$subfamily{$beesubname}},$beetribename;
	next if $seen{$beesubname}++;
	push @{$family{$beefamname}},$beesubname;
	next if $seen{$beefamname}++;
	push @{$superfamily{$beesupername}},$beefamname;

}}

## Print the information to file.
print "Opening $outfile...";
open (OUT,">$outfile");

print OUT "Taxon";
for (my $x=0;$x<scalar@names;$x++) {
	print OUT "\t$names[$x]";}

while (my ($key,@spit)= each (%superfamily)) {
print "...Printing SUPERFAMILY $key";
	print OUT "\n$key";
	for (my $x=0;$x<scalar@names;$x++) {
		print OUT "\t$counts{$names[$x].$key}";}
for (my $i=0;$i<scalar@{$superfamily{$key}};$i++) {
	my $beefamily=${$superfamily{$key}}[$i];

	unless ($beefamily=~/^Unk/) { print OUT "\n $beefamily";
		for (my $x=0;$x<scalar@names;$x++) {
			print OUT "\t$counts{$names[$x].$beefamily}";}}
	for (my $j=0; $j<scalar@{$family{$beefamily}};$j++) {
		my $subfam=${$family{$beefamily}}[$j];

		unless ($subfam=~/^Unk/) {print OUT "\n  $subfam";
		for (my $x=0;$x<scalar@names;$x++) {
			print OUT "\t$counts{$names[$x].$subfam}";}}
	for (my $m=0; $m<scalar@{$subfamily{$subfam}};$m++) {
		my $beetribe=${$subfamily{$subfam}}[$m];

		unless ($beetribe=~/^Unk/) {print OUT "\n  $beetribe";
			for (my $x=0;$x<scalar@names;$x++) {
			print OUT "\t$counts{$names[$x].$beetribe}";}}
	for (my $k=0;$k<scalar@{$tribe{$beetribe}};$k++) {
		my $beegenus=${$tribe{$beetribe}}[$k];
			print OUT "\n   $beegenus";
			for (my $x=0;$x<scalar@names;$x++) {
			print OUT "\t$counts{$names[$x].$beegenus}";}
	for (my $l=0;$l<scalar@{$genus{$beegenus}};$l++) {
		my $beesp=${$genus{$beegenus}}[$l];
			print OUT "\n    $beesp";
			for (my $x=0;$x<scalar@names;$x++) {
			print OUT "\t$counts{$names[$x].$beesp}";}
}}}}}}
print "\nDone.\n";
close ERROR;
close OUT;
exit;
