#!/usr/bin/perl

#########################################
# Script for downloading sequences      #
# from the NCBI database given a file   #
# containing GenBank numbers.           #
# S. Hedtke Sept 2011                   #
#########################################

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

# To use this script, type:
#      perl SHretrieve.pl filename
# Your filename should be a tab-delimited text file with a list of species names followed by one or more accession numbers.
# Your accession numbers should be separated from each other and from the species name by tabs.
# For example:
#     Apis mellifera   AY7638429   AY7394827  QQ834728
#     Poopsie doopsie   TT666728   NN82837
# You may need to convert your Excel or Word document into a .txt file if the script doesn't work.
# If you are running perl on a Mac, this is particularly annoying, as the line breaks MUST be saved as UNIX line breaks.
# This program generates the following files:
# 	infofile : contains information about each sequence pulled to a .info file.
#	errorfile : contains information about any sequences that aren't accessed correctly
#	genefiles : each sequence is placed in a separate file by gene name, as that gene is specified in the GenBank record.
#	            This means you may get, for example, gene files "COI" and "cytochrome oxidase I" even though they are
#		    clearly the same gene, since this script relies on the GB annotation.
#		    These files are in a rough fasta format that can be loaded into an alignment program 
#		    (such as MESQUITE http://mesquiteproject.org/mesquite/mesquite.html): 
#			>species1_accession#
#			sequence
#			>species2_accession#
#			sequence
#

# If you want the program to strip the sequences of introns, change this to Y for yes or N for no.
# Note that this only works for genes that have been properly annotated. Some accessions with introns aren't labeled as such.
my $stripintrons="N";

# Here is where information on each gene retrieved is dumped. This can be sorted by a separate script later.
my $directory='';	# if you want your files to be placed in a separate folder, add the path here. e.g., '/home/genes/';
my $infofile=$directory.'infofile.xls';    # you can name your info file here. e.g., 'BeeGeneInfo.xls';
my $errorfile=$directory.'Errors.txt';	   # you can re-name the error file here. e.g., 'BeeGeneErrors.txt';

# You don't need to make any changes below, unless you want to.

use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::AnnotationI;

open(ERROR,">>$errorfile");
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon++;
print ERROR "\n$mon/$mday/$year $hour:$min\n";


# These are variables for the program to use.
my @numbers;

### This loads the accession numbers from the file.
while (my $t=<>) {
  my @poo=split("\t",$t);
  my $lastpoo=scalar@poo-1;
  chomp $poo[$lastpoo];
  for (my $i=1; $i<scalar@poo; $i++) {
    push (@numbers,$poo[$i]);
}}
print "\nInformation loaded.\n";


### Now we retrieve the GenBank files (as objects) from the NCBI database.

for (my $i=0; $i<scalar@numbers; $i++) {
print "\n Grabbing GB $numbers[$i]... ";
    $accession= $numbers[$i];
    $genBank = new Bio::DB::GenBank;   # create GenBank object
    my $seq = $genBank->get_Seq_by_acc($accession);  #  grab the record by the accession number 
    if ($seq eq undef) {print ERROR "\nAccession no. $accession could not be retrieved."; next;}
    my $dna='';
    my $yes=0;
    my $species='';
### Get information from Features : gene, CDS, and species
    for my $features ($seq->get_SeqFeatures) {
      if ($features->has_tag('gene')) {
	for my $val ($features->get_tag_values('gene')) {$gene=$val;}
      } elsif ($features->has_tag('product')) {
	for my $val ($features->get_tag_values('product')) {$gene=$val;}
      } else {$gene='Unknown';}
      if ($features->primary_tag eq 'CDS' && $stripintrons eq 'Y') { $dna=$features->spliced_seq->seq; $yes++; }
      if ($features->has_tag('organism')) {
	for my $val ($features->get_tag_values('organism')) { $species = $val; }
    } elsif ($species eq '') { $species='Unknown'; }}

if ($yes==0) {$dna=$seq->seq(); print ERROR "\n Accession no. $accession, gene $gene, was not stripped of introns.";}

### Get information from Annotation : authors and reference citation

    my $anno_collection = $seq->annotation;
    my @reference=$anno_collection->get_Annotations('reference');
    my $ref=$reference[0];
    my $hash_ref=$ref->hash_tree;
    my $authors=$hash_ref->{'authors'};
    my $journalref=$hash_ref->{'location'};

### Get information about classification of the organism.
    my @classification = $seq->species->classification;

### Print the information to file.

# To the gene file:
my $genefilename = $directory.$gene.".txt";
    unless (-e $genefilename) { open (GENE,">$genefilename"); } else { open (GENE,">>$genefilename");}   # If the file exists, then append the info, otherwise create a new file.
my $id=$species."_".$accession;
print GENE ">$id\n$dna\n";
close GENE;

# To the .info file:
    if (-e $infofile) { open (INFO,">>$infofile"); } else { open (INFO,">$infofile"); }
print INFO "$species\t$gene\t$accession\t$classification[5] / $classification[4] / $classification[3] / $classification[2]\t$authors\t$journalref\n";
close INFO;
  }
close ERROR;

print "\nCOMPLETED!!\n";

exit;
