#!/usr/bin/perl
#
# ps2pdfcrop.pl
#
# written by Tyler W. Davis
# Catawba College
#
# 2018-10-23 -- created
# 2018-10-23 -- last updated
#
# ------------
# description:
# ------------
# This script converts a multi-page ps image file to a cropped PDF file.
#
#
##############################################################################
## IMPORT MODULES
##############################################################################
use warnings;
use strict;

my $dir = ".";
my @files;
my $file;


##############################################################################
## MAIN PROGRAM
##############################################################################
opendir DIR, $dir or die "Cannot open dir $dir: $!";
@files = grep {/.*ps$/} readdir DIR;
closedir DIR;

foreach $file (@files) {
    print "Converting $file to PDF\n";
    system("ps2pdf", "$file");
}

opendir DIR, $dir or die "Cannot open dir $dir: $!";
@files = grep {/.*pdf$/} readdir DIR;
closedir DIR;

foreach $file (@files) {
    print "Cropping $file\n";
    system("pdfcrop", "$file", '--margins "25 25 25 25"');
}
