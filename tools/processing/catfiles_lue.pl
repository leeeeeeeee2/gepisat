#!/usr/bin/perl
#
# catfiles_lue.pl
#
# written by Tyler W. Davis
# University of Pittsburgh
#
# 2012-01-15 -- created
# 2014-12-01 -- last updated
#
# ------------
# description:
# ------------
# This script will concatenate all the files listed
# in the working directory that end with "*.txt" into 
# a single file.
#
# NOTE: 
# Currently set to concatenate station LUE files together.
#
# ----------
# changelog:
# ----------
# 01. Updated for LUE files [13.11.18]
# 02. Added 'Station,' to headerline [14.10.22]
# 03. General housekeeping [14.12.01]
#
##############################################################################
## IMPORT MODULES
##############################################################################
use strict;
use warnings;

##############################################################################
## DEFINITIONS
##############################################################################
my $outname = "";     # Output file name prefix
my $output = "";      # Output file name (prefix + "_All_Data.txt"
my $ans = "";         # User response to if output file already exists
my $directory = ".";  # File directory for reading / writing
my @files;            # Array for holding files found in directory
my @filelist;         # Sorted files list
my $headerline = "";  # Headerline for output file
my $file = "";        # Individual file names (from filelist)
my $stationid = "";   # Station ID read from file name
my @lines;            # Array of lines read from file
my $line = "";        # Individual line (from lines array)

##############################################################################
## MAIN PROGRAM
##############################################################################
# Prompt for name (used to name the outfile):
print "Enter output name: ";
$outname = <>;
chomp( $outname );

# Name the output file:
$output = $outname . "_All_Data.txt";

# Check to see if that file already exists:
$ans = "y";
if (-e $output) {
    print "That file already exists, continue (y/n)? ";
    $ans = <>;
    chomp( $ans );
}

# Continue if user says it's okay:
if ( lc( $ans ) eq "y" ) 
{
    # Open the directory and read all files matching regex
    # and sort the files you found by their filenames:
    opendir(DIR, $directory) or die $!;
    @files = grep  { /^\w{2}-\w{3}_LUE.txt?/ } readdir(DIR);
    @filelist = sort(@files);
    closedir(DIR);
    
    # Create the output file:
    open(OUT, ">$output") or die $!;
    
    # Open the first file in the working directory:
    open(ALTFILE,"<$filelist[1]") or die $!;
    
    # Print the headerline to the out file:
    $headerline = <ALTFILE>;
    $headerline = 'Station,' . $headerline;
    print OUT $headerline;
    close(ALTFILE);
    close(OUT);
    
    # Check to see if the array has any files:
    if (@files)  
    {
        # Open outfile for writing/appending:
	    open(OUT, ">>$output") or die $!;
        
        foreach $file (@filelist)
        {
            print "$file\n";
            $stationid = substr $file, 0, 6;
            
            # Open data file for reading:
            open(FILE, "<$file") or die $!;
            @lines = <FILE>;
            foreach $line (@lines)
            {
                # Should avoid blank lines.
                if ($line =~ /^\d+.*/)
                {
                    print OUT "$stationid,$line";
                }
            }
            close(FILE);
        }
        close(OUT);
    }
    else
    {
        # Array was empty.
        print "Found no files.\n";
    }
}
else
{
    # User opted to not overwrite existing file.
    print "Quitting.\n";
}
