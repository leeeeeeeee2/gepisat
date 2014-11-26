#!/usr/bin/perl
#
# file_handler.pl
#
# written by Tyler Davis
# 
# 2013-10-28 -- created
# 2014-01-29 -- last updated
#
# ------------
# description:
# ------------
# This script is meant to speed up file handling in windows.
# 1. Creates subdirectories if they do not already exist
# 2. Moves files into their corresponding subdirectory
#
# ----------
# changelog:
# ----------
# 01. added grep to readdir to read only text files [14.01.29]
#
##############################################################################
## MODULES
##############################################################################
use warnings;
use strict;
use File::Path;
use File::Copy;

##############################################################################
## DEFINITIONS
##############################################################################
#my $obs_dir = "C:\\Users\\Tyler\\Desktop\\ro\\";
my $obs_dir = "C:\\Users\\Tyler\\Projects\\gepisat\\data\\fluxdata\\synth_hourly_allvars_csv\\";
my @files;
my $file;
my $stationid;
my $my_subdir;
my $obs_dir_new;

##############################################################################
## MAIN
##############################################################################
opendir DIR, $obs_dir or die "Cannot open dir $obs_dir: $!";
@files = grep {/.*txt/} readdir DIR;
closedir DIR;

foreach $file (@files) {
	# Pull the directory name from file name:
	$stationid = substr($file, 0, 6);
	$my_subdir = lc($stationid);
	#
	# Append subdir to directory path:
	$obs_dir_new = $obs_dir . $my_subdir;
	#
	# Check to see that directory exists:
	if (! -d $obs_dir_new) {
		my $dirs = eval { mkpath($obs_dir_new) };
		die "Failed to create $obs_dir_new: $@\n" unless $dirs;
	}
	#
	# Move the file to its subdirectory:
	my $from_file = $obs_dir . $file;
	$obs_dir_new = $obs_dir_new . "\\";
	move($from_file, $obs_dir_new);
}
