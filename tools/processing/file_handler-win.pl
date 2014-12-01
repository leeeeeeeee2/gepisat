#!/usr/bin/perl
#
# file_handler-win.pl
#
# written by Tyler W. Davis
# 
# 2013-10-28 -- created
# 2014-12-01 -- last updated
#
# ------------
# description:
# ------------
# This script is meant to speed up file handling in Windows.
# 1. Creates subdirectories if they do not already exist
# 2. Moves files into their corresponding subdirectory
#
# ----------
# changelog:
# ----------
# 01. added grep to readdir to read only text files [14.01.29]
# 02. general housekeepling [14.12.01]
#
##############################################################################
## IMPORT MODULES
##############################################################################
use warnings;
use strict;
use File::Path;
use File::Copy;

##############################################################################
## DEFINITIONS
##############################################################################

### USER DEFINED VARIABLE ###
my $obs_dir = "C:\\Users\\Tyler\\Desktop\\ro\\";

my @files;
my $file;
my $stationid;
my $my_subdir;
my $obs_dir_new;

##############################################################################
## MAIN PROGRAM
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
