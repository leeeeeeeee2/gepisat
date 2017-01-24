#!/usr/bin/python
#
# file_handler-any.py
#
# written by Tyler W. Davis
# last updated: 2017-01-23
#
# ------------
# description:
# ------------
# This script is meant to speed up file handling for GePiSaT output files.
# 1. Creates subdirectories based on station names if they do not already exist
# 2. Moves files into their corresponding subdirectory
#
##############################################################################
# IMPORT MODULES
##############################################################################
import errno
import glob
import logging
import os
import re
import shutil
import sys


##############################################################################
# FUNCTIONS
##############################################################################
def mkdir_p(path):
    """
    Name:     mkdir_p
    Inputs:   str, directory path (path)
    Outputs:  None.
    Features: Makes directories, including intermediate directories as
              required (i.e., directory tree)
    Ref:      tzot (2009) "mkdir -p functionality in python,"
              StackOverflow, Online:
              http://stackoverflow.com/questions/600268/
                mkdir-p-functionality-in-python
    """
    try:
        logging.debug('building directory tree')
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            logging.debug('directory exists')
            pass
        else:
            logging.error('failed to create directory!')
            raise

##############################################################################
# MAIN
##############################################################################
if __name__ == '__main__':
    usage_str = "python file_handler.py [directory]"

    # Check for command line argument:
    if len(sys.argv) > 1:
        my_dir = sys.argv[1]
        if not os.path.isdir(my_dir):
            raise ValueError(
                "Directory '%s' does not exist! Script usage: %s" % (
                    my_dir, usage_str))
    else:
        my_dir = "."

    # Search for plain text files (change if necessary):
    s_pattern = "*.txt"
    s_str = os.path.join(my_dir, s_pattern)
    my_paths = glob.glob(s_str)
    n_files = len(my_paths)
    if n_files == 0:
        print("No files in '%s'" % (my_dir))
    else:
        for my_path in my_paths:
            my_file = os.path.basename(my_path)

            # Search for station name (could go awry!):
            result = re.search('^(\S{6})_', my_file)
            if result:
                my_station = result.group(1)
                out_dir = os.path.join(my_dir, my_station.lower())
                to_continue = True

                # Create the new subdirectory:
                if not os.path.isdir(out_dir):
                    try:
                        mkdir_p(out_dir)
                    except:
                        print("Failed to create directory '%s'" % (out_dir))
                        to_continue = False
                    else:
                        to_continue = True

                if to_continue:
                    out_path = os.path.join(out_dir, my_file)
                    try:
                        # Move file into subdirectory:
                        shutil.move(my_path, out_path)
                    except:
                        print("Failed to move file '%s' to '%s'" % (
                            my_path, out_dir))
