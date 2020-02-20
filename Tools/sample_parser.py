#!/usr/bin/env python
import re
import subprocess
import glob
import sys

filepath = sys.argv[1] 
LIST_FILE = glob.glob(filepath + "/*.fastq.gz") # get file names
LIST_NAME = [] # create list to put samples names in later
for name in LIST_FILE:  # Regex magic to get the sample names, only works for illumina files tho, but so does all the pipeline ...
    removeDir = re.sub(filepath + "/", '', name)
    removeDir = re.sub(".fastq.gz", '', removeDir)
    if re.sub("_L\w+_R1_001", '', removeDir) != removeDir:
        name = re.sub("_L\w\w\w_R1_\w\w\w", '', removeDir)
        LIST_NAME.append(name)
LIST_NAME = list(set(LIST_NAME)) # make it all pretty
# write the magic into a txt
with open('samples.txt', 'w') as f:
    for item in LIST_NAME:
        f.write("%s\n" % item)