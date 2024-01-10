#!/usr/bin/env python

import subprocess
from argparse import ArgumentParser
import shutil
import os
from Bio import SeqIO

# arguments are transcriptome fasta and fastq file











# check fasta header is unique
def checkHeader(fasta):
   headers = []
   with open(fasta, "r") as f:
     for record in SeqIO.parse(f, "fasta"):
       headers.append(record.description)
   flag = len(set(headers)) == len(headers)
   return flag



# create directory structure
def makebin(rootpath):
  os.mkdir(rootpath)
  resultpath = os.path.join(rootpath, "results")
  os.mkdir(resultpath)
  inputpath = os.path.join(rootpath, "results")
  os.mkdir(inputpath)


# move files
def move(oldPath,newPath):
  shutil.move(oldPath,newPath)