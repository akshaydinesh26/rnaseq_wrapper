#!/usr/bin/env python

import subprocess
import argparse
import shutil
import os
from Bio import SeqIO
import tempfile

# arguments are transcriptome fasta and fastq file

# validate input file path
def validate_file(file):
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError("{0} does not exist".format(file))
    return file


# Input filepath arguments
parser = argparse.ArgumentParser(description='Transcriptome based RNAseq abundance estimation')
parser.add_argument('-T', dest='transcriptomefile', required=True,
                    help="Reference transcriptome file",type=validate_file)
parser.add_argument('-R1', dest='read1', required=True,
                    help="Fastq reads R1",type=validate_file)
parser.add_argument('-R2', dest='read2', required=True,
                    help="Fastq reads R2",type=validate_file)
args = parser.parse_args()

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
  resultpath = os.path.join(rootpath, "inputData")
  os.mkdir(resultpath)
  inputpath = os.path.join(rootpath, "results")
  os.mkdir(inputpath)


temp_dir = tempfile.TemporaryDirectory()
# temp_dir.cleanup()
# To replace strings in nextflow script
def replacestring(inFilepath,outFilepath,toReplace,replaceWith):
   with open(inFilepath, "rt") as fin:
      with open(outFilepath, "wt") as fout:
         for line in fin:
            fout.write(line.replace(toReplace , replaceWith))


def initialize(): 
   if (checkHeader(args.transcriptomefile)):
      makebin(temp_dir)
      inputpath=os.path.join=(temp_dir,"inputData")
      shutil.copy(args.transcriptomefile,inputpath)
      shutil.copy(args.read1,inputpath)
      shutil.copy(args.read2,inputpath)



  
   








