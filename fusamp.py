#!/usr/bin/env python3
from subprocess import run, PIPE
import os, sys, argparse

parser = argparse.ArgumentParser(description="Run pipeline for targeted methylseq (alignment and call)")
parser.add_argument("--fasta", help="path to ref fasta file")
parser.add_argument("--fastq1", help="path to fastq1 file")
parser.add_argument("--fastq2", help="path to fastq2 file")
parser.add_argument("--bam", help="path to ref fasta file")
parser.add_argument("--adapters", help="path to file with two adapter sequences")
parser.add_argument("-p", "--panel", help="path to bed file")
parser.add_argument("--controls", help="path file with control genes")
parser.add_argument("-o", "--output", help="name of output directory to put results")


args = parser.parse_args()
try:
	(fasta, fastq1, fastq2, bam, adapters, panel, controls, output)=(os.path.realpath(i) for i in vars(args).values())
except TypeError:
	parser.print_help()
	sys.exit(1)


scripts=os.path.dirname(__file__)+'/scripts/'

## 1. Align
run('python3 %sfusions_preprocess.py --fasta %s --fastq1 %s --fastq2 %s --bam %s --adapters %s' % (scripts, fasta, fastq1, fastq2, bam, adapters), shell=True, executable='/bin/bash')
## 2. Call C>T and calculate mean metrics
run('python3 %sfusions_assess.py  -p %s --controls %s --bam %s -o %s' % (scripts, panel, controls, bam, output), shell=True, executable='/bin/bash')
