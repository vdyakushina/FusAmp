#!/usr/bin/env python3

from subprocess import run, PIPE
import re, os, re, argparse, pathlib

parser = argparse.ArgumentParser(description="Analyse fusion panel")
parser.add_argument("--fasta", help="path to ref fasta file. e.g. ./WG_IAD154041_THCAF2/WG_IAD154041_THCAF_v3.20181217._Reference.fasta")
parser.add_argument("--fastq1", help="path to fastq1 file")
parser.add_argument("--fastq2", help="path to fastq2 file")
parser.add_argument("--bam", help="path to ref fasta file")
parser.add_argument("--adapters", help="path to file with two adapter sequences, each sequnce in new line. e.g. adapters_Illumina.fa")

args = parser.parse_args()
(ref, fastq1, fastq2, bam, adapters)=(os.path.realpath(vars(args)[i]) for i in ['fasta', 'fastq1', 'fastq2', 'bam', 'adapters'])
print(f'fastq1: {fastq1}\nfastq2: {fastq2}')
outdir=os.path.dirname(bam)+'/'
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 

print('Trimm adapters')
ptrn=re.findall(r'.f[ast]*q.gz', fastq1)[0]
adapters=open(adapters, 'r').read().strip().split('\n')

fastq1_trimm=outdir+os.path.basename(fastq1).strip(ptrn)+'_trimmed.fq.gz'
fastq2_trimm=outdir+os.path.basename(fastq2).strip(ptrn)+'_trimmed.fq.gz'
run_cutadapt=f"cutadapt -a {adapters[0]} -A {adapters[1]} -o {fastq1_trimm} -p {fastq2_trimm} {fastq1} {fastq2} -m 50:50 -z"
run(run_cutadapt, shell=True)

print('Align')
run_bwa=f"bwa mem -t 20 {ref} {fastq1_trimm} {fastq2_trimm} | samtools view -bS - | samtools sort - | samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' - -o {bam} -@ 20"
run(run_bwa, shell=True)
run('samtools index %s'%bam, shell=True)
run('rm %s %s' %(fastq1_trimm, fastq2_trimm), shell=True)