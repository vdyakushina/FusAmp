# FusAmp

Analize fusions, 3'/5' imbalance, expression markers from targeted amplification based sequencing.

## Results:
1. bam aligned with bwa met<br/>
2. Files with qc and read counts (will be writen to --output directory):<br/>

   2.1. qc.tsv - some qc metrics:
      - total	- total number of raw reads;
      - unmapped - number of reads unmapped;
      - total_with_MAPQ60	number of reads with MAPQ 60;
      - on_target - % of reads aligned on target. Calculation: sum of values received with bedtools coverage against designed bed, devided by total;
      - controls_07_sum	- number of reads aligned on control genes under condition that target covered by read at least to 0.7 and read covered with target at least to 70%. Calculation: sum of values received with bedtools coverage -f 0.7 -F 0.7;
      - controls_integrity - % reads aligned on control genes under condition -f 0.7 -F 0.7 from total number of reads aligned on control genes;
      - controls_mean_fraction - mean value from frctions of targets (accordinf ded file) covered. Calulated from columun containind "% of A at depth" from bedtools coverage resulting file.<br/>

   The main QC metric - controls_07_sum represents how RNA molecules were degradated <br/>
      
   2.2. fusions.tsv - raw read counts from fusions, calculated with bedtools coverage. Contains non confident calls, filtered fusions are in results.tsv.<br/>
   
   2.3. imbalance_five_p.tsv - raw read counts from 5' targets, calculated with bedtools coverage under condition -f 0.7 -F 0.7. Normalised with controls_07_sum.<br/>
   
   2.4. imbalance_three_p.tsv - raw read counts from 3' targets, calculated with bedtools coverage under condition -f 0.7 -F 0.7. Normalised with controls_07_sum.<br/>
   
   2.5. imbalance.tsv Calculated as difference of counts 3' - 5' In case there ara several targets from 5' and 3' ends, the sum over tagets is calculated for each end.<br/>

4. results.tsv
   - fusions- fusions selected from fusions.tsv with algorithm (will describe later)
   - imbalance - imbalance target selected if it is outlier from distribution. Not implemented yet.
   - expression - expression targets selected if it is outlier from distribution. Not implemented yet.
     
## USAGE

### Install

wget https://github.com/vdyakushina/FusAmp/archive/refs/heads/main.zip<br/>
unzip main.zip<br/>
cd bwa-meth-master/<br/>

### Run

  fusamp.py --fasta [reference.fa] --fastq1 [R1.fq.gz] --fastq2 [R2.fq.gz] --bam [bamfile] --adapters [adapters.file] --panel [Designed.bed] --controls [controls.file] --output [output directory]

  help:<br/>

  fusamp.py -h<br/>
  --fasta FASTA         path to ref fasta file<br/>
  --fastq1 FASTQ1       path to fastq1 file<br/>
  --fastq2 FASTQ2       path to fastq2 file<br/>
  --bam BAM             path to ref fasta file<br/>
  --adapters ADAPTERS   path to file with two adapter sequences<br/>
  -p PANEL, --panel PANEL
                        path to bed file<br/>
  --controls CONTROLS   path file with control genes<br/>
  -o OUTPUT, --output OUTPUT
                        name of output directory to put results<br/>

## Prerequisite<br/>
  python3<br/>
  samtools<br/>
  bedtools<br/>
  standart python libraries: pandas, numpy, subprocess, re, os, argparse, pathlib

