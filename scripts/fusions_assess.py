#!/usr/bin/env python3

import pandas as pd
import numpy as np
from subprocess import run, PIPE
import re, os, argparse, pathlib

parser = argparse.ArgumentParser(description="Analyse fusion panel")
parser.add_argument("-p", "--panel", help="path to design bed file. e.g. ./WG_IAD154041_THCAF_v3.20181217._Designed.bed")
parser.add_argument("--controls", help="path file with control genes. Each gene in new line. e.g. ./WG_IAD154041_THCAF2/control_genes.txt")
parser.add_argument("--bam", help="path to bam file")
parser.add_argument("-o", "--output", help="name of output directory to put results")


args = parser.parse_args()
(panel, controls, bam, outdir)=(os.path.realpath(vars(args)[i]) for i in ['panel', 'controls', 'bam', 'output'])

pathlib.Path(outdir).mkdir(parents=True, exist_ok=True) 
control_genes=open(controls, 'r').read().strip().split('\n')

bedtools_str=f'bedtools coverage -b {bam} -a {panel} > {outdir}/coverage.tsv'
bedtools_str_07=f'bedtools coverage -f 0.7 -F 0.7 -b {bam} -a {panel} > {outdir}/coverage_07.tsv'

run(bedtools_str, shell=True, executable='/bin/bash')
run(bedtools_str_07, shell=True, executable='/bin/bash')
total=int(run('samtools view -c %s' % bam, shell=True, executable='/bin/bash', stdout=PIPE).stdout.decode().strip())
total_60=int(run('samtools view -q 60 -c %s' % bam, shell=True, executable='/bin/bash', stdout=PIPE).stdout.decode().strip())
unmapped=int(run('samtools view -f 4 -c %s' % bam, shell=True, executable='/bin/bash', stdout=PIPE).stdout.decode().strip())

coverage=pd.read_csv(f'{outdir}/coverage.tsv', header=None, sep='\t')
coverage_07=pd.read_csv(f'{outdir}/coverage_07.tsv', header=None, sep='\t')

on_target=sum(coverage[6])/total

controls=coverage[coverage[0].str.contains('|'.join(control_genes), regex=True)][[0,6]].set_index(0).to_dict()[6] # slice table to dict with control genes counts
controls_fraction=coverage_07[coverage_07[0].str.contains('|'.join(control_genes), regex=True)][[0,9]].set_index(0).to_dict()[9] # slice table to dict with control genes fraction
controls_07=coverage_07[coverage_07[0].str.contains('|'.join(control_genes), regex=True)][[0,6]].set_index(0).to_dict()[6] # slice table to dict with control genes

controls_sum = sum(list(controls.values()))
controls_07_sum = sum(list(controls_07.values()))
controls_mean_fraction = np.mean(list(controls_fraction.values()))
controls_integrity = controls_07_sum/controls_sum
qc={'total':total, 'unmapped':unmapped, 'total_with_MAPQ60':total_60, 'on_target':on_target, 'controls_07_sum':controls_07_sum,'controls_integrity':controls_integrity, 'controls_mean_fraction': controls_mean_fraction}

result={type:{None:None} for type in ['fusions','imbalance','expression']}

fusions=coverage_07[coverage_07[5].str.contains('TYPE=Fusion')][[0,6]].set_index(0).to_dict()[6]

# Check if the reads support fusion appropriately

def get_match_finish(tbl, read):
	strt,CIGAR=tbl[tbl[0]==read][[3,5]].iloc[0]
	mtchs=[int(i) for j,i in enumerate(re.findall(r'[0-9]+', CIGAR)) if re.findall(r'[A-Z]', CIGAR)[j] not in ['S','H']]
	match_finish=int(strt)+sum(mtchs)
	return(int(strt), match_finish)

def get_dis_val(tbl1, tbl2, read, strt, stp, brk):
	match_finish=max(get_match_finish(tbl1, read)[1], get_match_finish(tbl2, read)[1])
	act_start=min(get_match_finish(tbl1, read)[0], get_match_finish(tbl2, read)[0])
	dis_val=min((brk-act_start)/(brk-strt), (match_finish-brk)/(stp-brk))
	return(dis_val)

for n,v in fusions.items():
	if v>=2:
		frd1,frd2,rvr1,rvr2=[run('samtools view -f%s %s %s' % (flag, bam, n), shell=True, executable='/bin/bash', stdout=PIPE).stdout.decode().strip().split('\n') for flag in ['98', '146', '82', '162']]
		if not all([any(i) for i in [frd1, frd2, rvr1, rvr2]]):
			continue
		else:	
			strt, stp, brk = int(coverage_07[coverage_07[0]==n][1].iloc[0]), int(coverage_07[coverage_07[0]==n][2].iloc[0]), int(coverage_07[coverage_07[0]==n][5].str.split('BREAKPOINT=',expand=True)[1].str.split(';',expand=True)[0].iloc[0])
			brk=strt+brk
			frd1, frd2, rvr1, rvr2=[pd.DataFrame(drn)[0].str.split('\t', expand=True) for drn in [frd1, frd2, rvr1, rvr2]]
			check_results_frd,check_results_rvr=[],[]
			for (check_results, tbl1, tbl2) in zip([check_results_frd, check_results_rvr], [frd1, rvr1],[frd2, rvr2]):
				for read in tbl1[0]:
					if read in list(tbl2[0]):
						dis_val=get_dis_val(tbl1, tbl2, read, strt, stp, brk)
						if dis_val>=0.8:			
							check_results.append(dis_val>=0.8)
					else:
						continue
			if all([check_results_frd, check_results_rvr]):
				result['fusions'][n]=len(check_results_frd+check_results_rvr)
			else:
				continue


five_p=coverage_07[(coverage_07[5].str.contains('TYPE=5p3pAssay'))&(coverage[0].str.contains('5p'))][[0,6]].set_index(0).to_dict()[6]
five_p={i:j/controls_07_sum for (i,j) in five_p.items()}

three_p=coverage_07[(coverage_07[5].str.contains('TYPE=5p3pAssay'))&(coverage[0].str.contains('3p'))][[0,6]].set_index(0).to_dict()[6]
three_p={i:j/controls_07_sum for (i,j) in three_p.items()}

imbalance={name:sum([v for t,v in three_p.items() if name in t])-sum([v for t,v in five_p.items() if name in t]) for name in [n.split('_')[0] for n in three_p.keys()]} # calculate difference 3p - 5p. In case several regions for each end, sum forregions is implimented

expressions=coverage_07[(coverage_07[5].str.contains('TYPE=GeneExpression'))&(~coverage[0].str.contains('|'.join(control_genes), regex=True))][[0,6]].set_index(0).to_dict()[6]
expressions={i:j/controls_07_sum for (i,j) in expressions.items()}

### Write results

for fn,dct in zip(['qc.tsv', 'fusions.tsv', 'imbalance_five_p.tsv', 'imbalance_three_p.tsv', 'imbalance.tsv', 'expression.tsv', 'controls.tsv'], [qc, fusions, five_p, three_p, imbalance, expressions, controls_07]):
	with open(f'{outdir}/{fn}', 'w') as outf:
		outf.write('\n'.join([n+'\t'+str(v) for n,v in dct.items()])+'\n')

final_t=pd.DataFrame(columns=['type', 'name', 'value'])
for key in result.keys():
	if len(result[key])>1:
		result[key].pop(None)
	tmp_t=pd.DataFrame(result[key].items(), columns=['name', 'value'])
	tmp_t['type']=key
	final_t=pd.concat([final_t, tmp_t], axis=0)

final_t.to_csv(f'{outdir}/results.tsv', sep='\t', index=False)
