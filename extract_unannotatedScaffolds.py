#!/usr/bin/env python

#date June 2020
# @Hasiba Asma
import os
import sys
import argparse
from Bio import SeqIO

#This script is used to extract those scaffolds from the Masked fasta file that does not have any annotated gene in its gff, as we know that from preflight output
#example command
#./extract_unannotateedScaffold.py -po preflightOutput -f nameOfMaskedFasta

def main():
	#global d1
	parser=argparse.ArgumentParser()
	parser.add_argument('-po','--preflightOutput',help='preflight output',required=True)
	parser.add_argument('-f','--fasta',help='fasta file',required=True)
	
	args = parser.parse_args()
	preflightOutPath=os.path.abspath(args.preflightOutput)
	fastaFile=os.path.abspath(args.fasta)
	fastaFileName=args.fasta
	annotatedScaffolds=[]
	started=''
	with open(preflightOutPath) as pf:
		for line in pf:
			if line.startswith('SEQUENCE-REGIONS'):
				#pf.next()
				started='True'
				
			elif started=='True' and line != '\n':
				chrs=line.split('\t')[0]
				annotatedScaffolds.append(chrs)
			elif line=='\n':
				started='False'

	print('number of chromosomes that are annotated',len(annotatedScaffolds))
	#print('Here is the list',annotatedScaffolds)
	fasta_sequence = SeqIO.parse(open(fastaFile),'fasta')
	#open fasta file and exclude those that are unannotated
	newFastaFileName='mappedOnly_' + fastaFileName
	with open(newFastaFileName,'w') as f:
		for seq in fasta_sequence:
			name = seq.id
			nuc = seq.seq
			if name in annotatedScaffolds:
				print(name)
				SeqIO.write([seq], f, "fasta")

main()
