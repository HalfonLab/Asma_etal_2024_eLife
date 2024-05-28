#!/usr/bin/env python

#(c) Hasiba Asma
#December 2022


#This script can be used to find hits per locus and size of locus
#input directory should be the directory which contain all the output of bedtools closest command (make sure there is nothing else in there)

# e.g.,
#./checkSameLocus_ForSimulation.py shuffled_files/

import os
import sys
import argparse
import pprint
import re
import csv




#file=sys.argv[1]
#takes directory instead
directory_in_str=sys.argv[1]

a=os.getcwd()
subdirectory=a+'/'+directory_in_str.strip('/')
#print(subdirectory)

#looping over all the files in the subdirectory
for root, dirs, files in os.walk(subdirectory):
	for filename in files:
		#print(os.path.join(root, filename))
		print(filename)
		#initializing dictionaries and lists

		diction={} #this dictionary would save all locus such that leftFlankedGene-RightFlankedGene would be the key, and value would be all the hits in the locus
		sizeOfLocusDict={} #this dictionary will save size of locus as their values, and locus itself would be the keys
		intronDict={} #for each locus this dictionay would save if the hit is within intron 
		locusList=[] #saving all possible/unique locus as a list

		
		
		with open(os.path.join(root, filename),'r') as infile:
			for line in infile:
				#print(line)
				cols=line.split('\t')
				schr=cols[0] #SCRM chr name
				sstart=cols[1] # SCRM start base pair
				send=cols[2] # SCRM end base pair
				leftF=cols[14] #name of left flanking gene
				rightF=cols[19] #name of right flanking gene
		
				sizeOfSCRM=int(send)-int(sstart)
				#sizeOfLocus=sizeOfSCRM+abs(int(cols[15]))+abs(int(cols[20]))
				sizeOfLocus= (int(cols[17]))-(int(cols[13]))
		
				#exception case where SCRM is overlapping two genes, locus size would be end of right gene - start of left gene
				if (int(sstart) < int(cols[13])) and (int(sstart) < int(cols[17])) and (int(send) > int(cols[13])) and (int(send) > int(cols[17])):
					sizeOfLocus=(int(cols[18]))-(int(cols[12]))
				
				#saving SCRM hit coordinate to coord variable
				coord=schr+':'+sstart+'-'+send
				#left and right flanked gene names- saving them together and calling this locus as 'flankedGenes'
				flankedGenes=leftF+'-'+rightF
				#creating locuslist
				if flankedGenes not in locusList:
					locusList.append(flankedGenes)
					#locusList.append(flankedGenes+':'+str(abs(sizeOfLocus)))
					#creating dictionary to save size of locus, flanked gene as a key
					if flankedGenes not in sizeOfLocusDict:
						sizeOfLocusDict[flankedGenes]=str(abs(sizeOfLocus))
			
					#adding intron information
					if leftF==rightF:
						intronDict[flankedGenes]='intronic'
					else:
						intronDict[flankedGenes]='not-intronic'

				#if locus not in dictionary add it now
				if flankedGenes not in diction:
					#print('flankedgenes not in dictionary yet')
					diction[flankedGenes]=coord
			
				#and if it is already there add it in
				else:
					if coord not in diction[flankedGenes]:
						#print('flanking genes are ',flankedGenes,' butcoord ',coord,' not in dict ',diction[flankedGenes])
						diction[flankedGenes]=diction[flankedGenes]+','+coord
						#print('added now? ',diction[flankedGenes] )



		#pprint.pprint(diction)
		#print(locusList)

		with open('HitsAndSizePerLocus_'+filename,'w') as outfile:
			for item in diction.keys():
				#print(item)
				scrms=diction[item].split(',')
				#print(item,' ',str(len(scrms)),' ',sizeOfLocusDict[item])
				outfile.write(item+'\t'+str(len(scrms))+'\t'+sizeOfLocusDict[item]+'\t'+intronDict[item]+'\t'+filename+'\n')
