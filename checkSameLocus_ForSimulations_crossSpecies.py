#!/usr/bin/env python

#(c) Hasiba Asma
#December 2022


#This script can be used to find hits per locus and size of locus
#input directory should be the directory which contain all the output of bedtools closest command (make sure there is nothing else in there)


#this function takes in ortholog dict file and species name and returns the dictionary compOrthoDict[speciesName][]
def orthologsList(orthoFile,speciesName):
	compOrthoDict ={}
	with open(orthoFile) as orf:
		for line in orf:
			geneIDList=[]
			geneID=''
			line=line.strip('\n')
			cols=line.split('\t')
			orthologID=cols[0]
			if cols[1].find(',')!=-1:
				geneIDList=cols[1].split(',')
				if speciesName not in compOrthoDict:
					compOrthoDict[speciesName] = {}
				if orthologID not in compOrthoDict[speciesName]:
					compOrthoDict[speciesName][orthologID] = geneIDList
			#if there is one gene associated with ortholog
			else:
				geneID=cols[1]
				if speciesName not in compOrthoDict:
					compOrthoDict[speciesName] = {}
				if orthologID not in compOrthoDict[speciesName]:
					compOrthoDict[speciesName][orthologID]=geneID
	return(compOrthoDict)

def dictionaryLocus(spex,locusName,scrmCoord):
	#diction={}
	#if locus not in dictionary add it now
	if spex not in diction:
		diction[spex]={}
	if locusName not in diction[spex]:
		#print('flankedgenes not in dictionary yet')
		diction[spex][locusName]=scrmCoord

	#and if it is already there add it in
	else:
		if scrmCoord not in diction[spex][locusName]:
			#print('flanking genes are ',flankedGenes,' butcoord ',coord,' not in dict ',diction[flankedGenes])
			diction[spex][locusName]=diction[spex][locusName]+','+scrmCoord
								#print('added now? ',diction[flankedGenes] )
	return(diction)
# e.g.,
#./checkSameLocus_ForSimulation.py shuffled_files/

import os
import sys
import argparse
import pprint
import re
import csv
import glob

diction={} 
def main():
	speciesList=[]
	countOfSpecies={}
	
	#file=sys.argv[1]
	#takes directory instead
	directory_in_str=sys.argv[1]
	fbgnidFile=sys.argv[2]
	#ortho file
#	OrthodmelOFile=sys.argv[2]
#	OrthoaaegOFile=sys.argv[3]
#	OrthoagamOFile=sys.argv[4]

	a=os.getcwd()
	subdirectory=a+'/'+directory_in_str.strip('/')
	#print(subdirectory)

	#dmelCompleteOrtholog= orthologsList(OrthoDmelOFile,'Dmel')
	#looping over all the files in the subdirectory
	#pprint.pprint(dmelCompleteOrtholog)
	for root, dirs, files in os.walk(subdirectory):
		for filename in files:
			#print(os.path.join(root, filename))
			#print(filename)
		
		
			#diction={} #this dictionary would save all locus such that leftFlankedGene-RightFlankedGene would be the key, and value would be all the hits in the locus
			sizeOfLocusDict={} #this dictionary will save size of locus as their values, and locus itself would be the keys
			intronDict={} #for each locus this dictionay would save if the hit is within intron 
			locusList=[] #saving all possible/unique locus as a list			
			
			#speciesName=filename.split('_')[2].rstrip('.bed')
			
			string_split = filename.split(".")
			last_part = string_split[-2]
			#last_part_split = last_part.split("_")
			speciesName=last_part	
	
			#print(speciesName)
			if speciesName not in speciesList:
				speciesList.append(speciesName)
			#glob.glob('Ortho'+speciesName+'OFile')
			#dictFile='Ortho'+speciesName+'OFile'
			#print(dictFile)
			completeOrtholog= orthologsList(speciesName+'_ortholog_dictfile.txt',speciesName)
			#read bed file and create output file with HitsAndSizePerLocus_ beginning

	#reading all files one by one
			with open(os.path.join(root, filename),'r') as infile:
				for line in infile:
					#print(line)
					cols=line.split('\t')
					schr=cols[0]
					sstart=cols[1]
					send=cols[2]
					leftF=cols[14]
					rightF=cols[19]
					
					
					if (leftF!= '.') and (rightF!='.'):
						#print(leftF)
						#print(dmelCompleteOrtholog['Dmel'].items())
						for keyO,valO in completeOrtholog[speciesName].items():
							#print(leftF,keyO)
							if leftF == keyO:
								#print(leftF,keyO)
								#print('present')
								#print(dmelCompleteOrtholog[speciesName][leftF])
								leftF=completeOrtholog[speciesName][leftF]
								break

						for keyO,valO in completeOrtholog[speciesName].items():
							#print(leftF,keyO)
							if rightF==keyO:
								#print('present')
								#print(dmelCompleteOrtholog[speciesName][leftF])
								rightF=completeOrtholog[speciesName][rightF]
								break
						#print('none found')
						#try 2:
						#trying to see if there is a fly ortholog then replace the name with that
						
# 			 			for keyO,valO in completeOrtholog[speciesName].items():
# 							#print(leftF,keyO)
# 							if leftF in keyO:
# 								#print('present')
# 								#print(dmelCompleteOrtholog[speciesName][leftF])
# 								leftF=completeOrtholog[speciesName][leftF]
# 								break
# 
# 						for keyO,valO in completeOrtholog[speciesName].items():
# 							#print(leftF,keyO)
# 							if rightF in keyO:
# 								#print('present')
# 								#print(dmelCompleteOrtholog[speciesName][leftF])
# 								rightF=completeOrtholog[speciesName][rightF]
# 								break
						
						
						#if leftF in dmelCompleteOrtholog:
						#	print('yes')
						#	print(OrthoDmelOFile[leftF])
						#else:
						
					
						sizeOfSCRM=int(send)-int(sstart)
						#sizeOfLocus=sizeOfSCRM+abs(int(cols[15]))+abs(int(cols[20]))
						sizeOfLocus= (int(cols[17]))-(int(cols[13]))
		
						#exception case where SCRM is overlapping two genes, locus size would be end of right gene - start of left gene
						if (int(sstart) < int(cols[13])) and (int(sstart) < int(cols[17])) and (int(send) > int(cols[13])) and (int(send) > int(cols[17])):
							sizeOfLocus=(int(cols[18]))-(int(cols[12]))

						coord=schr+':'+sstart+'-'+send
						#left and right flanked gene names- saving them together and calling this locus as 'flankedGenes'
						#print(leftF,rightF)
						flankedGenes=str(leftF)+'-'+str(rightF)
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

						
						#add this to a function.. create dictionary for each species i guess??
						diction=dictionaryLocus(speciesName,flankedGenes,coord)


			#pprint.pprint(diction)
			#exit(0)
			#print(locusList)
			#if there is need to create a file for  size and hits per locus thing
			# with open('HitsAndSizePerLocus_'+filename,'w') as outfile:
# 				for item in diction[speciesName].keys():
# 					#print(item)
# 					scrms=diction[speciesName][item].split(',')
# 					#print(scrms)
# 					#print(item,' ',str(len(scrms)),' ',sizeOfLocusDict[item])
# 					outfile.write(item+'\t'+str(len(scrms))+'\t'+sizeOfLocusDict[item]+'\t'+intronDict[item]+'\t'+filename+'\n')
# 					
	#go through all the dmel genes one by one and scanning all species locus directory, if present add + 1 to count dictionary
	with open(fbgnidFile, 'r') as fb:#, open('finalOutput_orthoPara_test.txt', 'a') as out:
		for line in fb:
			if line.startswith('#') or line == '\n':
				continue

			else:

				geneFBsymb = line.strip('\n')
				geneID = geneFBsymb
				#print('geneid',geneID)
				if geneID not in countOfSpecies:
					countOfSpecies[geneID]=0
				for species in speciesList:
					#print(species)
					for locus in diction[species]:
						out=0
						#print('locus',locus)
						genes=locus.split('-')
						for gene in genes:
							#print('gene',gene)
							#print('geneid',geneID)
							if geneID==gene:
								#print(geneID,gene)
								#print('break')
								countOfSpecies[geneID]+=1
								out=1
								break
						if out==1:
							break
						
								
		#pprint.pprint(countOfSpecies)
		
		#
		#geneListGreaterThan1=[word for word, occurrences in countOfSpecies.items() if occurrences >= 1]
		geneListGreaterThan5=[word for word, occurrences in countOfSpecies.items() if occurrences >= 5]
		geneListGreaterThan10=[word for word, occurrences in countOfSpecies.items() if occurrences >= 10]
		geneListGreaterThan11=[word for word, occurrences in countOfSpecies.items() if occurrences >= 11]
		geneListGreaterThan12=[word for word, occurrences in countOfSpecies.items() if occurrences >= 12]
		geneListGreaterThan13=[word for word, occurrences in countOfSpecies.items() if occurrences >= 13]
		geneListGreaterThan14=[word for word, occurrences in countOfSpecies.items() if occurrences >= 14]
		geneListGreaterThan15=[word for word, occurrences in countOfSpecies.items() if occurrences >= 15]
		geneListGreaterThan16=[word for word, occurrences in countOfSpecies.items() if occurrences >= 16]
		#print(geneListGreaterThan16)
		print(len(geneListGreaterThan5),len(geneListGreaterThan10),len(geneListGreaterThan11),len(geneListGreaterThan12),len(geneListGreaterThan13),len(geneListGreaterThan14),len(geneListGreaterThan15),len(geneListGreaterThan16))
		#['who', 'joey']
				# for k, v in diction.items():
# 					print('v',v)
# 					for x in v.items():
# 				
# 						print('x',x)
# 						print('xkey',x.keys())
		
main()
