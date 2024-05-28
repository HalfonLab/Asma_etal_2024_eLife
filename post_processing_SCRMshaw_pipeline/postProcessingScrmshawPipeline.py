#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

####################################################################################-------BRIEF DESCRIPTION------###################################################################################################																																																					#
# Halfon Lab
# Date: Sept 2018																																															#
#Purpose: This script is written to do some post processing on the SCRMshaw multiple offsets output and runs peaks calling algorithm (MACs) in order to get more robust CRM predictions						#
#Input: It needs two input 																																															#
# (i) combined SCRMshaw_offset's output (including multiple training sets and multiple methods)																														#
# (ii) number of predictions to extract from each of the offset's prediction default = 5000																															#
#Outputs of this script are individual BED formated PEAKs file for each of the training set and each method, which can be concatenated to one file and used as an input for the evaluation pipeline script			#
####################################################################################################################################################################################################################
import os
import pybedtools
import statistics
import argparse
import sys
import re
import pprint
import shutil
from scipy import stats
import scipy.stats
import pprint
import csv
import subprocess
from collections import Counter
import numpy as np
import numpy.matlib as npm
import pandas as pd


##############################################------------MODULES/FUNCTION---------------##########################################
# This function will parse combined file of multiple outputs of scrmshaw to individual unique files for each of the training set and statistical method used and will return three methods dictionaries containing all the training sets as their keys
def parse_output(outfile,numparse): 
	# creating three separate dictionaries based on methods used  
	d_hexmcd={}
	d_imm={}
	d_pac={}
	x={''}
	global d_hexmcd_val
	global d_imm_val
	global d_pac_val
	d_hexmcd_val=['']
	d_imm_val=['']
	d_pac_val=['']

	with open (outfile) as file:
		rows=(line2.split('\t') for line2 in file)
		for row in rows:
		#based on the 14th column(names of different data sets) and 15th column (statistical method used) of scrmshaw_joined_output file giving values to each of the three method`s dictionaries
			if (row[15]=='hexmcd') and (int(row[16]) <= int(numparse)):
				#print(row[16])
				#print(numparse)
				if row[14] not in d_hexmcd:
					myRow = [] # create a new list to use
					myRow.append(row) # add my new row to our new list
					d_hexmcd[row[14]] = myRow  #create a new entry if it isn't in the dictionary already
				else:
					d_hexmcd.get(row[14]).append(row)
					#count_hexmcd=count_hexmcd+1
			elif (row[15]=='imm')and (int(row[16]) <= int(numparse)):
				if row[14] not in d_imm:
					myRow = []
					myRow.append(row)
					d_imm[row[14]] = myRow
				else:
					d_imm.get(row[14]).append(row)
			elif (row[15]=='pac') and (int(row[16]) <= int(numparse) ):
				if row[14] not in d_pac:
					myRow = []
					myRow.append(row)
					d_pac[row[14]] = myRow
				else:
					d_pac.get(row[14]).append(row)

			
		#calculating number of keys(datasets) each dictionary ends up having		
		for key in d_hexmcd.keys():
			d_hexmcd_val.append(key)
		for key in d_imm.keys():
			d_imm_val.append(key)
		for key in d_pac.keys():
			d_pac_val.append(key)
			
	#creating separate files for each method w.r.t datasets, using the above newly created three dictionaries and moving them to tmp (temporary folder).
	
	#These individual unique files(based on methods and their training sets) will go through downstream processes(calling peaks through MACs) one by one
	for key in d_hexmcd.keys():
		noOflines=len(d_hexmcd[key])
		
		with open(os.path.join(subdirectory,'hexmcd_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_hexmcd[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
	
	for key in d_imm.keys():
		with open(os.path.join(subdirectory,'imm_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_imm[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")

	for key in d_pac.keys():
		with open(os.path.join(subdirectory,'pac_'+key+'_fullLength.bed'),'w') as outfile:
			for i in d_pac[key]:
				line=str(i)
				line=line.replace('\'','')
				line=line.strip('[')
				line=line.strip(']')
				line=line.strip(':\\n')
				outfile.write(" ".join(line.split()).replace(', ', '\t'))
				outfile.write("\n")
				
	return(d_imm_val,d_hexmcd_val,d_pac_val)


#------------------------------------------------------------------------------------------------------

# This function will extract user specified number of predictions from each of the offset for the given training set
def extract_topN_scrms_scoreCurve2(fullLengthFilePath,cutoff,method,TSET):
	#extract num of scrms from offset
	#extractedFileName=str(cutoff)+'.'+method+'_'+TSET+'.bed'
	extractedFileName2='all_'+str(cutoff)+'.'+method+'_'+TSET+'.bed'
	i=1
	with open(fullLengthFilePath,'r') as infile, open(os.path.join(subdirectory,extractedFileName2),'w') as outfile:
		prev_line=''
		lines=[]
		for current_line in infile:
			#current_line2=current_line.strip()
			#lines.append(current_line2)
			prev_line2=prev_line.strip()
			lines.append(prev_line2)
			col=current_line.split('\t')
			#print(col[10])
			rank=col[16].strip('\n')

			if int(rank) <= int(cutoff):
				outfile.write(current_line)
			###del
			if prev_line!='':
				prev_rank=prev_line.split('\t')[16].strip('\n')
				#print('rank of this line ',rank,' rank of prev line ',prev_rank)
				if int(prev_rank) > int(rank):
					#print("------------New file------------")
					with open(os.path.join(subdirectory,'file'+str(i)+'_'+method+'_'+TSET+'.txt'), 'w') as f:
						for item in lines:
							f.write("%s\n" % item)
					lines=[]
					i+=1
				
			prev_line=current_line
		else:
			with open(os.path.join(subdirectory,'file'+str(i)+'_'+method+'_'+TSET+'.txt'), 'w') as f:
				for item in lines:
					f.write("%s\n" % item)
			###del 
	#print('number of files generated ',i)
	fileX=str(cutoff)+'.'+method+'_'+TSET+'.bed'
	#tmp2= open(os.path.join(subdirectory,'cutoff_concatenated.txt'),'a')
	tmp2= open(os.path.join(subdirectory,fileX),'a')
	#print(tmp2)
	for j in range(1,i+1):
		#print("----------------------File no. ",j)
		valuesScore=[]
		fileName= os.path.join(subdirectory,'file'+str(j)+'_'+method+'_'+TSET+'.txt')
		#print(fileName)
		with open(fileName,'r') as s:
			for line in s:
				if line!= '\n':
					#print(line)
					row=line.split('\t')
					valuesScore.append(float(row[3]))

		#print('list of scrmshae score of this file is: ',valuesScore)
		# pull out the list from pandas frame
		valuesScore=list(valuesScore)
		valuesScore=sorted(valuesScore,key=float,reverse=True)
		#valuesAmp=sorted(valuesAmp,key=float,reverse=True)

		#for scores cutoff

		#get coordinates of all the points
		nPointsScore = len(valuesScore)
		allCoordScore = np.vstack((range(nPointsScore), valuesScore)).T
		#np.array([range(nPoints), values])

		# get the first point
		firstPointScore = allCoordScore[0]
		# get vector between first and last point - this is the line
		lineVecScore = allCoordScore[-1] - allCoordScore[0]
		lineVecNormScore = lineVecScore / np.sqrt(np.sum(lineVecScore**2))

		# find the distance from each point to the line:
		# vector between all points and first point
		vecFromFirstScore = allCoordScore - firstPointScore

		# To calculate the distance to the line, we split vecFromFirst into two 
		# components, one that is parallel to the line and one that is perpendicular 
		# Then, we take the norm of the part that is perpendicular to the line and 
		# get the distance.
		# We find the vector parallel to the line by projecting vecFromFirst onto 
		# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
		# We project vecFromFirst by taking the scalar product of the vector with 
		# the unit vector that points in the direction of the line (this gives us 
		# the length of the projection of vecFromFirst onto the line). If we 
		# multiply the scalar product by the unit vector, we have vecFromFirstParallel
		scalarProductScore = np.sum(vecFromFirstScore * npm.repmat(lineVecNormScore, nPointsScore, 1), axis=1)
		vecFromFirstParallelScore = np.outer(scalarProductScore, lineVecNormScore)
		vecToLineScore = vecFromFirstScore - vecFromFirstParallelScore
		# distance to line is the norm of vecToLine
		distToLineScore = np.sqrt(np.sum(vecToLineScore ** 2, axis=1))
		# knee/elbow is the point with max distance value
		idxOfBestPointScore = np.argmax(distToLineScore)

		#print("value of score at score elbow")
		#print(valuesScore[idxOfBestPointScore])
		scoreAtElbow=valuesScore[idxOfBestPointScore]
				
		#print('------File no. ',j,' score elbow is at ',scoreAtElbow)
		
		########## extracting the predictions above the score elbow
		#count=0
		fileName2= 'cutoff_AllFiles_'+fileName
		#tmp2= open(fileName2,'a')
		with open(fileName,'r') as fileA:
			for lineA in fileA:
				if lineA!= '\n':
					rowA=lineA.split('\t')
					#print(line2)
					#if the peak's score value is equal or above the score cutoff only then it will write to file
					if float(rowA[3]) >= scoreAtElbow:
						#print(row[16])
						tmp2.write(lineA)
						#count+=1
	tmp2.close()
	#path=os.path.abspath(extractedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==fileX:
				path=os.path.abspath(os.path.join(root,name))
	
	return path
	
#------------------------------------------------------------------------------------------------------

# This function will extract user specified number of predictions from each of the offset for the given training set
def extract_topN_scrms_scoreCurve(fullLengthFilePath,cutoff,method,TSET):
	#extract num of scrms from offset
	#extractedFileName=str(cutoff)+'.'+method+'_'+TSET+'.bed'
	#extractedFileName2='all_'+str(cutoff)+'.'+method+'_'+TSET+'.bed'
	#i=1
	# with open(fullLengthFilePath,'r') as infile, open(os.path.join(subdirectory,extractedFileName2),'w') as outfile:
# 		prev_line=''
# 		lines=[]
# 		for current_line in infile:
# 			#current_line2=current_line.strip()
# 			#lines.append(current_line2)
# 			prev_line2=prev_line.strip()
# 			lines.append(prev_line2)
# 			col=current_line.split('\t')
# 			#print(col[10])
# 			rank=col[16].strip('\n')
# 
# 			if int(rank) <= int(cutoff):
# 				outfile.write(current_line)
# 			###del
# 			if prev_line!='':
# 				prev_rank=prev_line.split('\t')[16].strip('\n')
# 				#print('rank of this line ',rank,' rank of prev line ',prev_rank)
# 				if int(prev_rank) > int(rank):
# 					#print("------------New file------------")
# 					with open(os.path.join(subdirectory,'file'+str(i)+'_'+method+'_'+TSET+'.txt'), 'w') as f:
# 						for item in lines:
# 							f.write("%s\n" % item)
# 					lines=[]
# 					i+=1
# 				
# 			prev_line=current_line
# 		else:
# 			with open(os.path.join(subdirectory,'file'+str(i)+'_'+method+'_'+TSET+'.txt'), 'w') as f:
# 				for item in lines:
# 					f.write("%s\n" % item)
# 			###del 
# 	#print('number of files generated ',i)
	#fileX=str(cutoff)+'.'+method+'_'+TSET+'.bed'
#	tmp2= open(os.path.join(subdirectory,'cutoff_concatenated.txt'),'a')
 	#tmp2= open(os.path.join(subdirectory,fileX),'w')
# 	#print(tmp2)
# 	for j in range(1,i+1):
# 		#print("----------------------File no. ",j)
	valuesScore=[]
# 		fileName= os.path.join(subdirectory,'file'+str(j)+'_'+method+'_'+TSET+'.txt')
# 		#print(fileName)
	with open(fullLengthFilePath,'r') as s:
		for line in s:
			if line!= '\n':
				#print(line)
				row=line.split('\t')
				valuesScore.append(float(row[3]))

		#print('list of scrmshae score of this file is: ',valuesScore)
	# pull out the list from pandas frame
	valuesScore=list(valuesScore)
	valuesScore=sorted(valuesScore,key=float,reverse=True)
		#valuesAmp=sorted(valuesAmp,key=float,reverse=True)

	#for scores cutoff
	#get coordinates of all the points
	nPointsScore = len(valuesScore)
	allCoordScore = np.vstack((range(nPointsScore), valuesScore)).T
	#np.array([range(nPoints), values])

	# get the first point
	firstPointScore = allCoordScore[0]
	# get vector between first and last point - this is the line
	lineVecScore = allCoordScore[-1] - allCoordScore[0]
	lineVecNormScore = lineVecScore / np.sqrt(np.sum(lineVecScore**2))

		# find the distance from each point to the line:
		# vector between all points and first point
	vecFromFirstScore = allCoordScore - firstPointScore

		# To calculate the distance to the line, we split vecFromFirst into two 
		# components, one that is parallel to the line and one that is perpendicular 
		# Then, we take the norm of the part that is perpendicular to the line and 
		# get the distance.
		# We find the vector parallel to the line by projecting vecFromFirst onto 
		# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
		# We project vecFromFirst by taking the scalar product of the vector with 
		# the unit vector that points in the direction of the line (this gives us 
		# the length of the projection of vecFromFirst onto the line). If we 
		# multiply the scalar product by the unit vector, we have vecFromFirstParallel
	scalarProductScore = np.sum(vecFromFirstScore * npm.repmat(lineVecNormScore, nPointsScore, 1), axis=1)
	vecFromFirstParallelScore = np.outer(scalarProductScore, lineVecNormScore)
	vecToLineScore = vecFromFirstScore - vecFromFirstParallelScore
	# distance to line is the norm of vecToLine
	distToLineScore = np.sqrt(np.sum(vecToLineScore ** 2, axis=1))
	# knee/elbow is the point with max distance value
	idxOfBestPointScore = np.argmax(distToLineScore)

		#print("value of score at score elbow")
		#print(valuesScore[idxOfBestPointScore])
	scoreAtElbow=valuesScore[idxOfBestPointScore]
				
	#print('------concatenated file score elbow is at ',scoreAtElbow)
	#sys.exit()
	fileX=str(cutoff)+'.scoreOf_'+str(scoreAtElbow)+'.'+method+'_'+TSET+'.bed'
	tmp2= open(os.path.join(subdirectory,fileX),'w')
		########## extracting the predictions above the score elbow
		#count=0
	#fileName2= 'cutoff_AllFiles_'+fileName
		#tmp2= open(fileName2,'a')
	with open(fullLengthFilePath,'r') as fileA:
		for lineA in fileA:
			if lineA!= '\n':
				rowA=lineA.split('\t')
				#print(line2)
				#if the peak's score value is equal or above the score cutoff only then it will write to file
				if float(rowA[3]) >= scoreAtElbow:
					#print(row[16])
					tmp2.write(lineA)
						#count+=1
	tmp2.close()
	#path=os.path.abspath(extractedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==fileX:
				path=os.path.abspath(os.path.join(root,name))
	
	return path,scoreAtElbow
#---------------------------------------------------------------------------------------------------------------------------
#This function takes the tab delimited file and returns bed version of it to perform the functions of bed

def bed_conversion(tab_delimited_path):
	
	bed_tabdelimited=pybedtools.BedTool(tab_delimited_path)
		
	return bed_tabdelimited

#---------------------------------------------------------------------------------------------------------------------------
#This function will take in tab delimited file path and convert it into py bed version of that file(on which bedtools functions can be applied like sort/merge etc) and return its path
def bedtools_sorting(extractedScrmsPathBED,sortedFileName):
	
	extractedScrmsPathBED.sort().saveas(subdirectory+'/'+sortedFileName)	
	
	#path=os.path.abspath(sortedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==sortedFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	return path
	
	
	
#---------------------------------------------------------------------------------------------------------------------------
#This function will take in BED file and sum up the score for each of the 10 bp overlapping window and return every 10 bp window with its score(i.e summed)
def sum_of_score(sortedFilePath,sortedFileName,sumOfScoreFileName):
	#calculating sum of score and saving it to a dictionary
	diction={}

	with open(sortedFilePath,'r') as infile:

		for line in infile:
			col=line.split('\t')
			chrName=col[0]
			startCoord=col[1]
			endCoord=col[2]
			score=float(col[3])

			j=int(startCoord)
			tempI=int(startCoord)+10
			bp=range(tempI,int(endCoord)+10,10)

			for i in bp:
				keyName=chrName+':'+str(j)+'-'+str(i)
				#print(keyName)
				if keyName not in diction:
					#print('No')
					diction[keyName]=score

				else:
					#print('yes')
					diction[keyName]=diction[keyName]+score	
		
				j+=10	


	file2=sortedFileName.strip('bed')+'csv'

	#converting dictionary to csv file
	outfile=csv.writer(open('sumScored'+file2,'w'))
	for key, val in diction.items():
		outfile.writerow([key,val])			
	#pprint.pprint(diction)	


	#from csv to bed file
	file1='sumScored'+file2

	with open(file1,'r') as infile, open(os.path.join(subdirectory,sumOfScoreFileName),'w') as outfile: 

		for line in infile:
			col=line.split(',')
			coord=col[0]
			score=col[1]
			chrAndCoords=coord.split(':')
			chrName=chrAndCoords[0]
			Coords=chrAndCoords[1]

			bothCoords=Coords.split('-')
			startCoord= bothCoords[0]
			endCoord=bothCoords[1]

			outfile.write(chrName+'\t'+startCoord+'\t'+endCoord+'\t'+score)

	
	#path=os.path.abspath(sumOfScoreFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==sumOfScoreFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function will take in BED file and return whole genome coverage including any of the missing coordinates with the score value of 0.
def filling_missing_coords(sortedSumOfScoreFilePath,wholeGenomeFileName):
	with open(sortedSumOfScoreFilePath,'r') as infile,open(os.path.join(subdirectory,wholeGenomeFileName),'w') as outfile:
		prevEnd=0
		chrNames={}
		for line in infile:
			col=line.split("\t")

			chrName=col[0]
			start=col[1]
			end=col[2]
			score=col[3]
			if chrName not in chrNames:
				chrNames[chrName]=''
				prevEnd=0
			if prevEnd==start:
				outfile.write(line)

			elif int(start) > int(prevEnd) :
				outfile.write(chrName+'\t'+str(prevEnd)+'\t'+start+'\t'+str('0.0')+'\n')
				outfile.write(line)


			prevEnd=end
	
	#path=os.path.abspath(wholeGenomeFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==wholeGenomeFileName:
				path=os.path.abspath(os.path.join(root,name))
	
	
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function will run MACS program through system call.
def callingMACS(wholeGenomeFilePath,macsOutputName,cutoff):
	subprocess.call(["macs2","bdgpeakcall","-i",wholeGenomeFilePath,"-c",str(cutoff),"-o",str(macsOutputName)])
	
	path=os.path.abspath(macsOutputName)
	
	return path
#---------------------------------------------------------------------------------------------------------------------------
#This function is converting the output of MACs(narrowPeak) file to a BED version of it (more like a SCRMshaw output)
def peaksToScrms(macsOutputPath,peaksToScrmsName,TSET):

	with open(macsOutputPath,'r') as infile, open(os.path.join(subdirectory,peaksToScrmsName),'w') as outfile:
		for line in infile:
			if not line.startswith('track'):
				col=line.split('\t')
				if '_' not in TSET:
					setAndmeth=col[3].split('_')
					set=setAndmeth[0]
					meth=setAndmeth[1]
					amplitude=str(int(col[4])/10)
			
			#mapping1.adult_mesoderm_imm_narrowPeak5
				else:
					setAndmeth=col[3].split('_')
					set=setAndmeth[0]+'_'+setAndmeth[1]
					meth=setAndmeth[2]
					amplitude=str(int(col[4])/10)	
					
				outfile.write(col[0]+'\t'+col[1]+'\t'+col[2]+'\t'+amplitude+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+'0'+'\t'+set+'\t'+meth+'\n')
				
	#path=os.path.abspath(peaksToScrmsName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==peaksToScrmsName:
				path=os.path.abspath(os.path.join(root,name))
	return path
	

#---------------------------------------------------------------------------------------------------------------------------
#This function will retrieve the SCRMshaw score of the peaks by intersecting it with SCRMshaw prediction file
def intersect_peaks_and_scrms(peaksToScrmsPathBED,extractedScrmsPathBED,intersectedFileName):
	peaksToScrmsPathBED.intersect(extractedScrmsPathBED,loj=True).merge(c=[4,21,22,23,24,25,26,27,28,29,30,31,32,33],o=['max','max','distinct','distinct','max','distinct','distinct','distinct','distinct','max','distinct','distinct','distinct','distinct']).saveas(subdirectory+'/'+intersectedFileName)
	
#	path=os.path.abspath(intersectedFileName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==intersectedFileName:
				path=os.path.abspath(os.path.join(root,name))

	return path
	
#---------------------------------------------------------------------------------------------------------------------------
#This function will sort the peaks output file according to the amplitude of their peaks and rank it based on that.
def sortAndRank_basedOnAmplitude(intersectedFilePath,finalPeaksFileName):
	data=pd.read_csv(intersectedFilePath,delimiter='\t',header=None)
	numOfPeaks=len(data)+1
	ranks=pd.Series(range(1,numOfPeaks))
	dataSorted=data.sort_values(by=3,ascending=False)
	dataSorted=dataSorted.reset_index(drop=True)
	dataSorted[18]=ranks
	#finalName='peaksFinal_'+TSET+'_'+method
	dataSorted.to_csv(finalPeaksFileName,sep='\t',index=False,header=False)
	
	path=os.path.abspath(finalPeaksFileName)
	return path,numOfPeaks

#---------------------------------------------------------------------------------------------------------------------------
#This function will extract top N peaks based on user's choice of selection from amplitude curve 
def extract_topN_scrms_amplitudeCurve(finalExtractedPeaksFileName,finalPeaksFilePath,topNmethod):

	row=[]
	valuesAmp=[]

	with open(finalPeaksFilePath) as filePeaks:
		for line in filePeaks:
			row=line.split('\t')
			valuesAmp.append(float(row[3]))
	
	#---
	#Finding ELBOW point on Amplitude curve
	valuesAmp=list(valuesAmp)
	#print(valuesAmp)
	#get coordinates of all the points
	nPointsAmp = len(valuesAmp)
	allCoordAmp = np.vstack((range(nPointsAmp), valuesAmp)).T
	#np.array([range(nPoints), values])

	# get the first point
	firstPointAmp = allCoordAmp[0]
	# get vector between first and last point - this is the line
	lineVecAmp = allCoordAmp[-1] - allCoordAmp[0]
	lineVecNormAmp = lineVecAmp / np.sqrt(np.sum(lineVecAmp**2))

	# find the distance from each point to the line:
	# vector between all points and first point
	vecFromFirstAmp = allCoordAmp - firstPointAmp
	scalarProductAmp = np.sum(vecFromFirstAmp * np.matlib.repmat(lineVecNormAmp, nPointsAmp, 1), axis=1)
	vecFromFirstParallelAmp = np.outer(scalarProductAmp, lineVecNormAmp)
	vecToLineAmp = vecFromFirstAmp - vecFromFirstParallelAmp
	# distance to line is the norm of vecToLine
	distToLineAmp = np.sqrt(np.sum(vecToLineAmp ** 2, axis=1))

	# knee/elbow is the point with max distance value
	idxOfBestPointAmp = np.argmax(distToLineAmp)
	#print("value at 843")
	#print(valuesAmp[idxOfBestPointAmp])
	ampAtElbow=valuesAmp[idxOfBestPointAmp]
	#print(idxOfBestPointScore,idxOfBestPointAmp)
	
	print("elbow amplitude",ampAtElbow)
	print("rank of elbow amplitude",idxOfBestPointAmp)
	#---
	#Finding MEDIAN at amplitude curve
	medianAmp=statistics.median(valuesAmp)
	print("median of amplitude curve is ",medianAmp)
	#find the closest rank for the mean of the amplitude curve
	closestValTomedian=min(valuesAmp,key=lambda x:abs(x-medianAmp))
	##print("closest value of amplitude to median of amplitude curve is ",closestValTomedian)
	closestRankTomedian=valuesAmp.index(closestValTomedian)+1
	print("rank of the closest value of amplitude to median of amplitude curve is ",closestRankTomedian)

	#print("sending this value of rank of median outside this function",closestRankTomedian)
	elbowPointRank=idxOfBestPointAmp
	medianPointRank=closestRankTomedian
	##return idxOfBestPointAmp,scoreAtElbow,ampAtElbow
	##return closestRankTomedian,scoreAtElbow,ampAtElbow

	if topNmethod == 'elbow':
		countE=0
		fN=finalExtractedPeaksFileName+'_ElbowPointAmplitudeCurve_'+str(elbowPointRank)+'_peaks.bed'
		with open(finalPeaksFilePath) as fileRead, open(fN,'w') as fileWrite:
			for line in fileRead:
				row=line.split('\t')
				lastCol=len(row)-1
				if int(row[lastCol]) <= elbowPointRank:
					fileWrite.write(line)
					countE+=1
		numberOfFinalPeaks=countE
	elif topNmethod == 'median':
		countM=0
		fN=finalExtractedPeaksFileName+'_MedianPointAmplitudeCurve_'+str(medianPointRank)+'_peaks.bed'
		with open(finalPeaksFilePath) as fileRead, open(fN,'w') as fileWrite:
			for line in fileRead:
				row=line.split('\t')
				lastCol=len(row)-1
				if int(row[lastCol]) <= medianPointRank:
					fileWrite.write(line)
					countM+=1
		numberOfFinalPeaks=countM
	pathF=os.path.abspath(fN)
	return pathF,numberOfFinalPeaks,fN

#---------------------------------------------------------------------------------------------------------------------------
#this funnction will take out gene id from the line of gff;  could vary a little bit
# 
# def extract_gene_id(input_string):
#     # Define a regular expression pattern to match 'ID=gene-' followed by any characters until a semicolon or space
#     pattern = r'ID=gene-[^; ]+'
# 
#     # Use the findall function to extract all matching substrings
#     gene_ids = re.findall(pattern, input_string)
# 
#     # If gene_ids is not empty, extract the first match
#     if gene_ids:
#         return gene_ids[0].replace('ID=', '')
#     else:
#         return None  # Return None if no match is found


def extract_gene_id(input_string):
    # Try to match 'ID=' followed by any characters until a semicolon or space
    id_match = re.search(r'ID=([^; ]+)', input_string)
    
    if id_match:
        gene_id = id_match.group(1)
        #if not gene_id.startswith('gene-'):
         #   gene_id = 'gene-' + gene_id
        return gene_id
    else:
        return None  # Return None if no 'ID=' is found
#---------------------------------------------------------------------------------------------------------------------------
#This function will read in gff file and save all genes information (coordinates) into a dictionary
def parse_gff(gff_file):
	gene_dict = {}
	with open(gff_file, 'r') as gff:
		for line in gff:
			if not line.startswith('#'):
				fields = line.strip().split('\t')
				if len(fields) >= 9 and fields[2] == 'gene':
					#gene_id = fields[8].split(';')[0].replace('ID=', '')
					gene_id =extract_gene_id(fields[8])
					gene_dict[gene_id] = (fields[0], int(fields[3]), int(fields[4]))
	return gene_dict
	
#---------------------------------------------------------------------------------------------------------------------------
#This function
def find_flanking_genes(scrm_gene_list, gene_dict, scrm_chr, scrm_start, scrm_end):
	left_flank_gene = None
	right_flank_gene = None
	#print(scrm_gene_list)
	for gene in scrm_gene_list:
		if gene in gene_dict:
			#print('present in gene_dict')
			gene_chr, gene_start, gene_end = gene_dict[gene]
			if gene_chr == scrm_chr:
				# SCRM is completely contained within the gene
				if scrm_start >= gene_start and scrm_end <= gene_end:
					left_flank_gene = right_flank_gene = gene
				elif scrm_start < gene_start < scrm_end:  # Overlaps on the left side
					left_flank_gene = gene
				#elif scrm_start < gene_start < scrm_end and (left_flank_gene is not None) and right_flank_gene is None:  #  both genes are overlapping on the left side
				#	right_flank_gene = gene
				elif scrm_start < gene_end < scrm_end and right_flank_gene is None:  # Overlaps on the right side
					right_flank_gene = gene
				#elif scrm_start < gene_end < scrm_end and (right_flank_gene is not None) and left_flank_gene is None:  # both genes are overlapping on the right side
				#	left_flank_gene = gene
				#
				elif gene_end > scrm_start and (left_flank_gene is None or gene_end > gene_dict[left_flank_gene][2]):
					left_flank_gene = gene
				elif gene_start < scrm_end and (right_flank_gene is None or gene_start < gene_dict[right_flank_gene][1]):
					right_flank_gene = gene
		else:
			print(gene,'not found in gene_dict')
	return left_flank_gene, right_flank_gene
	#log_file.close()

#----
def update_scrmsOut(scrmsOut_file, gene_dict,path_log_file):
	updated_lines = []
	with open(scrmsOut_file, 'r') as scrmsOut, open(path_log_file, "a") as logfile:
		for line in scrmsOut:
			#print(line)
			fields = line.strip().split('\t')
			genes_col6 = fields[5].split(',')
			genes_col11 = fields[10].split(',')
			# save_log line for any  case where right or left flanks have more than one gene
			if len(genes_col6) > 1 or len(genes_col11) >1:
				logfile.write(line)

			if len(genes_col6) == 2 and set(genes_col6) == set(genes_col11):
				#print(line)
				#logfile.write(line)
				scrm_chr, scrm_start, scrm_end = fields[0], int(fields[1]), int(fields[2])
				left_flank_gene, right_flank_gene = find_flanking_genes(genes_col6, gene_dict, scrm_chr, scrm_start, scrm_end)
				
				#print("Here------------------------")
				#print('left',left_flank_gene)
				#print('right',right_flank_gene)
				if left_flank_gene and right_flank_gene:
					#print('updated',left_flank_gene,'-',right_flank_gene)
					fields[5] = left_flank_gene
					fields[6] = left_flank_gene
					fields[10] = right_flank_gene
					fields[11] = right_flank_gene
				else:
					#print("missing a flanking gene")
					#print(left_flank_gene)
					#print(right_flank_gene)
					#print("updated",genes_col6[0],'-',genes_col11[1])
					fields[5] = genes_col6[0]
					fields[6] = genes_col6[0]
					fields[10] = genes_col11[1]
					fields[11] = genes_col11[1]

			#if right flanking gene has more than one gene and left has just 1 gene, we are gonna loose the third one
			elif len(genes_col6) > 1 and set(genes_col6) != set(genes_col11) and len(genes_col11) ==1:
				#logfile.write(line)
				fields[5]=fields[5].split(',')[0]
				fields[6] = fields[6].split(',')[0]
			# if left flanking gene has more than one gene and right has just 1 gene, we are gonna loose the third one
			elif len(genes_col11) > 1 and set(genes_col6) != set(genes_col11) and len(genes_col6) ==1:
				#logfile.write(line)
				fields[10]=fields[10].split(',')[0]
				fields[11] = fields[11].split(',')[0]
			elif len(genes_col11) > 1 and set(genes_col6) != set(genes_col11) and len(genes_col6) >1:
				fields[5]=fields[5].split(',')[0]
				fields[6] = fields[6].split(',')[0]
				fields[10]=fields[10].split(',')[0]
				fields[11] = fields[11].split(',')[0]

			##if both right and left flanking gene has more than one gene then we are gonna loose everything except 2
			elif len(genes_col6) > 2 and set(genes_col6) == set(genes_col11):
				#logfile.write(line)
				fields[5]=fields[5].split(',')[0]
				fields[6] = fields[6].split(',')[0]
				fields[10]=fields[10].split(',')[1]
				fields[11] = fields[11].split(',')[1]
			updated_lines.append('\t'.join(fields))

	return updated_lines

#----
def write_updated_scrmsOut(scrmsOut_file, updated_lines):
	with open(scrmsOut_file, 'w') as scrmsOut:
		scrmsOut.write('\n'.join(updated_lines))
		scrmsOut.write('\n')  #adds newline at end of final line so file concatenate properly later on
	#pathFinal=os.path.abspath(scrmsOut)
	#path=os.path.abspath(peaksToScrmsName)
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name==scrmsOut_file:
				pathFinal=os.path.abspath(os.path.join(root,name))
	return pathFinal

#############################################-------MAIN FUNCTION-----##########################################################

def main():
	#some of the variables are declared global because these are being used in main program as well as some functions
	global d  
	global d2
	global countd
	global subdirectory
	global patternRecovery
	global totalNumberOfCrmsKnownToCauseExpression
	
	totalNumberOfCrmsKnownToCauseExpression=0
	t=1 #count of sets done
	#Temporary directory in which all the intermediate files will be moved
	subdirectory='tmp' 

	#command line parsing files from users
	parser=argparse.ArgumentParser()
	parser.add_argument('-so','--scrmJoinedOutputFile',help='Scrmshaw Output file concatenated ',required=True)
	parser.add_argument('-num','--numOfScrms',help='Number of Scrms to start from, default is 5000',default=5000)
	parser.add_argument('-topN','--topNmethod',help='point of extraction for top peaks i.e Elbow or Median or None on amplitiude curve',default='elbow')
	parser.add_argument('-gff','--gffFile',help='GFF file used for running SCRMshaw ',required=True)
	args = parser.parse_args()
	scrmJoinedOutputFile=args.scrmJoinedOutputFile
	gffFile=args.gffFile
	numOfScrms=args.numOfScrms
	numOfScrms=int(numOfScrms)
	topNmethod=str(args.topNmethod)
	my_path=os.getcwd()
	if not os.path.isdir(my_path+'/'+subdirectory):
		os.makedirs(my_path+'/'+subdirectory)
	
	#checking if topNmethod value is correct from the options 
	topNmethod = topNmethod.lower()
	if topNmethod not in ['none', 'elbow', 'median']:
		print('incorrect value of topN for extracting top peaks from amplitude curve')
		exit()
	#iterating through each keys (different training sets) of the three method's dictionary:
	methods=['imm','hexmcd','pac']
	
	#Creating list to iterate through three methods 
	methods_val=[None,None,None]
	num=0
	x=0
	
	#gff file path
	gffFile=os.path.abspath(gffFile)
	#saving all the genes coorddiates from GFF
	gene_dict = parse_gff(gffFile)
	#pprint.pprint(gene_dict)
	#Parsing the output file into separate files for each training set and each method via creating three dictionaries for each method: keys being the names of training sets associated with that method 
	scrmJoinedOutputFile=os.path.abspath(scrmJoinedOutputFile)
	methods_val[0],methods_val[1],methods_val[2]=parse_output(scrmJoinedOutputFile,35000)

	#this loop is used to iterate through three methods dictionaries
	for num in range(len(methods)):
	
		#removing empty string from list
		while '' in methods_val[num]:
  	  		methods_val[num].remove('')
		
		print("Now method:"+methods[num])
		print("methods_Val of num")
		print(methods_val[num])
		#pprint.pprint(methods_val[num])
		
		#loop is used to iterate through all the training sets in the method		
		for x in methods_val[num]:
			TSET=x
			method= methods[num]
			print("Tset: "+TSET)
			print("Method: "+method)
			
			#individual training set's full length scrmshaw predictions file
			file= method+"_"+TSET+"_fullLength.bed"

			#getting the path of this file
			for root, dirs, files in os.walk(os.getcwd()):
				for name in files:
					if name==file:
						#print(name)
						my_path=os.path.abspath(os.path.join(root,name))
		
			#print("Path of your file: is"+ my_path)
			fileFullLengthPath=my_path		
		
			#extract num of scrms from offsets combined output
			extractedScrmsPath,minScore= extract_topN_scrms_scoreCurve(fileFullLengthPath,numOfScrms,method,TSET)
						
			#calculating the min score to use as cutoff for macs from above file .
# 			process = subprocess.Popen(["sort",'-V',"-k4",extractedScrmsPath], stdout=subprocess.PIPE)
# 			output = process.communicate()[0]
# 			tmp=output.split('\n')
# 			tmp1=tmp[0].split('\t')
# 			minScore= tmp1[3]
# 			print(minScore +' minScore')
# 			with open("sortedTEST",'w') as f:
# 				processRev = subprocess.Popen(["sort",'-V','-r',"-k4",extractedScrmsPath], stdout=subprocess.PIPE)
# 				outputRev = processRev.communicate()[0]
# 				f.write(outputRev)
				
			#sorting the extracted individual offset combined file.
			#first converting it to bed
			extractedScrmsPathBED=bed_conversion(extractedScrmsPath)
			sortedFileName='sorted_'+str(numOfScrms)+'.'+method+'_'+TSET+'.bed'
			sortedFilePath=bedtools_sorting(extractedScrmsPathBED,sortedFileName)
			sumOfScoreFileName=str(numOfScrms)+'.'+'original_'+method+'_'+TSET+'.bdg'
			
			# 	calculating sum of score for each 10bp window 
			sumOfScoreFilePath=sum_of_score(sortedFilePath,sortedFileName,sumOfScoreFileName)

			#sorting the above created file
			sumOfScoreFilePathBED=bed_conversion(sumOfScoreFilePath)
			sumOfScoreFileName='sorted_'+sumOfScoreFileName
			sortedSumOfScoreFilePath=bedtools_sorting(sumOfScoreFilePathBED,sumOfScoreFileName)
			
			#filling out the missing windows across the genome with the score with 0.0 (required step for MACs)
			wholeGenomeFileName='whole_genome_coverage_'+sumOfScoreFileName
			wholeGenomeFilePath=filling_missing_coords(sortedSumOfScoreFilePath,wholeGenomeFileName)
			
			#calling MACs
			macsOutputName=TSET+'_'+method
			macsOutputPath=callingMACS(wholeGenomeFilePath,macsOutputName,minScore)
			

			#converting MACS output narrowPeak to BED file (like SCRMshaw prediction file)
			peaksToScrmsName='peaks_'+macsOutputName
			peaksToScrmsPath=peaksToScrms(macsOutputPath,peaksToScrmsName,TSET)
	

			#intersecting peaks output file with the SCRMshaw prediction file and get the SCRMshaw score for each of the peak.

			#converting both to pybed version to perform bedtools functions
			peaksToScrmsPathBED=bed_conversion(peaksToScrmsPath)
			extractedScrmsPathBED=bed_conversion(extractedScrmsPath)
			#intersecting
			intersectedFileName='Intersected_ScrmsAndPeaks_'+TSET+'_'+method+'.bed'
			intersectedFilePath=intersect_peaks_and_scrms(peaksToScrmsPathBED,extractedScrmsPathBED,intersectedFileName)
			
			#reading intersected file as pandas dataframe to sort based on amplitude	
			numOfPeaks=peaksToScrmsPathBED.count()
			finalPeaksFileNameA='scrmshawOutput_peaksCalled_'+TSET+'_'+method+'_'+str(numOfPeaks)+'_peaks.bed'
			finalPeaksFilePath,numOfpeaks2=sortAndRank_basedOnAmplitude(intersectedFilePath,finalPeaksFileNameA)
			print('All number of peaks on the amplitude curve for the set '+TSET+'_'+method+': '+str(numOfpeaks2))
			
			#adding chunk from filtering Flanking Genes--
			#finalPeaksFileName='scrmshawOutput_peaksCalled_'+TSET+'_'+method+'_'+str(numOfPeaks)+'_peaks.bed'
			
			#log_file = open("log_flankingMoreThanOneGenes.txt", "w")
			#log_file_path = "log_flankingMoreThanOneGenes_"+TSET+'_'+method+".txt"
			log_file_path = "log_flankingMoreThanOneGenesFromAllSets.txt"
			log_file = open(log_file_path, "a")
			path_log_file = os.path.abspath(log_file_path)
			
			
			if (topNmethod != 'none'):
				
				finalExtractedPeaksFileNameB='scrmshawOutput_peaksCalled_'+TSET+'_'+method#+'_peaks.bed'
				finalPeaksFilePath,numOfpeaks3,finalExtractedPeaksFileNameB2=extract_topN_scrms_amplitudeCurve(finalExtractedPeaksFileNameB,finalPeaksFilePath,topNmethod)
				print('Top number of peaks on the amplitude curve for the set '+TSET+'_'+method+': '+str(numOfpeaks3))
				shutil.move(finalPeaksFileNameA, 'tmp/')
				updated_lines = update_scrmsOut(finalPeaksFilePath, gene_dict,path_log_file)
				finalPeaksFilePath=write_updated_scrmsOut(finalExtractedPeaksFileNameB2, updated_lines)
				#shutil.move(finalExtractedPeaksFileNameB, 'tmp/')
			else:
				updated_lines = update_scrmsOut(finalPeaksFilePath, gene_dict,path_log_file)
				finalPeaksFilePath=write_updated_scrmsOut(finalPeaksFileNameA, updated_lines)
			# adding chunk from filtering Flanking Genes--
# 			finalPeaksFileName='scrmshawOutput_peaksCalled_'+TSET+'_'+method+'_'+str(numOfPeaks)+'_peaks.bed'
# 			
# 			log_file = open("log_flankingMoreThanOneGenes.txt", "w")
# 			path_log_file=os.path.abspath(log_file)
# 			log_file_path = "log_flankingMoreThanOneGenes.txt"
# 			log_file = open(log_file_path, "w")
# 			path_log_file = os.path.abspath(log_file_path)
#
#			updated_lines = update_scrmsOut(finalPeaksFilePath, gene_dict,path_log_file)
#			finalPeaksFilePath=write_updated_scrmsOut(finalPeaksFileName, updated_lines)
			log_file.close()
			#shutil.move(finalExtractedPeaksFileNameB, 'tmp/')
			#moving the extra files to tmp directory
			for root, dirs, files in os.walk(os.getcwd(),topdown=False):
				for name in files:
					if name==TSET+'_'+method or name=='sumScoredsorted_'+str(numOfScrms)+'.'+method+'_'+TSET+'.csv':
						shutil.move(name, 'tmp/')
		
			
main()										
