#!/bin/bash

i=0
for task in task_offset_0_1 task_offset_10_2 task_offset_20_3 task_offset_30_4 task_offset_40_5 task_offset_50_6 task_offset_60_7 task_offset_70_8 task_offset_80_9 task_offset_90_10 task_offset_100_11 task_offset_110_12 task_offset_120_13 task_offset_130_14 task_offset_140_15 task_offset_150_16 task_offset_160_17 task_offset_170_18 task_offset_180_19 task_offset_190_20 task_offset_200_21 task_offset_210_22 task_offset_220_23 task_offset_230_24 task_offset_240_25

do 
	echo $task
	echo $i
	perl /projects/academic/mshalfon/Scripts/Generate_top_N_SCRMhits.pl -d $task -n 5000 -o scrmshawOutput_offset_$i.5000scrms
	let i=i+10


done

#concatenating individual offset orig files to one file, which is going to be used as an input for the next pipeline
cat scrmshawOutput_offset_* > scrmshawOutput_offset_0to240.bed

echo "Moving Temporary files"
#making directory to move temporary created files i.e individual 0 to 240 offset files
mkdir scrmsIndividualHits_0to240offset

mv scrmshawOutput_offset*5000scrms scrmsIndividualHits_0to240offset/

