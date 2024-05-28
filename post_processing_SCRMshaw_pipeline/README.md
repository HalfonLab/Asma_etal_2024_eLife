

################################################################################												
   								  			
	`````````SCRMshaw post processing pipeline`````````````
				Halfon Lab			
  		     Date: Sept 2018
	 	     Updated: May 2024
  		      										
################################################################################

INTRODUCTION:
This script is written to do post processing steps on the SCRMshaw multiple offsets output including but not limited to use of peaks calling algorithm.


	1. INPUT
	2. USAGE
	3. PARAMETERS
	4. OUTPUT


1. INPUT:

The only file required for this script is the concatenated output(set of predictions) of SCRMshaw-HD (10 bp offset) and the number of predictions to take if more than 5000 (default), and gff annotation file. The concatened output can be obtained using the script "concatenatingOffsetsResults.sh". However, a better solution is to run the complete post-processing pipeline using the script "post_processing_complete.sh", which will do the concatenation, run the post processing, and clean up the directories.

2.  USAGE:

For running on the CCR cluster at University at Buffalo, use the Slurm script "Slurm_postProc-complete.sh" as follows:

sbatch Slurm_postProc-complete.sh <path-to-GFF-annotation-file>

If running elsewhere, note that the pipeline requires MACs2 installed on the system.It will also check the status of MACs2 installation or will terminate with the error.
The following Python modules are also required to run this script. Please make sure these modules have already been properly installed and are recognizable:
pybedtools, statistics, scipy, numpy, pandas, csv, subprocess, shutil and itertools.

Following is an example of command line execution.
>python postProcessingScrmshawPipeline.py  -so SCRMshawConcatenatedOutputFile  -num 5000  -topN Median -gff GFFfile

3. PARAMETERS:

	-so	<str>	The output file for concatenation shell scripySCRMshaw\
	-num	<int>	Number of predictions to processs
	-topN	<str>	Median/Elbow/All
	-gff	<str>	GFF file path

5. OUTPUT:

Depending on user's input, if the input consists of multiple training sets the pipeline will create peaks file for each individual training set.
And a temporary folder which will have the files created during execution.

