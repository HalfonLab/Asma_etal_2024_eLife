#!/usr/bin/perl -w
use strict;
use Getopt::Long;


#Author - Kushal Suryamohan and Marc S. Halfon
# Creation date: October 6th, 2014

# Utility program to post-process SCRMshaw data
# Creates a BED-compatible file with the N top-ranked hits for each training set and method from a SCRMshaw run
#List all contents of hits folder for each SCRMshaw method, and parse the ranked hits file within each tissue and extract FASTA Headers
#Input - Path to SCRMshaw results directory, number of desired hits per training set and method, Name of output file that will be of the following format - 
##Output:
#Chromosome	CRM-start	CRM-end	Ng1	5'or3'	dist	Ng2	5'or3'	dist	local rank	global rank	Training set name

#Revision history:
#version 2: Marc Halfon, December 23 2016
#revised to ensure that all output has 14 fields, or warn if otherwise
#cleaned up code a bit

#version 3: MSH, Jan 3 2017
#add code to include source for each top hit

#version 3.1: MSH May 11 2017
#implement Getopt::Long, minor bug fixes

#version 3.3: MSH July 11, 2017
#adds 'restict_chr' option and subroutine

#version 3.4: MSH Oct 26, 2017
#fixes hard coding 'mapping' in training set name for appending file name

#version 3.5 MSH Dec 20 2017
#add a column with the rank of each hit

my $version = 3.5;   #program version

Usage() unless @ARGV;

#define options
my $folder = "";  #path to SCRMshaw output directory
my $outfilename = ""; #name of file for output
my $N = 250; #number of hits for each training set/method
my $help = 0;
my $show_version = 0; 
my $restrict_chr = ""; #file for restricting to specific chromosomes
	#the 'restrict_chr' file should be in the form of a regular expression containing
	#the chromosomes that are ALLOWED. Anything else will be excluded if this option is
	#set

GetOptions(	"directory=s" => \$folder,
			"outfile=s" => \$outfilename,
			"n=i" => \$N,
			"r=s" => \$restrict_chr,
			"version" => \$show_version,
			"help" => \$help,
		);	

die "$0 --version $version\n" if $show_version;
Usage() if $help;
Usage() unless ($folder && $outfilename && $N);

my $keep_chr = "";
if ($restrict_chr){		#gets regular expression for checking chromosomes later
	$keep_chr = Chr_Restriction();
}	

my @files = ();
#my $folder = $ARGV[0];
#my $outfilename = $ARGV[1];
#my $N = $ARGV[2];

print "$0 v$version\n";

print "$folder \n";
open (OUT, '>', $outfilename) || die "Couldn't open outputfile\n";

list_recursively($folder);

close(OUT);

#############Subroutine#####################

#List contents of a folder recursively and find files with suffix ".hits.ranked" and parse its contents

############################################

sub Usage{
	die "\nUsage: $0 -d directory -o out_file_name -n number_of_hits [default 250] -r restrict_chr [file, optional] -v show_version
		\"directory\" is directory containing SCRMshaw results
		\"number_of_hits\" is the desired number of hits for EACH trainging set and method (i.e., imm, hexmcd, pac) and defaults to 250
		\"restrict_chr\" is a file containing a regular expression of allowable chromosome names
		\n\n"
}		

############################################

sub Chr_Restriction {
	
	open(CRES, '<', $restrict_chr) or die "Can't open chr restriction file: $!";
	my $re = "";
		while(<CRES>){
			$re = $_;
			chomp($re);
		}
		return $re;
	}		

############################################

sub list_recursively {

	my @tmparray = (); #added 12/23/2016 #hold NULL fields to keep lines equal length
	my $linecount=0;   #added 12/23/2016 #keeps track of output line numbers

	my($directory) = @_; #@_ contains values passed as arguments to subroutine from main program
	my @files = ();
	unless (opendir(DIR, $directory))
	{
		print "Cannot open directory $directory!\n";
		exit;
	}
	#Read directory, ignoring special entries "." and ".."
	@files = grep (!/^\.\.?$/,readdir(DIR));
	closedir(DIR);
	
	#If file, print name, else if a directory/subdirectory, then recurse through the contents of this subdir/dir
	
	## ****NOTE- prepend directory name to gain access
	foreach my $file(@files)
	{
		if (-f "$directory/$file") #if /$directory/$file is a file
		{
			#print "File is $file\n";
			#Grab only files with suffix ".hits"
			if ($file =~ /.hits.ranked$/ ) 
			{ #found file, parse it now
				#print "File is $file\n";
				my $F = $directory."/".$file; #Append directory path to gain access
				print "processing $F\n";
				open (IN, '<', $F) || die "Cannot open $file to parse\n";
				my %hash = ();
				my $key;
				my $top = 0;
				while(<IN>)
				{
    					chomp;
					if (/^>/) {
						$key = $_; #FASTA Header
						my $line = $key; 
						$line =~ s/^>//g;   #remove > at beginning of line	
						my @fields = split /[:,\s\/]+/, $line; #Split the different fields
						
						if ($restrict_chr){
							#my $keep_chr = Chr_Restriction();
							next unless $fields[0] =~ /$keep_chr/;
						}	
						
						my $end = $fields[1] + 500;
						splice @fields, 2, 0, $end; #Inserting End coordinate into line 
						
						
						$linecount++;
						
						####fixing number of fields per line####
						if (scalar(@fields) < 14){
							my $k = 14 - scalar(@fields);
 							for my $i (0..($k-1)){
 								$tmparray[$i] = "NULL";
 							}
							warn "number of fields < 14, appending $k fields: line $linecount\n";	
						}	
							
						if (scalar(@fields) > 14){
							warn "PROBLEM (NOT FIXED): too many fields at $linecount\n";
						}	
							
						push(@fields, @tmparray); 
														
						@tmparray = ();

						#append directory and file information
						#parse file name:
						$file =~ /(\w.*)\.hits/;
						my $training = $1;
						#parse directory, kludgy but simple:
						my $method = "undef";
						if($directory =~ /hits\/imm/){$method = "imm";}
						if($directory =~ /hits\/pac/){$method = "pac";}
						if($directory =~ /hits\/hexmcd/){$method = "hexmcd";}
						
						my $rank = $top + 1;				
						print OUT join("\t", @fields, $training, $method, $rank), "\n" ;
													
						if($top==($N-1)){   #subtract 1, else will go too many times
							last;
						}
						$top+=1;
					} 
					else {
						$hash{$key} .= $_;
					}
				}
			
			}
		}
		elsif (-d "$directory/$file")
		{
			#Call list_recursively sub - recursive call
			list_recursively("$directory/$file");
		}
	}	
}

exit;

