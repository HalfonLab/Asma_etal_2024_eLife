#!/usr/bin/perl -w

################################################################
# preflight-scrmshaw.pl  v1.5
#
# (c) Marc S. Halfon September 2019
#
# use before running SCRMshaw to check that annotation gff and
# genome FASTA files are properly formatted
#
# check the 'name/id' output from the gff to see if this will
# be suitably parsed for desired SCRMshaw output
###############################################################


=head1 Description

This file checks that a set of GFF and FASTA files are compatible with SCRMshaw. It generates a log file confirming compatibility and/or describing any issues. It also provides example output of SCRMshaw's 'gene' and 'exon' files.



=head1 Usage

perl preflight-scrmshaw.pl --gff gff3_annotation_file --fasta genome_FASTA_file [options] 

	--gff/-g <str>		A gff3-formatted genome annotation file. Mandatory argument.
	--fasta/-f <str>	A FASTA-formatted genome sequence. Mandatory argument. 
	
	--suppress_log/-s	Suppresses log file generation and sends all output to STDOUT
	--print_diff/-p		Takes value 'gff' or 'fasta', prints seqid's found in one file but not the other
	--number/-n <int>	Number of lines of example gff3.pl output to print, default=5
	--help/-h		Prints this help text

=cut

use strict;
use Getopt::Long;
use File::Basename;
#use Data::Dumper;


my $Gff = "";
my $genome = "";
my $out = 0;
my $print_diff = "";
my $checkcount = 5;
my $help = "";
my $stoppoint = 0;

GetOptions(
    "fasta=s"=>\$genome,
    "suppress_log"=>\$out,
	"gff=s"=>\$Gff,
	"print_diff:s"=>\$print_diff,
    "number:i"=>\$checkcount,   #number of lines to output, default = 5
    "endpoint:i"=>\$stoppoint,  # sets finite number of gff lines to parse
    "help"=>\$help,
   ); 

## Check exit conditions
die `pod2text $0` if (!$genome || !$Gff || $help);

my($day, $month, $year) = (localtime)[3,4,5];
$month = $month+1;
$year = $year + 1900;
my $datestring = "$month-$day-$year";


# print $fh "$Gff\n";
# print $fh "$genome\n";


my($gffname, $gffpath, $gffsuffix) = fileparse($Gff);
my($genomename, $genomepath, $genomesuffix) = fileparse($genome);

my $logfilename = "preflight_$genomename" . "_$gffname" . "_$datestring";


## global hashes and arrays
my %seqid = (); #hash of seqids
my %seqidHOA = (); #hash of arrays for full gene table
my %geneid = ();
my %exonid = ();
my %chr = (); #for fasta file chromosomes
my %annot_types = (); #hash of 'types' from 3rd field of gff
my %headerchr = (); #hash for chr ids from ##sequence-region lines
my $tmpchrsize = 0; # tmp store of scaffold/chr length

my @genearray = ();
my @exonarray = ();

#open out file: see https://www.itworld.com/article/2816641/redirecting-standard-error-in-perl.html
my $fh;
if ($out){
	open ($fh, ">&STDOUT") or die "Problem directing to STDOUT for writing: $!";
} else {
	open ($fh, ">$logfilename") or die "Can't open $logfilename for writing: $!";
}		

print $fh "Validation report for $gffname, $genomename: $datestring\n\n";


##======================== validate the GFF file ====================##
#note that this is not a true GFF3 validation but simply checks
#that the format is good enough for what we need for SCRMshaw

# Open gff file
open(IFILE, $Gff) or die("Couldn't open gff file: $Gff for reading\n");

my $gfflinecount = 0;
my $numfielderr = 0;
my $whiteerr =0;
my $coorderr = 0;
my $stranderr = 0;
my $typeerr = 0;
my $hasgene = 0;
my $hasexon = 0;

warn "\nprocessing GFF3 file: $Gff\n";

while(<IFILE>){
	$gfflinecount++;
	
	##############################################
	#debug (hidden command line option)
	if (($stoppoint > 0) && ($gfflinecount==$stoppoint))
		{ print $fh "\n**stopping at line $gfflinecount**\n\n"; 
		last;
	}
	##############################################
		
	
	chomp;
    
    last if $_ eq '##FASTA'; #we don't want this part of the file
    

### this part updated v1.3 to see to solve problem of RefSeq format GFF where  "##sequence-region" != seqID
   if ($_ =~ /##sequence-region/){
    	my @line = split;
   
   		$headerchr{$line[1]}{seen} = 0; #establishes the hash key=>value
   		$headerchr{$line[1]}{length} = $line[3];  #gets the length
    	$tmpchrsize = $line[3]; #holds the length in tmp for when seqID is different
    	next;
    }	
###
 
    next if $_ =~ /^#/;
	
	#my ($seqid-0,  $source-1, $type-2, 
	#		$start-3,  $end-4,    $score-5,	
	#		$strand-6, $phase-7,	$attr-8)	
	
	my @line = split('\t', $_);  #assumes tab-delimited file
	#but scrmshaw code uses any white space in gff3.pl
	
	
### this part this part updated v1.3 to see to solve problem of RefSeq format GFF where  "##sequence-region" != seqID
# assumes that the feature lines will immediately follow the ##sequence-region if not
# all grouped together at the beginning

	$headerchr{$line[0]}{length} = $tmpchrsize unless $headerchr{$line[0]};
	$headerchr{$line[0]}{seen}++;
###	
	
	
	#track the 'types' used in this GFF file:
	$annot_types{$line[2]}++;
		
	#only want 'gene' or 'exon' types:
	next unless ( ($line[2] eq 'gene') || ($line[2] eq 'exon') );
	
	#first check: right number of fields for gff3
	if (scalar(@line) != 9) {
		print $fh "GFF3: WARNING: number of fields is incorrect at line $gfflinecount\n";
		$numfielderr++;
		print $fh "$_\n";
		
	}
	
	#next check: does the SCRMshaw gff3.pl avoid white-space issues?
	 my ($chr_id,  $source, $type, $start,  $end,  $score,  $strand, $phase, $attr)  = split('\s+', $_);
	 
	 if(($chr_id ne $line[0]) ||
	 	($source ne $line[1]) ||
	 	($type ne $line[2])	  ||
	 	($start != $line[3])  ||
	 	($end != $line[4])    ||
	 	($score ne $line[5])  ||
	 	($strand ne $line[6]) ||
	 	($phase ne $line[7]) ) {
	 	#($attr ne $line[8]) )   {
	 		print $fh "GFF3: WARNING: white space problem at line $gfflinecount\n";
	 		$whiteerr++;
	 	}	
	
	
	#next check: numeric start and stop coordinates
	if (($line[3] =~ /\D/) || ($line[4] =~ /\D/)){
		print $fh "GFF3: WARNING: non-numeric data in start/stop fields at line $gfflinecount\n";
	$coorderr++;
	}	
	
	#check strand:
	if ($line[6] !~ /^[+-]$/) {	
		print $fh "GFF3: WARNING: strand field is not +/- only at line $gfflinecount\n"; 
	$stranderr++;
	}		
	
	#get seqid (chr)
	$seqid{$line[0]}++ unless $seqid{$line[0]};
			
	#get gene and exon data, including storing all gene BED-like to get distances later
	if ($line[2] eq 'gene'){
		$hasgene++;
		
		#this part of code from SCRMshaw 'gff3.pl'
		my $id = "";
		if($line[8] =~ /;/){
			$id  = $1 if $line[8] =~ /[ID|Name]=(\S+?);/;
		} else {
			$id  = $1 if $line[8] =~ /[ID|Name]=(\S+)$/;
		}
		
		$geneid{$id}++;
		
		my $tmpline = "$line[0]\t$line[3]\t$line[4]\t$line[6]\t$id\n";
		
		push( @{$seqidHOA{$line[0]}}, $tmpline) ; #all the gene lines sorted by chr
		#for analysis, might be useful to add to end of above line: 'unless $hasgene > $checkcount'
		
		if($hasgene <= $checkcount){
			
			push(@genearray, "original: $line[8]");
			push(@genearray, "name/id: $id");
			push(@genearray, "parsed: $tmpline");
			#print $fh "debug: here $hasgene\n";
		}	
	
	} elsif ($line[2] eq 'exon'){
		$hasexon++;
	
		#this part of code from SCRMshaw 'gff3.pl'
		my $id = "";
		if($line[8] =~ /;/){
			$id  = $1 if $line[8] =~ /[ID|Name]=(\S+?);/;
		} else {
			$id  = $1 if $line[8] =~ /[ID|Name]=(\S+)$/;
		}
		
		$exonid{$id}++;
		
		if($hasexon <= $checkcount){
			my $tmpline = "$line[0]\t$line[3]\t$line[4]\t$line[6]\t$id\n";
			push(@exonarray, "original: $line[8]");
			push(@exonarray, "name/id: $id");
			push(@exonarray, "parsed: $tmpline");
		}	
	
	
	} else {
		print $fh "GFF: WARNING: unknown type at line $gfflinecount\n";
		$typeerr++;
		next; #should never really get here, but...
	}	
		

}
close(IFILE);

#want to make sure that have both 'gene' and 'exon' types (updated v1.5):
if (undef $annot_types{'gene'}){
	print $fh "GFF: WARNING: No 'gene' set as a feature type\n";
	$typeerr++;
}

if (undef $annot_types{'exon'}){
	print $fh "GFF: WARNING: No 'exon' set as a feature type\n";
	$typeerr++;
}

if($numfielderr == 0){
	print $fh "GFF3: all lines have correct number of fields\n";
}

if($whiteerr == 0){
	print $fh "GFF3: no white space problems in fields\n";
}
	
if($typeerr == 0){
	print $fh "GFF3: no type errors with gene/exon\n";
}		
	
if($coorderr == 0){
	print $fh "GFF3: all lines have numeric start/stop coordinates\n";
}
if($stranderr == 0){
	print $fh "GFF3: all lines have correct strand designations\n";
}	

print $fh "GFF3: number of genes: $hasgene\n";
print $fh "GFF3: number of exons: $hasexon\n";	

print $fh "\n------------------------------------------------------------\n";
print $fh "\nGENE DATA\n";
foreach my $elem (@genearray){
	print $fh "$elem\n";
}	
	
print $fh "\n------------------------------------------------------------\n";print $fh "\nEXON DATA\n";
foreach my $elem (@exonarray){
	print $fh "$elem\n";
}	
print $fh "------------------------------------------------------------\n";



##======================== validate the FASTA file ====================##
#again, not a complete validation but checks what we need

#check sequence ATGCNatgcn
#get headers into hash

my $currentheader = "";
my $seqerr1 = 0;
my $seqerr2 = 0;

open(GENEIN, "<$genome") or die "Can't open genome file: $!";

warn "processing genome file: $genome\n";

while(<GENEIN>){
	chomp;
	if ($_ =~ /^>/){   #FASTA header line
		$currentheader = $_;
		
		if($currentheader =~ /\|/){   #pipes cause problem later in scrm.pl with call to kmers.pl
			$seqerr1++;
		}
		
		my @tmpchr = split('\s+',$currentheader);
		$tmpchr[0] =~ s/>//;  #remove the leading '>'
		$chr{$tmpchr[0]}++;
		
		} else {
		if ($_ =~ /[^ATGCNatgcn]/){   # (^ negates the class)
			$seqerr2++;
		}	
	}
}
close(GENEIN);
if($seqerr1 > 0 && $seqerr2 > 0){
	if($seqerr1 > 0){
		print $fh "Problem with FASTA header: need to remove pipes ('|'), suggest replace with underscore ('_')\n";
	}
	if($seqerr2 > 0){
		print $fh "FASTA: improper character in sequence $currentheader\n";
	}
} else {
	print $fh "FASTA: all sequences have proper characters\n\n";
}




##================= compare the seqids between the files ==============##

#seqids in FASTA but not in GFF: (from PERL Cookbook)
my @fastaonly = ();
foreach (keys %chr){
	push(@fastaonly, $_) unless exists $seqid{$_};
}

my @gffonly = ();
foreach (keys %seqid){
	push(@gffonly, $_) unless exists $chr{$_};
}

print $fh "------------------------------------------------------------\n";
my $not_in_gff = scalar(@fastaonly);
if( $not_in_gff > 0){
	print $fh "There are $not_in_gff seqids not in the GFF file (use the '-p gff' option to output)\n";
	print $fh map {$_, "\n"} @fastaonly if ($print_diff eq 'gff');
} else {
	print $fh "All seqids in FASTA are also in GFF\n";
}	
	
my $not_in_fasta = scalar(@gffonly);
if ($not_in_fasta > 0){
	print $fh "There are $not_in_fasta seqids not in the FASTA file (use the '-p fasta' option to output)\n";
	print $fh map {$_, "\n"} @gffonly if ($print_diff eq 'fasta');
} else {
	print $fh "All seqids in GFF are also in FASTA\n";	
}					
print $fh "See below for list of valid seqids and lengths\n";
print $fh "\n";
print $fh "------------------------------------------------------------\n";


##======================== get intergenic differences ====================##
# in this implementation, this is just an approximate average of the
# intergenic differences because we don't fully account for nested genes. 
# We are assuming this number is small
# relative the whole and won't have a big effect.
# We calculate as follows:
# 		normal case: dist = start[n] - end[n-1]   (treat the first end as zero)
# 		overlapping: dist = start[n+1] - end[n] => so don't add in distance, but do update end[n-1]
# 		nested: dist = start[n+1] - end[n-1] => so don't add in distance, don't update end[n-1]



# "Shwartzian transform":
# must remember: effectively, first index of relevant array to sort is now 1 not 0
# my @sorted = map {$_->[0]}
#              sort { $a->[1] cmp $b->[1] ||
#                     $a->[2] <=> $b->[2] ||
#                     $a->[3] <=> $b->[3] }
#              map {chomp;[$_,split(/,/)]} <DATA>;

my $totaldist = 0;
my $maxdist = 0;
my $totalgenes = 0;
my $totaldistgenes = 0;
my $negativedist = 0;
my $genelength = 0;

print $fh "\nIntergenic distance information\n\n";

foreach my $key (keys %seqidHOA){
	my @sorted = ();
	@sorted = map {$_->[0]}
             sort { $a->[2] <=> $b->[2] }
             map {chomp;[$_,split(/\t/)]} @{$seqidHOA{$key}};

	#chr[0], start[1], end[2]

	my $oldend = 0;
	
	foreach my $gene (@sorted){
		
		$totalgenes++;
		
		my @line = split('\t', $gene);
		my $dist = $line[1] - $oldend;
			#print $fh "eval $line[2], $oldline\n";	
					
		if($dist < 0){
			$negativedist++;
		} else {
			$totaldist += $dist;
			$maxdist = $dist if $dist > $maxdist; #track the maximum distance
			
			
			print $fh "huge: $gene\n" if $dist > 700000;
			
			
			$totaldistgenes++;   #only add if used for taking average
		}	
		
		$genelength += ($line[2] - $line[1]);
		
		$oldend = $line[2] unless ($line[2] < $oldend);  #reset 'end' unless nested gene
		
	} #end foreach @sorted
} #end for each HOA

my $avglen = $genelength/$totalgenes;
printf $fh "Average gene length is %.2f bp\n", $avglen;

my $average = $totaldist/$totaldistgenes;
printf $fh "Average (pseudo) intergenic distance is %.2f bp (max is %.0f bp)\n", $average, $maxdist;

print $fh "There are $negativedist pairs with negative distance (omitted for average)\n";

print $fh "total genes: $totalgenes  ($totaldistgenes used for avg distance)\n";

my $size = keys %seqidHOA;
print $fh "total chr: $size\n";

print $fh "\n";

print $fh "\n";

##======================== output annotation 'type' data ====================##
print $fh "------------------------------------------------------------\n";
print $fh "ANNOTATION 'TYPE' DATA FROM THIS GFF FILE:\n";
foreach my $key (sort keys %annot_types){
	print $fh "$key\n";
	}

print $fh "\n";	

##============ list chromosomes and sizes from GFF  header====================##
print $fh "\n";
print $fh "------------------------------------------------------------\n";
print $fh "SEQUENCE-REGIONS WITH GENES/EXONS (bp) -- does not assess other types:\n";

my $undeffrag = 0;
my $maxundef = 0;

##this part updated v1.3 to fix problem when "##sequence-region" != seqID
foreach my $key (sort keys %headerchr){
	next if $headerchr{$key}{seen} == 0;
	if (exists $seqid{$key}){
		print $fh "$key\t$headerchr{$key}{length}\n";
	} else {
		$undeffrag++;
		$maxundef = $headerchr{$key}{length} unless $maxundef > $headerchr{$key}{length};
	}	
}

my $avg = 0; 
#this prevents a divide-by-zero error if no ##sequence-region defined
if ($undeffrag > 0){
	my $avg = $maxundef/$undeffrag;
}

printf $fh "\nThere are %.0f sequence-regions with no genes/exons (largest is %.0f bp, average %.2f bp)\n", $undeffrag, $maxundef, $avg;
print $fh "\n";	



close($fh);

__END__

some notes: need to think about ID vs Name and which we want; can this be selected at run time via modification of scrm.pl and gff.pl. This preflight for instance could be used as a screening tool to decide which to use and then the proper selection can be made.

