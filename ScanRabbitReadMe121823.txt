## ScanRabbit-- Assembly-indepdendent filtering of NGS Datasets for a well defined protein motif
##              Default settings are for Candidate Obelisk RNAs as per I.Zheludev profiles; These can be modified for any peptide profile
##              Version an08 from 12-18-2023
## ->Syntax
##  - python ScanRabbit<ver>.py  Substrate=<file,file..>  msafile="COR_ORF_1.msa" <Other parameters>
##
## Further details of the program and it's paramters are below:
##
## ->SCANNING
##   The scan is carried out in two steps, first to identify a "core" sequence that nominates an individual segment as a potential match, second uses a profile to obtain a probability score
##   We define a "Core" as a contiguous region of the target protein family that will be used to nominate invidual reads or segments for further analysis
##   and a "Block" region containing this that will be used to calculate the match score based on the profile
##   The geometry of the block we are scanning for is --------BBBCCCCCCCCCCBBBBB----- where C elements are the "Core" and B+C elements form the "Block"
##   The user provides a multiple sequence alignment (MSA) file along with
##       A value indicating where the Block starts in the MSA (Offset)
##       A value indicating where the Core starts within the block (BlockToCoreOFfset)
##       A value for the length of the Core
##       A value for the length of the block
##   For a multiple sequence alignment file where the context of the Core and block is as followd --------BBBCCCCCCCCCCBBBBB-----
##       Offset = 8, BlockToCoreOffset=3, BlockLength=18, CoreLength=10
##   All offsets are zero-based
##   Two kinds of prescan are possible, a motif based scan and a profile based scan
##   The standard motif '[RYFGHPWY][KRISHCV]RGY.[DE].[GNST][ST]' captures 461/505 instances of Obelisk ORF1 in the most conserved region (91.3%)
##   This is not a full regular expression syntax: Each position can have one or more amino acids (multiple possibilities in brackets), with periods indicating no specification
##   The complexity of the dictionary used for the second search filter comes from the number of combinations here
##   Fully unspecified positions don't add to the size, but positions with many possible amino acids do.  So a position where there is little constraint
##   Should generally be left as a '.'
##   There is no allowance here for gaps in the conserved block.
##   B sequences upstream and downstream are in the block used to calculate a score based on the multiple sequence alignment
##   B sequences should be in a region with no gaps in the alignment but can have some flexibility
##  -   Here are the parameters used to specify the scan
##  - MSAFile = <'COR_ORF_1.msa'> : this is the multiple sequence alignment file ('COR_ORF_1.msa' for the case of  Obelisk ORF1)
##  - CoreSpecification = <'[RYFGHPWY][KRISHCV]RGY.[DE].[GNST][ST]'> ## A profile that will be used as an intermediate search prefilter (Default is for the CORE or COR Obleisks)
##  - ProfilePreFilter = <False> ## This replaces the explicit pattern matching with "CoreSpecification1" with implicit use of probability from the profile
##  -   CoreLength = <10>    # Length of the prefiler corerecommend 8,9,10,12
##  -   BlockToCoreOffset1 = 0 ## How many amino acids are upstream of the core but still in block
##  -   Offset1 = <176>  ## where the conserved block starts in the  msa (zero-based and amino acid counts)
##  -   BlockLength1 = <18> ## how many amino acids in the conserved block (including core).
##  -   MinPreFilterScore=<8.0> : Minimum score in "core" sequence-based prefiltering (in log2 units)
##  -   MinMatchScore=<10.0> : Minimum match score for reporting for full block (in log2 units)
##  -   PositionAdjust = <False> ## Setting to true allows the user to skip columns in the MSA that have fewer than half of active rows filled in
##
## ->INPUT
##  - Substrate=<FilesSpecification>
##       This can be a single file, a group of files (separated by commas, a wildcard specification, or a file of filenames
##       For the file-of-filenames option, the file should be named with a ".files" file extension and have a file name at the beginning of each row
##       Other data can follow as long as there is a tab after the filename
##       (As always, please try to avoid spaces and commas in file names)
##       An SRR/DRR/ERR entry can also be used directly (e.g. SRR1234567) as a Substrate and Rabbit will try to find fast(er)q-dump to download
##  - FirstRead=<0> ## A positive integer n will ignore the first n reads on each file (default is 0: start at beginning)
##  - LastRead=<0> ## A positive integer n will stop after n reads (LastRead=0 will read until end-of-file)
##  - PairedReads=<True> ## This instructs Rabbit to look for paired R1 and R2 files (will still process single reads if there are no obvious R1/R2 pairs)
##
## ->OUTPUT
##  - ReadByReadOutput <True> ## Provides a Fastq or FastA file with individual reads
##  -   MaxFlankingLength <300> ## Maximum length flanking the match to provide in fasta/fastq 
##  - MatchByMatchOutput <False> ## Provides a MatchByMatch list of matches
##  -   MatchByMatchFullSequence <False> ## Provides a full sequence for each read in the Match output table
##  - FileByFileOutput <True> ## Provides a file-by-file summary table of matches
##  -   MatchSpeciesToReport = <5> ## Zero reports all sequences matching criteria; positive integer shows sequence for the most common n+highest scoring 5 matches (only 5 if these sets ovrerlap)
##  - LogOutput <True> ## Provides a Log File for troubleshooting and workflow archiving 
##  -   ReportGranularity=<1000000> ## How often to report progress (default every 1000000 reads)
##
## ->MULTIPROCESSING
##  - Threads = <1> ## For large numbers of files or single very large files, setting threads to the number of used cpu cores (e.g. 16) may speed up the scan dramatically
##  - Interleave = False ## Setting this to true sets up within-file multithreading for one or a few large files, you will still need to provide a number for "Threads"
##       for a single pair of files (or single file) that is gzip compressed, there is a modest return-on-investment for multiprocessing through interleave
##       It appears that the optimal value for threads here is ~4 on a MacBookAir with about a 2x sped up processing.
##
## ->Other Notes
##  - For SRR/ERR/DRR datasets not present locally, ScanRabbit will try to download the files from NCBI (linux/mac only)
##  - All features should work on Linux/Mac/Windows with Python 3.7+, but the program will run much more quickly with the compiled python tool PyPy
##  - The pypy interpreter (www.pypy.org) is available "pypy.org".
##  -    Explicit installation of PyPy is not required Instead of a dedicated Pypy installation, you can just download, unzip, and run
##  -      e.g. with the command line: <Folder name for Pypy>/bin/pypy3 ScanRabbit<version>.py Substrate=<files> <other parameters>
