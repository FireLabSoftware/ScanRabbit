#!/usr/bin/env python
## 
## ######################
## ScanRabbit-- Assembly-indepdendent filtering of NGS Datasets for a well defined protein motif
##              Default settings are for Candidate Obelisk RNAs as per I.Zheludev profiles; These can be modified for any peptide profile
##              Version cc02 from 01-15-2025
## ->Syntax
##  for searching with a multiple sequence alignment query (msa file)
##  - python ScanRabbit<ver>.py  Substrate=<file,file..>  msafile="<msa_file_name> <Other parameters>
##  for a simple query (single amino acid sequence), the command line can be
##  - python ScanRabbit<ver>.py  Substrate=<file,file..>  Query=<aa sequence or fasta file> <Other parameters>
##  for an obelisk search as carried out in Zheludev et al you can load the defaults used in that study with the command line item "obelisk"
##  - python ScanRabbit<ver>.py  Substrate=<file,file..>  obelisk
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
##   For a multiple sequence alignment file where the context of the Core and block is as follows --------BBBCCCCCCCCCCBBBBB-----
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
## ->DISTANCE METRICS 
##  - Several different approaches can be used to determine membership in the class of interest
##  -   By setting "bayesScore=True" the user can specify a bayesian approach that uses amino acid frequencies at each position
##  -       in the specified alignment.  Each aligned sequence is given one "vote" and to avoid division by zero, additional votes
##  -       added that are split equally between the 20 amino acids and stop codon (so-called regularization).
##  -     BayesRegularPedestal adds the equivalent of a specified number of reads to the multiple sequence alignment with equal representation of each amino acid
##  -     BayesRegularProportional adds additional regularization based on a proportion of the "real" sequences in the alignment
##  -  Alternatively, the scoringMatrix setting can be used to choose a matrix of amino acid relatedness values to calculate a score from.
##  -     Choices of matrix are
##  -       blosum45, blosum62, blosum80: Standard scoring matrix based on protein block alignments with increasing stringency
##  -         values are from NCBI (source is their C Toolkit Cross Reference "c/​data/​BLOSUMn")
##  -       aaEqual (counts numbers of identies to aligned reference sequences, with equal score for each aa)
##  -       baseEqual (scores identities to reference sequences using the infrequency of codons to upweight amino acids with fewer codons (e.g. M or W)
##  -  Reported score values are in units of log2, so a value of 20 should be obtained from only 1 in a million random sequences
##  -    There are many assumptions and approximations in score calculation and the scores will change substantially depending on the
##  -       choice of matrix and other parameters.
##  -    minMatchScore: The default score cutoff for reporting is 15.0, but this may need to be adjusted for specific purposes
##  -       (in general the 15.0 default is somewhat too generous, but this may not be as evident for blosum comparisons)
##  -  A second cutoff, minPrefilterScore (default=8) is used in profile-based prefiltering to set the cutoff below which
##  -    the remainder of a sequence will not be examined in detail.
##  -  In some cases, scanRabbit will report matches that fail in the core region but spuriously were picked due to third base ambiguity
##  -    Such cases may actually be 'real' matches caught fortuitously-- but they can be turned off if desired by setting
##  -     enforceCore=True
##  - 
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
## ###############
## End Help
SubstrateFiles1 = '' ##  Can be explicit lists of files, or wildcards.  Rabbit will look for R1 and R2 files that can be paired
OutFileBase1 = 'default' ## A base file name for R1 and R2 output and log file
ReportGranularity1 = 1000000 ## how often to report results
minMatchScore1 = 15.0 ## minimum score for match (in log2 units)
Threads1 = 1
Interleave1 = False
myThread1 = -1
Delimiter1 = '\r'
ProvideReads1 = True ## Setting this to False just provides a FileByFile summary (no individual reads)
MatchSpeciesToReport1 = 5 ## Zero reports all sequences meeting criteria.  Setting to a positive integer n (e.g. =5) shows sequence for only the n most common blocks 
AbbreviateSummary1 = 'default' ## Setting this to False outputs each matching sequence onto a separate line, true mashes all matches into a single semicolon-delimited list; default is True for single file runs, otherwise false
PositionAdjust1 = False ## Setting this to true allows the user to skip columns in the MSA that have fewer than half of active rows filled in
## Some additional user-specifiable values
FastQDumpProgram1 = 'fasterq-dump' ## Full path of fastq-dump/fasterq-dump program to download from SRA if needed.

## Some parameters that would change with a different motif search
msaFn1 = '' ## For Obelisk search based on Zheludeve et al, use 'COR_ORF_1.msa'
query1 = '' ## can be a fasta amino acid sequence file name (single reference) or an amino acid sequence
scoringmatrix1 = 'blosum62' ## Matrix used for direct sequence queries. choices are blosum45, blosum62, blosum80, aaEqual, and baseEqual
bayesScore1 = False ## Uses a bayesian posterior probability instead of a summed score
bayesRegularPedestal1 = 2.0 ##1.0 ##2.0    ## A regularization applied when calculating match scores from an alignment--
## A value of 2 effectively adding the equivalent of two aligned reference sequences with random amino acid composition to the profile
bayesRegularProportion1 = 0 ## 0.05   ## Similar to BayesRegularPedestal1 but adding random-amino-acid reads as a proportion of the total number of sequences in the alignment
enforceCore1 = True ## setting this to true enforces that the Core needs to match (with ambiguity but no mismatches
##   A standard multiple sequence alignment in msa format (based on fasta)
##   For 'Rabbit' to work for this you will need to identify a highly conserved non-gapped portion of the alignment which is at least 10 aa
CoreSpecification1 = '' ## for Obelisk search based on Zheludeve et al, use'[RYFGHPWY][KRISHCV]RGY.[DE].[GNST][ST]' ## This is a profile that will be used as an intermediate search prefilter.
CoreLength1 = 10    # recommend 8,9,10,12
ProfilePreFilter1 = False ## setting this to true uses an alternative profile-based filtering method
## if a core specification is not provided, a profile will be built using the provided MSA file.  Note that this will then use a core that is defined
## as starting at position Offset+BlockToCoreOffset in the alignment with length CoreLength.  Effective recovery of useful sequences may also
## involve setting MinPrefilterScore to a larger value (default is 8)
Offset1 = 0 ## for Obelisk search based on Zheludeve et al, use 176  ## where the conserved block starts in the  msa (zero-based and amino acid counts)
BlockLength1 = 0 ## for Obelisk search based on Zheludeve et al, use 18  ##how many amino acids in the conserved block (including core).
##   The entire block needs to be contained in a single read, so making this block too long may lose instances where only a part is contained
BlockToCoreOffset1 = 0 ## How many amino acids are upstream of the core but still in block
ReadsToProcess1 = 0 ## Setting this to a positive integer limits the number of reads (read pairs) processed for each file
MatchByMatchOutput1 = False ## Provides a list of matches on a line-by-line basis
SummaryOutput1 = True ## Provides a file-by-file output that describes "greatest hits" in each file
LogOutput1 = True ## Generates a log output file for debugging

## An alternative pre-filtering scheme that ignores the regular-expression-like motif and uses the multiple sequence alignment
## A Core Length is needed for this-- the relevant core section will start at position Offset1+BlockToCoreOffset1 in the multiple sequence alignment
## And move forward for CoreLength amino acids.  
MinPrefilterScore1 = 8.0 ## If a profile score (rather than pattern match) is used as a prefilter, this is the score required for analysis.
                         ## Units are log(2) so a value of 8 would mean that ~1/2^8 k-mers would meet the prefilter criteria
MaxFlankingLength1 = 300 ## This value caps the distance away from the core that is reported for any individual match
MatchByMatchWholeSequence1 = False ## This provides the entire sequence of the instances in match-by-match output  
StartRead1 = 0 ## Where to start in the substrate file (reads)
EndRead1 = 0 ## Where to stop in the substrate file (reads)
StartChar1 = 0 ## Where to start in the substrate file (bytes) ## This is for specific multithreading tasks and generally not useful
EndChar1 = 0 ## Where to stop in the substrate file (bytes) ## This is for specific multithreading tasks and generally not useful
PairedReads1 = True ## Try to Pair R1 and R2 files if possible
confirm_gzip1 = True ## This will (fairly quickly) confirm that gzip files are intact and skip files that raise an error
import gzip
from glob import glob
import os,sys,re
from time import time, strftime, localtime
import subprocess
from itertools import chain, product, zip_longest
from collections import Counter
from math import log
from random import randint


for a0 in sys.argv[1:]:
    if a0.startswith('=') or a0.endswith('='):
        print('Encountered equals sign in unexpected position in command line; syntax should be free of spaces')
        exit()
    if a0.startswith('h') or a0.startswith('-h'):
        for L1 in open(sys.argv[0],mode='rt').readlines():
            if L1.startswith('## End Help'): break
            if L1.startswith('##'):
                print(L1.strip('#').strip())
        exit()
    if a0.lower().startswith('obelisk'): ## This will load defaults for an Obelisk Motif search as done in Zheludev et al
        msaFn1 = 'COR_ORF_1.msa'
        CoreSpecification1 = '[RYFGHPWY][KRISHCV]RGY.[DE].[GNST][ST]' 
        Offset1 = 176  
        BlockLength1 = 18
        continue
    s1 = a0.split('=')[0].strip().strip("'").strip('"')
    v1 = a0.split('=')[1].strip().strip("'").strip('"')
    if s1.lower().startswith('substrate'):
        SubstrateFiles1 =  v1
    if s1.lower().startswith('fastq') or s1.lower().startswith('fasterq'):
        FastQDumpProgram1 =  v1
    if s1.lower().startswith('scoring'):
        scoringmatrix1 = v1
    if s1.lower().startswith('que'):
        query1 = v1
    if s1.lower().startswith('msa'):
        msaFn1 = v1
    if s1.lower().startswith('corespecification'):
        CoreSpecification1 = str(v1)
    elif s1.lower().startswith('report'):
        ReportGranularity1 =  int(v1)
    elif s1.lower().startswith('maxflanking'):
        MaxFlankingLength1 =  int(v1)
    elif s1.lower().startswith('threads'):
        Threads1 =  int(v1)        
    elif s1.lower().startswith('offset'):
        Offset1 =  int(v1)        
    elif s1.lower().startswith('mythread'):
        myThread1 =  int(v1)        
    elif s1.lower().startswith('startchar'):
        StartChar1 =  int(v1)        
    elif s1.lower().startswith('endchar'):
        EndChar1 =  int(v1)        
    elif s1.lower().startswith('startread'):
        StartRead1 =  int(v1)        
    elif s1.lower().startswith('endread'):
        EndRead1 =  int(v1)        
    elif s1.lower().startswith('blocklength'):
        BlockLength1 =  int(v1)        
    elif s1.lower().startswith('blocktocore'):
        BlockToCoreOffset1 =  int(v1)        
    elif s1.lower().startswith('corelength'):
        CoreLength1 =  int(v1)        
    elif s1.lower().startswith('bayesregularped'):
        bayesRegularPedestal1 =  float(v1)        
    elif s1.lower().startswith('bayesregularpro'):
        bayesRegularProportion1 = float(v1)        
    elif s1.lower().startswith('providereads') or s1.lower().startswith('readbyreadout'):
        if v1.lower()[0] in 'fn-':
            ProvideReads1 = False
        else:
            ProvideReads1 = True       
    elif s1.lower().startswith('profilepre'):
        if v1.lower()[0] in 'fn-':
            ProfilePreFilter1 = False
        else:
            ProfilePreFilter1 = True       
    elif s1.lower().startswith('confirm'):
        if v1.lower()[0] in 'fn-':
            confirm_gzip1 = False
        else:
            confirm_gzip1 = True       
    elif s1.lower().startswith('abbreviate'):
        if v1.lower()[0] in 'fn-':
            AbbreviateSummary1 = False
        else:
            AbbreviateSummary1 = True       
    elif s1.lower().startswith('pairedreads'):
        if v1.lower()[0] in 'fn-':
            PairedReads1 = False
        else:
            PairedReads1 = True
    elif s1.lower().startswith('enforcecore'):
        if v1.lower()[0] in 'fn-':
            enforceCore1 = False
        else:
            enforceCore1 = True
    elif s1.lower().startswith('interlea'):
        if v1.lower()[0] in 'fn-':
            Interleave1 = False
        else:
            Interleave1 = True       
    elif s1.lower().startswith('positionadjust'):
        if v1.lower()[0] in 'fn-':
            PositionAdjust1 = False
        else:
            PositionAdjust1 = True       
    elif s1.lower().startswith('bayesscore'):
        if v1.lower()[0] in 'fn-':
            bayesScore1 = False
        else:
            bayesScore1 = True
    elif s1.lower().startswith('matchbymatchout'):
        if v1.lower()[0] in 'fn-':
            MatchByMatchOutput1 = False
        else:
            MatchByMatchOutput1 = True
    elif s1.lower().startswith('logout'):
        if v1.lower()[0] in 'fn-':
            LogOutput1 = False
        else:
            LogOutput1 = True
    elif s1.lower().startswith('matchbymatchwhole'):
        if v1.lower()[0] in 'fn-':
            MatchByMatchWholeSequence1 = False
        else:
            MatchByMatchWholeSequence1 = True
    elif s1.lower().startswith('summaryout') or s1.lower().startswith('filebyfileout'):
        if v1.lower()[0] in 'fn-':
            SummaryOutput1 = False
        else:
            SummaryOutput1 = True
    elif s1.lower().startswith('minmatchscore'):
        minMatchScore1 =  float(v1)        
    elif s1.lower().startswith('minprefilterscore'):
        MinPrefilterScore1 =  float(v1)        
    elif s1.lower().startswith('matchspeciestoreport'):
        MatchSpeciesToReport1 =  float(v1)        
    elif s1.lower().startswith('out'):
        OutFileBase1 =  str(v1)        
t0 = time()
vnow = strftime("%m%d%y_%H%M%S",localtime())

def pypyreminder1(context):
    if not('pypy' in sys.version.lower()) and myThread1<1:
        print()
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        if context=='before':
            print('WARNING: ScanRabbit will be running in standard Python, not PyPy')
        else:
            print('WARNING: ScanRabbit has just run in standard Python, not PyPy')
        print('Rabbit will be very slow with standard Python (up to 10x slower than PyPy)')
        print('Download current native PyPy from PyPy.org for full-speed function')
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print()
pypyreminder1('before')


if Threads1>1:
    if Interleave1:
        ThreadMessage1 = 'Interleaf_'+str(myThread1)+'_of_'+str(Threads1)+'_'
    else:
        ThreadMessage1 = 'Thread_'+str(myThread1)+'_of_'+str(Threads1)+'_'        
else:
    ThreadMessage1 = ''

if OutFileBase1 == 'default':
    OutFileBase1 = 'ScanRabbit_'+vnow
if SummaryOutput1:
    SummaryFileName1 = OutFileBase1+'_Summary.tsv'
    SummaryFile1 = open(SummaryFileName1, mode='w')
if MatchByMatchOutput1:
    MatchByMatchFileName1 = OutFileBase1+'_MatchByMatch.tsv'
    MatchByMatchFile1 = open(MatchByMatchFileName1, mode='w')
if LogOutput1:
    LogFileName1 = OutFileBase1+'_LogFile.tsv'

BaseArray1 = []
for i in range(256):
    if chr(i).upper() in 'GATCN':
        BaseArray1.append(chr(i).upper())
    elif chr(i) in 'SDFHJKLQWERYIOPZXVBM*':
        BaseArray1.append('N')
    elif chr(i) in 'U':
        BaseArray1.append('T')
    else:
        BaseArray1.append('')
def FileInfo1(FileID):
    if type(FileID)==str:
        Fn1=FileID
    else:
        Fn1=FileID.name
    s11=os.stat(Fn1)
    return ','.join([Fn1,
                     '#',
                      'Path='+os.path.abspath(Fn1),
                        'Size='+str(s11[6]),
                        'Accessed='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[7])),
                        'Modified='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[8])),
                        'Created='+strftime("%m_%d_%y_at_%H_%M_%S",localtime(s11[9])),
                        'FullInfo='+str(s11)])



BaseD1 = Counter({'g':0,'a':1,'t':2,'c':3,'G':0,'A':1,'T':2,'C':3,'n':0,'N':0})
BaseL1 = [0]*256
BaseA1 = [0]*256
aBaseL1 = ['G','A','T','C']
for b1 in BaseD1:
    BaseL1[ord(b1)] = BaseD1[b1]
    BaseA1[ord(b1)] = 3-BaseD1[b1]

## for debugging-- keep
def SeqToValue1(s):
    return sum([4**i*BaseL1[ord(c)] for (i,c) in enumerate(s[::-1])])
def ValueToSeq1(v,l):
    return ''.join([aBaseL1[(v>>(2*i)) & 3] for i in range(l-1,-1,-1)])

## antisense-- returns the reverse complement of argument string 
vfilterminus=''
vASB11="AaCcNn*nNgGtT"
for i in range(256):
    if chr(i) in vASB11:
        vfilterminus+=vASB11[12-vASB11.find(chr(i))]
    else:
        vfilterminus+='N'  ## switched from filterminus+='' 5/08/21 to avoid crashes from antisense and sense sequences having different lengths
def vantisense(s):
    '''return an antisense and filtered version of any sequence)'''
    return s.translate(vfilterminus)[::-1]


def LogNote1(*notes):
    note = ThreadMessage1+' '.join(map(str,notes))
    if LogOutput1:
        LogFile=open(LogFileName1,mode='a')
        if note.strip().strip('*')=='':
            LogFile.write(note+Delimiter1)
        else:
            LogFile.write(note+'\t'+'; t='+"{0:.3f}".format(time()-t0)+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' '+Delimiter1)
        LogFile.close()
    if note.strip().strip('*')=='':
        print(note)
    else:
        print(note.split('#')[0].replace('\t',' ').replace('\r','\n').strip(',') + '; t='+"{0:.3f}".format(time()-t0))

if myThread1==-1:
    LogNote1('Running ScanRabbit with parameters:\n  Version:','\n  '.join(sys.argv),'\n  Python Version:',sys.version,'\n\n')

class BlockInstance1():
    def __init__(self,core,score,ups,dwn,uid):
        self.core = core
        self.ups = [ups,]
        self.dwn = [dwn,]
        self.score = score
        self.uid = [uid]   ## uid is (readnumber, start position) 1-based. read num negative for R2, pos negative for antisense match
    def countme(self):
        self.count = len(set([abs(x[0]) for x in self.uid]))
    def assemble(self):
        plusL = []
        for u in self.ups:
            for i,c in enumerate(u[::-1]):
                if i>=len(plusL):
                    plusL.append(Counter())
                plusL[i][c] += 1
        assN = ''
        for c in plusL:
            nextc = c.most_common(1)[0][0]
            totac = sum(c.values())
            if c[nextc]/totac>0.9:
                assN+=nextc.upper()
            elif c[nextc]/totac>0.5:
                assN+=nextc.lower()
            else:
                assN+='x'
        self.assN = assN[::-1]
        minusL = []
        for d in self.dwn:
            for i,c in enumerate(d):
                if i>=len(minusL):
                    minusL.append(Counter())
                minusL[i][c] += 1
        assC = ''
        for c in minusL:
            nextc = c.most_common(1)[0][0]
            totac = sum(c.values())
            if c[nextc]/totac>0.9:
                assC+=nextc.upper()
            elif c[nextc]/totac>0.5:
                assC+=nextc.lower()
            else:
                assC+='x'
        self.assC = assC
        self.positions = [Counter() for i in range(4)]
        for uid in self.uid:
            if uid[0]>0 and uid[1]>0:
                self.positions[0][uid[1]] += 1
            elif uid[0]>0 and uid[1]<0:
                self.positions[1][uid[1]] += 1
            elif uid[0]<0 and uid[1]>0:
                self.positions[2][uid[1]] += 1
            elif uid[0]<0 and uid[1]<0:
                self.positions[3][uid[1]] += 1
        self.contexts = '['
        if self.positions[0]:
            self.contexts += 'R1p:'+str(sum(self.positions[0].values()))+'('+str(len(self.positions[0]))+'),'
        if self.positions[1]:
            self.contexts += 'R1m:'+str(sum(self.positions[1].values()))+'('+str(len(self.positions[1]))+'),'
        if self.positions[2]:
            self.contexts += 'R2p:'+str(sum(self.positions[2].values()))+'('+str(len(self.positions[2]))+'),'
        if self.positions[3]:
            self.contexts += 'R2m:'+str(sum(self.positions[3].values()))+'('+str(len(self.positions[3]))+'),'
        self.contexts = self.contexts[:-1]+']'
        
    def add(self,ups,dwn,uid):
        self.ups.append(ups)
        self.dwn.append(dwn)
        self.uid.append(uid)
    def report(self):
        r = self.assN+'_'+self.core+'_'+self.assC
        r += '__count='+str(self.count)
        r += '__score='+'%.3f'%self.score
        r += self.contexts
        return r
def find1(n):
    ''' find a file on this machine'''
    o = []
    for r,d,f in os.walk('/'): 
        if n in f:
            o.append(os.path.join(r, n))
    return o
 
def findfastqdump(candidate):
    ver = 0
    try:
        ver = subprocess.check_output([candidate,'-V'])
        return candidate
    except:
        pass
    if os.path.isfile('~/FastQDumpLocation.txt'):
        candidate = open('~/FastQDumpLocation.txt', mode='rt').read()
        try:
            ver = subprocess.check_output([candidate,'-V'])
            return candidate
        except:
            pass
    LogNote1('Looking for fasterq-dump-- if this fails, provide a location in the command line (fastqdump=<path>)')
    LogNote1('or reinstall and allow execution (chmod +X <path>)')
    LogNote1('Note finding fast(er)q-dump can take some real time, get a cup of coffee or tea')
    NewCands = find1('fasterq-dump')
    NewItems = []
    for candidate in NewCands:
        try:
            ver = subprocess.check_output([candidate,'-V'])
            NewItems.append([versioner1(ver),candidate])
        except:
            pass
    if NewItems:
        candidate = sorted(NewItems, reverse=True)[0][1]
        try:
            open('~/FastQDumpLocation.txt', mode='w').write(candidate)
        except:
            pass
        return candidate
    LogNote1('Unable to find fastq-dump.  If you intend to use thisRecommend that you reinstall this from ncbi or download fastq/fasta files directly')
    return '' 

if SubstrateFiles1.endswith('.files') and os.path.isfile(SubstrateFiles1):
    TempFileList1 = []
    for L0 in open(SubstrateFiles1,mode='rt'):
        TempFileList1.append(L0.split('\t')[0].strip())
    SubstrateFiles1 = FileList1
SubstrateFileList1 = []
for s11 in SubstrateFiles1.split(','):
    if not(s11): continue
    if os.path.isfile(s11):
        SubstrateFileList1.append(s11)
    elif list(glob(s11,recursive=True)):
        SubstrateFileList1.extend(sorted(list(glob(s11,recursive=True))))
    if not(SubstrateFileList1):
        SubstrateFileList1.append(s11)
if not(SubstrateFileList1):   
    LogNote1('**************************************************')
    LogNote1('MAJOR COMMAND LINE ERROR: No input substrate files found; Check your directory and filename specifications')
    LogNote1('**************************************************')
    exit()
SubstrateFileD1 = dict.fromkeys(SubstrateFileList1)

for Fn1 in list(SubstrateFileD1.keys()):
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)):
        if not(os.path.isfile(Fn1+'_1.fastq')):
            LogNote1(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
            FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,Fn1])
            LogNote1("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
        del(SubstrateFileD1[Fn1])
        if os.path.isfile(Fn1+'_1.fastq'):
            SubstrateFileD1[Fn1+'_1.fastq'] = 0
        else:
            LogNote1('Failed to download',Fn1,'- Will continue for now though')
        if os.path.isfile(Fn1+'_2.fastq'):
            SubstrateFileD1[Fn1+'_2.fastq'] = 0
SubstrateFileList1 = []
SubstrateFileList2 = []
UnitaryFileList1 = sorted(SubstrateFileD1.keys())
alreadygotthatone1 = False
AnyR2File1 = False
if PairedReads1:
    for f1,sfl1 in enumerate(UnitaryFileList1):
        if "Rabbit_Finds" in sfl1: continue
        if alreadygotthatone1:
            alreadygotthatone1=False
            continue
        if f1==len(UnitaryFileList1)-1:
            SubstrateFileList1.append(sfl1)
            SubstrateFileList2.append('')
        else:
            sfl2 = UnitaryFileList1[f1+1]
            if (sfl2.replace('_2.','_1.')==sfl1 or
                sfl2.replace('_R2.','_R1.')==sfl1 or
                sfl2.replace('_2','_1')==sfl1 or
                sfl2.replace('_R2','_R1')==sfl1 or
                sfl2.replace('R2','R1')==sfl1 or
                sfl2.replace('_2.','_1.',1)==sfl1 or
                sfl2.replace('_R2.','_R1.',1)==sfl1 or
                sfl2.replace('_2','_1',1)==sfl1 or
                sfl2.replace('_R2','_R1',1)==sfl1 or
                sfl2.replace('R2','R1',1)==sfl1 or
                sfl2[::-1].replace('.2_','.1_',1)==sfl1[::-1] or
                sfl2[::-1].replace('.2R_','.1R_',1)==sfl1[::-1] or
                sfl2[::-1].replace('2_','_1_',1)==sfl1[::-1] or
                sfl2[::-1].replace('2R_','1R_',1)==sfl1[::-1] or
                sfl2[::-1].replace('2R','1R',1)==sfl1[::-1]):
                SubstrateFileList1.append(sfl1)
                SubstrateFileList2.append(sfl2)
                alreadygotthatone1 = True
                AnyR2File1 = True
            else:
                SubstrateFileList1.append(sfl1)
                SubstrateFileList2.append('')
                alreadygotthatone1 = False
else:
    SubstrateFileList1 = UnitaryFileList1
    SubstrateFileList2 = ['',]*len(UnitaryFileList1)
if Threads1>1 and not(Interleave1):
    if Threads1>len(SubstrateFileList1): Threads1 = len(SubstrateFileList1)
    SubstrateFileList1 = SubstrateFileList1[myThread1-1::Threads1]            
    SubstrateFileList2 = SubstrateFileList2[myThread1-1::Threads1]            

TaskHeader1 =  '<!--Rabbit_Task_Header: '+OutFileBase1+'-->'+Delimiter1
TaskHeader1 += '<!--Command_Line: '+' '.join(sys.argv)+'-->'+Delimiter1
TaskHeader1 += '<!--PythonVersion: '+','.join(sys.version.splitlines())+'-->'+Delimiter1
TaskHeader1 += '<!--Rabbit_Version: '+FileInfo1(sys.argv[0])+'-->'+Delimiter1
TaskHeader1 += '<!--RunTime: '+vnow+'-->'+Delimiter1
TaskHeader1 += '<!--RunDirectory: '+os.getcwd()+'-->'+Delimiter1
TaskHeader1+=Delimiter1
AbbrevHeader1 = ''.join(TaskHeader1.splitlines()[:-1])+'<!--RabbitTableHeader-->'  ##ending with ':RabbitTableHeader' identifies a line as a row of table headers
def HeaderTranspose1(hT2):
    hT0 = '<!--\tOutput_Key\t\t-->'+Delimiter1
    hT0 += '<!--\tColumnNumber\tColumnHeader\t-->'+Delimiter1
    for iT2,nT2 in enumerate(hT2):
        nT2 = nT2.strip()
        if not(nT2.startswith('<!')):
            hT0 += '<!--\t'+str(iT2)+'\t'+nT2+'\t-->'+Delimiter1
    return hT0+'<!--\t\t\t-->'+Delimiter1

if AnyR2File1:
    NameA1 = ('R1File','R2File','directory','MatchCount','MaxScore','TotalReads','LastL1_Len','LastL2_Len','TopMatchList')
    NameA3 = ('R1File','R2File','directory','SequenceName','Strand','Read','Position','CoreCodingNA','UpstreamPep','CorePep','DownstreamPep','Score')
else:
    NameA1 = ('R1File','directory','MatchCount','MaxScore','TotalReads','LastL1_Len','TopMatchList')
    NameA3 = ('R1File','directory','SequenceName','Strand','Read','Position','CoreCodingNA','UpstreamPep','CorePep','DownstreamPep','Score')
    
Header1= '\t'.join(NameA1)+'\t'+AbbrevHeader1+'\t '+Delimiter1
if MatchByMatchWholeSequence1:
    NameA3= NameA3 + ('FullR1Sequence',)
    if AnyR2File1:
        NameA3= NameA3 + ('FullR2Sequence',)
Header3= '\t'.join(NameA3)+'\t'+AbbrevHeader1+'\t '+Delimiter1

if SummaryOutput1:
    SummaryFileName1 = OutFileBase1+'_Summary.tsv'
    SummaryFile1 = open(SummaryFileName1, mode='w')
if MatchByMatchOutput1:
    MatchByMatchFileName1 = OutFileBase1+'_MatchByMatch.tsv'
    MatchByMatchFile1 = open(MatchByMatchFileName1, mode='w')
if LogOutput1:
    LogFileName1 = OutFileBase1+'_LogFile.tsv'
if myThread1==-1:
    if SummaryOutput1:
        SummaryFile1.write(TaskHeader1)
        SummaryFile1.write(HeaderTranspose1(NameA1))
        SummaryFile1.write(Header1)
    if MatchByMatchOutput1:
        MatchByMatchFile1.write(TaskHeader1)
        MatchByMatchFile1.write(HeaderTranspose1(NameA3))
        MatchByMatchFile1.write(Header3)
        

if Threads1>1 and myThread1==-1:
    AllThreads1 = []
    ThreadUID1 = randint(100000000000,999999999999)
    for thrn1 in range(Threads1):
        AllThreads1.append(subprocess.Popen([sys.executable,]+sys.argv+['myThread='+str(thrn1+1),'out=Temp'+str(ThreadUID1)+vnow+str(10000+thrn1)]))
    for thisThread1 in AllThreads1:
        thisThread1.wait()
    if LogOutput1:
        LogFile1 = open(LogFileName1, mode='w')
    for thrn1 in range(Threads1):
        tn1 = 'Temp'+str(ThreadUID1)+vnow+str(10000+thrn1)
        if SummaryOutput1:
            sum1 = tn1 + '_Summary.tsv'
            Sum1 = open(sum1, mode='rt')
            for L1 in Sum1:
                SummaryFile1.write(L1)
            Sum1.close()
            os.remove(sum1)
        if MatchByMatchOutput1:
            hit1 = tn1 + '_MatchByMatch.tsv'
            Hit1 = open(hit1, mode='rt')
            for L1 in Hit1:
                MatchByMatchFile1.write(L1)
            Hit1.close()
            os.remove(hit1)
        if LogOutput1:
            log1 = tn1 + '_LogFile.tsv'
            Log1 = open(log1, mode='rt')
            for L1 in Log1:
                LogFile1.write(L1)
            Log1.close()
            os.remove(log1)
    if SummaryOutput1:
        SummaryFile1.close()
    if MatchByMatchOutput1:
        MatchByMatchFile1.close()
    if LogOutput1:
        LogFile1.close()
    if ProvideReads1:
        for thrn1 in range(Threads1):
            tn1 = 'Temp'+str(ThreadUID1)+vnow+str(10000+thrn1)
            reads1 = tn1+'_R1.fastq'
            AllReadsR1 = OutFileBase1+'_R1.fastq'
            if not(os.path.isfile(reads1)):
                reads1 = tn1+'_R1.fasta'
                AllReadsR1 = OutFileBase1+'_R1.fasta'
            Reads1 = open(reads1, mode='rt')
            OutFile1 = open(AllReadsR1, mode='a')
            for L1 in Reads1:
                OutFile1.write(L1)
            OutFile1.close()
            Reads1.close()
            os.remove(reads1)
            reads2 = tn1+'_R2.fastq'
            AllReadsR2 = OutFileBase1+'_R2.fastq'
            if not(os.path.isfile(reads2)):
                reads2 = tn1+'_R2.fasta'
                AllReadsR2 = OutFileBase1+'_R2.fastq'
            if os.path.isfile(reads2):
                Reads2 = open(reads2, mode='rt')
                OutFile2 = open(AllReadsR2, mode='a')
                for L2 in Reads2:
                    OutFile2.write(L2)
                Reads2.close()
                os.remove(reads2)
                OutFile2.close()
    exit()


## Three different reward/penalty matrices:
aalist1 = 'A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *'.split()
aaset1 = set('A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *'.split())
blosum1 = Counter()
blosum62M = '''   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -1 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1 -2  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  4 -3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4 -3  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -1 -3 -1 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0 -2  4 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1 -3  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0 -3  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3  3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4  3 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0 -3  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3  2 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3  0 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -3 -1 -1 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0 -2  0 -1 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -1 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -2 -2 -1 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -1 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3  2 -2 -1 -4
B -2 -1  4  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4 -3  0 -1 -4
J -1 -2 -3 -3 -1 -2 -3 -4 -3  3  3 -3  2  0 -3 -2 -1 -2 -1  2 -3  3 -3 -1 -4
Z -1  0  0  1 -3  4  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -2 -2 -2  0 -3  4 -1 -4
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1'''
blosum80M = '''   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
A  5 -2 -2 -2 -1 -1 -1  0 -2 -2 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -2 -1 -1 -6
R -2  6 -1 -2 -4  1 -1 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -4 -3 -3 -1 -3  0 -1 -6
N -2 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4  5 -4  0 -1 -6
D -2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -6 -4 -4  5 -5  1 -1 -6
C -1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1 -4 -2 -4 -1 -6
Q -1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3  0 -3  4 -1 -6
E -1 -1 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3  1 -4  5 -1 -6
G  0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4 -1 -5 -3 -1 -6
H -2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  2 -4 -1 -4  0 -1 -6
I -2 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3 -4  3 -4 -1 -6
L -2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -3 -3 -2 -2 -2  1 -4  3 -3 -1 -6
K -1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3 -1 -3  1 -1 -6
M -1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1 -3  2 -1 -1 -6
F -3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1 -4  0 -4 -1 -6
P -1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -3 -1 -3 -4  8 -1 -2 -5 -4 -3 -2 -4 -2 -1 -6
S  1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2  0 -3  0 -1 -6
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0 -1 -1 -1 -1 -6
W -3 -4 -4 -6 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -5 -4 -4 11  2 -3 -5 -3 -3 -1 -6
Y -2 -3 -3 -4 -3 -2 -3 -4  2 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2 -3 -2 -3 -1 -6
V  0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  4 -4  2 -3 -1 -6
B -2 -1  5  5 -4  0  1 -1 -1 -4 -4 -1 -3 -4 -2  0 -1 -5 -3 -4  5 -4  0 -1 -6
J -2 -3 -4 -5 -2 -3 -4 -5 -4  3  3 -3  2  0 -4 -3 -1 -3 -2  2 -4  3 -3 -1 -6
Z -1  0  0  1 -4  4  5 -3  0 -4 -3  1 -1 -4 -2  0 -1 -3 -3 -3  0 -3  5 -1 -6
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -6
* -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1'''
blosum45M = '''   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -2 -2  0 -1 -1 -1 -1 -5
R -2  7  0 -1 -3  1  0 -2  0 -3 -2  3 -1 -2 -2 -1 -1 -2 -1 -2 -1 -3  1 -1 -5
N -1  0  6  2 -2  0  0  0  1 -2 -3  0 -2 -2 -2  1  0 -4 -2 -3  5 -3  0 -1 -5
D -2 -1  2  7 -3  0  2 -1  0 -4 -3  0 -3 -4 -1  0 -1 -4 -2 -3  6 -3  1 -1 -5
C -1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -2 -2 -3 -1 -5
Q -1  1  0  0 -3  6  2 -2  1 -2 -2  1  0 -4 -1  0 -1 -2 -1 -3  0 -2  4 -1 -5
E -1  0  0  2 -3  2  6 -2  0 -3 -2  1 -2 -3  0  0 -1 -3 -2 -3  1 -3  5 -1 -5
G  0 -2  0 -1 -3 -2 -2  7 -2 -4 -3 -2 -2 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -5
H -2  0  1  0 -3  1  0 -2 10 -3 -2 -1  0 -2 -2 -1 -2 -3  2 -3  0 -2  0 -1 -5
I -1 -3 -2 -4 -3 -2 -3 -4 -3  5  2 -3  2  0 -2 -2 -1 -2  0  3 -3  4 -3 -1 -5
L -1 -2 -3 -3 -2 -2 -2 -3 -2  2  5 -3  2  1 -3 -3 -1 -2  0  1 -3  4 -2 -1 -5
K -1  3  0  0 -3  1  1 -2 -1 -3 -3  5 -1 -3 -1 -1 -1 -2 -1 -2  0 -3  1 -1 -5
M -1 -1 -2 -3 -2  0 -2 -2  0  2  2 -1  6  0 -2 -2 -1 -2  0  1 -2  2 -1 -1 -5
F -2 -2 -2 -4 -2 -4 -3 -3 -2  0  1 -3  0  8 -3 -2 -1  1  3  0 -3  1 -3 -1 -5
P -1 -2 -2 -1 -4 -1  0 -2 -2 -2 -3 -1 -2 -3  9 -1 -1 -3 -3 -3 -2 -3 -1 -1 -5
S  1 -1  1  0 -1  0  0  0 -1 -2 -3 -1 -2 -2 -1  4  2 -4 -2 -1  0 -2  0 -1 -5
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1  2  5 -3 -1  0  0 -1 -1 -1 -5
W -2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2  1 -3 -4 -3 15  3 -3 -4 -2 -2 -1 -5
Y -2 -1 -2 -2 -3 -1 -2 -3  2  0  0 -1  0  3 -3 -2 -1  3  8 -1 -2  0 -2 -1 -5
V  0 -2 -3 -3 -1 -3 -3 -3 -3  3  1 -2  1  0 -3 -1  0 -3 -1  5 -3  2 -3 -1 -5
B -1 -1  5  6 -2  0  1 -1  0 -3 -3  0 -2 -3 -2  0  0 -4 -2 -3  5 -3  1 -1 -5
J -1 -3 -3 -3 -2 -2 -3 -4 -2  4  4 -3  2  1 -3 -2 -1 -2  0  2 -3  4 -2 -1 -5
Z -1  1  0  1 -3  4  5 -2  0 -3 -2  1 -1 -3 -1  0 -1 -2 -2 -3  1 -2  5 -1 -5
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -5
* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1'''
if '62' in scoringmatrix1:
    blosumM1 = blosum62M
if '45' in scoringmatrix1:
    blosumM1 = blosum45M
if '80' in scoringmatrix1:
    blosumM1 = blosum80M

for i1,L1 in enumerate(blosumM1.splitlines()):
    if i1>0:
        L2 = map(int,L1.strip().split()[1:])
        myAA = L1.split()[0]
        for aa,val in zip(aalist1,L2):
            blosum1[(myAA,aa)] = val
aaEqual1 = Counter()
for a1 in aalist1:
    for a2 in aalist1:
        if a1==a2:
            aaEqual1[(a1,a2)] = 2*log(20,2)
        else:
            aaEqual1[(a1,a2)] = -1
baseEqual1 = Counter() ## half bit rewards for any give amino acid matching at random (2*log2(1/random-probability of aa))
baseEqual0 = Counter() ## raw probability of random match 
## reward penalty matrix based on equal base composition, presumed random codon choice, and requiring aa identity
aa2='GGGGEEDDVVVVAAAARRSSKKNNMIIITTTTW*CC**YYLLFFSSSSRRRRQQHHLLLLPPPP'
for a1 in aalist1:
    for a2 in aalist1:
        if a1 in aa2 and a1==a2:
            baseEqual1[(a1,a2)] = 2*log(64/aa2.count(a1),2)
            baseEqual0[a1] = aa2.count(a1)/64
        else:
            baseEqual1[(a1,a2)] = -1
SD1 = {}
s1 = []
profileCounts1 = []  ## Index is position in block, value is a counter that maps amino acid to occurences in reference alignment
profileWeights1 = [] ## Index is position in block, value is a counter that maps amino acid to weight (in 1/100 bit units)
if scoringmatrix1.lower().startswith('blo'):
    scoringmatrix1 = blosum1
elif scoringmatrix1.lower().startswith('aa'):
    scoringmatrix1 = aaEqual1
elif scoringmatrix1.lower().startswith('ba'):
    scoringmatrix1 = baseEqual1

if query1:
    if os.path.isfile(query1):
        query1 = ''.join([x.strip() for x in open(query1,mode='rt').readlines() if x[0]!='>']).upper()
    for i1,c1 in enumerate(query1):
        profileCounts1.append(Counter())
        profileWeights1.append(Counter())
        profileCounts1[i1][c1] += 1
        for c2 in 'ACDEFGHIKLMNPQRSTVWY*':
            profileWeights1[i1][c2] = int(50*scoringmatrix1[(c1,c2)])  ## scoring matrix is in half-bit units, profileWeights in 100th bit units
    if BlockLength1 == 0:
        BlockLength1 = len(query1)
else:
    msaF1 = chain(open(msaFn1,mode='rt').readlines(),['>'])
    for L0 in msaF1:
        if L0[0]=='>':
            if s1:
                SD1[n1] = ''.join(s1).upper()
                sL1 = len(SD1[n1]) ## length of first entry in file should be representative
            n1 = L0.strip()[1:]
            s1 = []
        else:
            s1.append(L0.strip())
    for i0,s1 in enumerate(SD1.values()):
        for i1,c1 in enumerate(s1):
            if i1>=len(profileCounts1):
                profileCounts1.append(Counter())
                profileWeights1.append(Counter())
            c1 = c1.upper()
            if c1 in aaset1:
                profileCounts1[i1][c1] += 1
    if not(profileCounts1[-1]):
        profileCounts1 = profileCounts1[:-1]
        profileWeights1 = profileWeights1[:-1]
    if BlockLength1 == 0:
        BlockLength1 = sL1

    for i1,pC1 in enumerate(profileCounts1):
        posCount1 = max(1,sum(pC1.values()))
        for refaa1 in 'ACDEFGHIKLMNPQRSTVWY*':
            if bayesScore1:
                reg1 = bayesRegularPedestal1+bayesRegularProportion1*posCount1
                regularizedRatioYes1 = (profileCounts1[i1][refaa1]+reg1/21)/(posCount1+reg1)
                regularizedRatioNo1 = baseEqual0[refaa1]
                profileWeights1[i1][refaa1] = 100*log(regularizedRatioYes1/regularizedRatioNo1,2)
            else:
                for samaa1 in 'ACDEFGHIKLMNPQRSTVWY*':
                    profileWeights1[i1][samaa1] += 50*scoringmatrix1[(refaa1,samaa1)]*(profileCounts1[i1][refaa1]/posCount1)

profileT1 = [sum(profileCounts1[i].values()) for i in range(len(profileCounts1))]
maxprofileT1 = max(profileT1[Offset1:Offset1+BlockLength1])

if PositionAdjust1:
    positionConvert1 = [] ## adjust position to account for alignment columns with no entry in most instances
    netoffset1 = 0
    for x1 in range(len(profileT1)):
        if x1<Offset1:
            positionConvert1.append(x1)
        else:
            if (profileT1[x1]<maxprofileT1/2) and max(profileT1[x1:])>2*profileT1[x1]:
                if len(positionConvert1)<Offset1+BlockLength1+1:
                    LogNote1('Adjusting alignment at downstream position '+str(x1)+' due to too few instances in alignment')
            else:
                positionConvert1.append(x1)
    profileCounts1 = [profileCounts1[x1] for x1 in positionConvert1]
    profileT1 = [profileCounts1[x1] for x1 in positionConvert1]
    profileWeights1 = [profileWeights1[x1] for x1 in positionConvert1]
profileWeightList1 = []
for p1 in range(BlockLength1):
    for v1 in range(64):
        profileWeightList1.append(int(round(profileWeights1[Offset1+p1][aa2[v1]])))


pList1 = []
for i11 in range(BlockLength1):
    pList1.append((i11*3,i11*3+1,i11*3+2,64*i11))
def FullValue1(s):
    v = 0
    for i11,i22,i33,o11 in pList1:
        v += profileWeightList1[o11+16*BaseL1[ord(s[i11])]+4*BaseL1[ord(s[i22])]+BaseL1[ord(s[i33])]]
    return v
pListA1 = []
for i11 in range(BlockLength1):
    pListA1.append((3*BlockLength1-1-i11*3,3*BlockLength1-1-i11*3-1,3*BlockLength1-1-i11*3-2,64*i11))
def FullValueA1(s):
    v = 0
    for i11,i22,i33,o11 in pListA1:
        v += profileWeightList1[o11+16*BaseA1[ord(s[i11])]+4*BaseA1[ord(s[i22])]+BaseA1[ord(s[i33])]]
    return v
# profileWeightListBinarized1 = []
# for p1 in range(BlockLength1-1,-1,-1):
#     for v1 in range(64):
#         profileWeightListBinarized1.append(profileWeights1[Offset1+p1][aa2[v1]])

# def FullValueBinarized1(n):
#     n0 = n
#     myV = 0
#     for i in range(BlockLength1):
#         myV += profileWeightListBinarized1[64*i+(n0&63)]
#         n0 >>= 6
#     return myV
##General Translation Code
aa1 = 'ACDEFGHIKLMNPQRSTVWY*'
aa2='GGGGEEDDVVVVAAAARRSSKKNNMIIITTTTW*CC**YYLLFFSSSSRRRRQQHHLLLLPPPP'
aaL2=list(aa2)
def translate1(s):
    return ''.join(aaL2[16*BaseL1[ord(s[i])]+4*BaseL1[ord(s[i+1])]+BaseL1[ord(s[i+2])]] for i in range(0,len(s)-2,3))
if CoreSpecification1 =='':
    ProfilePreFilter1 = True

if not(ProfilePreFilter1):
    m1 = []
    i1 = 0
    brac1 = ''
    CoreRE1 = re.compile(CoreSpecification1.replace('*','\\*'))
    while i1<len(CoreSpecification1):
        if CoreSpecification1[i1] in aa1 and not(brac1):
            m1.append(CoreSpecification1[i1])
        elif CoreSpecification1[i1]=='.':
            m1.append(CoreSpecification1[i1])
        elif CoreSpecification1[i1]==']':
            m1.append(''.join(brac1.strip()))
            brac1=''
        elif CoreSpecification1[i1]=='[':
            brac1 = ' '
        else:
            brac1 += CoreSpecification1[i1]
        i1+=1
    klen1 = 3*len(m1)
    m2 = []
    mask0 = 0
    ## create a binary mask that has 1 values at the first two positions of each triplet (e.g. b110110110110110 for a 5-mer protein core)
    ## create a dictionary of values for the binary representation that match the specification in CoreSpecification1
    for i1,v1 in enumerate(m1):
        if v1=='.':continue
        mask0 += (4*(64**(len(m1)-1-i1)))*15
        m2.append(set())
        for j1,c2 in enumerate(aa2):
            if c2 in v1:
                v2=j1>>2
                m2[-1].add(v2*(4*(64**(len(m1)-1-i1))))
    m2 = [list(m22) for m22 in m2]
    qs1 = set([sum(u1) for u1 in product(*m2)])  ## dictionary of kmer integer values matching the simple consensus
    ## create a binary "mask0" which is set to 1 for each bit that is the same in all binary representations that match the core specification
    ## create a binary "target0" which has one for all bit positions where all of the bits in qs1 "matches" are 1 (otherwise zeros)
    premask0 = 0
    pretarget0 = 0
    for i1 in range(2*klen1):
        e1 = 2**i1
        if e1&mask0:
            allyes1 = all([e1&qs11 for qs11 in qs1])
            allno1  = all([not(e1&qs11) for qs11 in qs1])
            if allyes1:
                premask0+=e1
                pretarget0+=e1
            elif allno1:
                premask0+=e1
    def PreFilter1(v):
        return (v&premask0==pretarget0) and (v&mask0 in qs1)

else:
    CoreRE1 = re.compile('.*')
    klen1 = 3*CoreLength1## Alternative Prefilter makes use of the Multiple sequence alignment for 8-10 amino acid core region
    if CoreLength1 in (4,6,8,10,14,16):
        SubCoreLength1 = 2  # Should be a divisor of CoreLength and higher values have a greater overhead (memory and start time)
                    # but may run faster.  Realistically should be 3 or less.  
    elif CoreLength1 in (3,9,12,15,18):
        SubCoreLength1 = 3        
    SuperProfile1 = []
    wrange1 = CoreLength1//SubCoreLength1
    for w1 in range(wrange1):
        SuperProfile1.append([])
        for x1 in range(64**SubCoreLength1):
            SuperProfile1[-1].append(0.0)
            for y1 in range(SubCoreLength1):
                aa00 = 63&(x1>>(6*(SubCoreLength1-y1-1)))
                SuperProfile1[-1][-1] += profileWeights1[Offset1+w1*SubCoreLength1+y1+BlockToCoreOffset1][aa2[aa00]]
    Mask4 = 64**SubCoreLength1-1
    yrange1 = [(6*SubCoreLength1*(CoreLength1//SubCoreLength1-w1-1)) for w1 in range(wrange1)]
    wwrange = list(range(wrange1))
    ystart1 = yrange1[0]
    ydel1 = 6*SubCoreLength1

    if CoreLength1==16:
        A0,A1,A2,A3,A4,A5,A6,A7 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A7[z1&4095]+A6[(z1>>12)&4095]+A5[(z1>>24)&4095]+A4[(z1>>36)&4095]+A3[(z1>>48)&4095]+A2[(z1>>60)&4095]+A1[(z1>>72)&4095]+A0[(z1>>84)&4095]
    elif CoreLength1==14:
        A0,A1,A2,A3,A4,A5,A6 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A6[z1&4095]+A5[(z1>>12)&4095]+A4[(z1>>24)&4095]+A3[(z1>>36)&4095]+A2[(z1>>48)&4095]+A1[(z1>>60)&4095]+A0[(z1>>72)&4095]
    elif CoreLength1==10:
        A0,A1,A2,A3,A4 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A4[z1&4095]+A3[(z1>>12)&4095]+A2[(z1>>24)&4095]+A1[(z1>>36)&4095]+A0[(z1>>48)&4095]
    elif CoreLength1==8:
        A0,A1,A2,A3 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A3[z1&4095]+A2[(z1>>12)&4095]+A1[(z1>>24)&4095]+A0[(z1>>36)&4095]
    elif CoreLength1==6:
        A0,A1,A2 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A2[z1&4095]+A1[(z1>>12)&4095]+A0[(z1>>24)&4095]
    elif CoreLength1==4:
        A0,A1 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A1[z1&4095]+A0[(z1>>12)&4095]
    elif CoreLength1==15:
        A0,A1,A2,A3,A4 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A4[z1&262143]+A3[(z1>>18)&262143]+A2[(z1>>36)&262143]+A1[(z1>>54)&262143]+A0[(z1>>72)&262143]
    elif CoreLength1==12:
        A0,A1,A2,A3 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A3[z1&262143]+A2[(z1>>18)&262143]+A1[(z1>>36)&262143]+A0[(z1>>54)&262143]
    elif CoreLength1==9:
        A0,A1,A2 = SuperProfile1
        def PreFilter1(z1):
            return MinPrefilterScore100<A2[z1&262143]+A1[(z1>>18)&262143]+A0[(z1>>36)&262143]
    elif CoreLength1==3:
        A0 = SuperProfile1[0]
        def PreFilter1(z1):
            return MinPrefilterScore100<A0[z1&262143]
    else:
        def PreFilter1(z1):
            v0 = 0.0
            for w1,y1 in enumerate(yrange1):
                v0 += SuperProfile1[w1][Mask4&(z1>>y1)]
            return MinPrefilterScore100<v0
        

minMatchScore100 = 100*minMatchScore1
MinPrefilterScore100 = 100*MinPrefilterScore1
BlockLength3 = 3*BlockLength1
BlockStartToCoreStart3 = 3*BlockToCoreOffset1
CoreEndToBlockEnd3 = 3*(BlockLength1-BlockToCoreOffset1)-klen1
StartCheckingSense1 = klen1+BlockStartToCoreStart3-1
StartCheckingAnti1 = klen1+CoreEndToBlockEnd3-1
mask1 = 4**(klen1-1)-1
ksam1 = 2*(klen1-1)

NR1 = 0 ## kept reads for run
TR1 = 0 ## total reads for run


if (SubstrateFileList1[0].lower().endswith('.fasta') or
      SubstrateFileList1[0].lower().endswith('.fa') or
      SubstrateFileList1[0].lower().endswith('.fna') or
      SubstrateFileList1[0].lower().endswith('.fa.gz') or
      SubstrateFileList1[0].lower().endswith('.fasta.gz') or
      SubstrateFileList1[0].lower().endswith('.fna.gz')):
    OutFileName1 = OutFileBase1+'_R1.fasta'
    OutFileName2 = OutFileBase1+'_R2.fasta'
else:
    OutFileName1 = OutFileBase1+'_R1.fastq'
    OutFileName2 = OutFileBase1+'_R2.fastq'
if ProvideReads1:
    OutFile1 = open(OutFileName1, mode='w')
    if AnyR2File1:
        OutFile2 = open(OutFileName2, mode='w')


if AbbreviateSummary1=='default':
    if len(SubstrateFileList1)==1 and Threads1==1:
        AbbreviateSummary1 = False
    else:
        AbbreviateSummary1 = True
for fn1,(F1,F2) in enumerate(zip(SubstrateFileList1,SubstrateFileList2)):
    F1n = os.path.basename(F1)
    F2n = os.path.basename(F2)
    LogNote1('ScanRabbit starting substrate file(s):',F1n,(F2n)*bool(F2n),'; file',fn1+1,'/',len(SubstrateFileList1))
    try:
        nr1 = 0 ## kept reads for file
        tr1 = 0  ## total processed reads for file
        tr0 = 0 ## total encountered reads for the file
        DataCycle1 = 2
        if (not(F1n.lower().endswith('fasta')) and
            not(F1n.lower().endswith('fa')) and
            not(F1n.lower().endswith('fasta.gz')) and
            not(F1n.lower().endswith('fa.gz')) and
            not(F1n.lower().endswith('fna.gz')) and
            not(F1n.lower().endswith('fna')) and
            'fastq' in F1n.lower().split('.',1)[-1]):
            DataCycle1 = 4
        if F1.lower().endswith('.gz'):
                F1o = gzip.open(F1, mode='rt')
                if F2n:
                    F2o = gzip.open(F2, mode='rt')
                else:
                    F2o = iter(())
        else:
            F1o = open(F1, mode='rt')
            if F2n:
                F2o = open(F2, mode='rt')
            else:
                F2o = iter(())
        BlockD1 = {}
        MaxScore1 = 0.0
        accu1 = []
        accu2 = []
        KeepRead1 = ''
        KeepRead3 = '' ## default is to toss
        if StartChar1:
            F1o.seek(StartChart1-1)
            if DataCycle1==4:
                for L1 in F1o:
                    if L1[0]=='@':
                        L01 = next(F1o)
                        if L01 == '@':
                            F1o.seek(-(len(L01),1))
                        else:
                            F1o.seek(-(len(L01)-len(L1),1))
                        break
            else:
                for L1 in F1o:
                    if L1[0] == '>':
                        F1o.seek(-(len(L1),1))
                        break
        NowChar1 = F1o.tell()
        for i1,(L1,L2) in enumerate(chain(zip_longest(F1o,F2o),[('>@','>@'),])):
            if not(L1): L1=' '
            if not(L2): L2=' '
            if DataCycle1==2 and i1%DataCycle1==0 and L1[0]!='>':
                DataCycle1=0
            if (DataCycle1 and i1%DataCycle1==0) or (DataCycle1==0 and L1[0]=='>'):
                tr0 += 1
                if EndRead1 and tr0>EndRead1: break
                if (not(Interleave1) or (tr0-1)%Threads1==(myThread1-1) or L1=='>@') and tr0>StartRead1:
                    if tr1 and tr1%ReportGranularity1==0:
                        LogNote1('ScanRabbit interim report, Substrate file(s):',F1n,F2n*bool(F2n),' Reads:',tr1, ' KeptReads:', nr1)
                    v0=0; a0=0; v3=0; a3=0
                    if accu1:
                        if DataCycle1==0:
                            accu1[1] = ''.join(accu1[1])
                            accu2[1] = ''.join(accu2[1])
                        s1 = accu1[1]
                        s3 = accu2[1]
                        if len(s1)>=BlockLength3:
                            for j1,c1 in enumerate(s1):
                                cd1 = BaseL1[ord(c1)]
                                v0 = ((v0&mask1)<<2)+cd1
                                a0 = (a0>>2)+((3-cd1)<<ksam1)
                                if (j1>=StartCheckingSense1) and PreFilter1(v0) and (j1+CoreEndToBlockEnd3<len(s1)):
                                    s0 = s1[j1-StartCheckingSense1:j1+CoreEndToBlockEnd3+1]
                                    score100 = FullValue1(s0)
                                    if score100>minMatchScore100:
                                        score1 = score100/100
                                        pep1 = translate1(s0)
                                        if enforceCore1 and not(re.match(CoreRE1,pep1[BlockToCoreOffset1:])): continue
                                        SM3 = max((j1+1)%3,j1-StartCheckingSense1-MaxFlankingLength1)
                                        SM4 = min(len(s1),j1+CoreEndToBlockEnd3+1+MaxFlankingLength1)
                                        ups1 = translate1(s1[SM3:j1-StartCheckingSense1])
                                        dwn1 = translate1(s1[j1+CoreEndToBlockEnd3+1:SM4])
                                        pro1 = ups1.lower()+pep1+dwn1.lower()
                                        KeepRead1 += '__Match_Candidate_from_'+F1n+'_sense_@base_'+str(j1-klen1-1)+'_Peptide='+pro1+'_Score='+'%.3f'%score1
                                        if not(pep1 in BlockD1):
                                            BlockD1[pep1] = BlockInstance1(pep1,score1,ups1,dwn1,(i1+1,j1+1))
                                        else:
                                            BlockD1[pep1].add(ups1,dwn1,(i1+1,j1+1))
                                        MaxScore1 = max(score1,MaxScore1)
                                        if MatchByMatchOutput1:
                                            MatchByMatchFile1.write(F1n+'\t'+
                                               (F2n+'\t')*bool(AnyR2File1)+
                                               os.path.dirname(os.path.abspath(F1))+'\t'+
                                               accu1[0].strip().lstrip('>')+'\t'+
                                               'R1'+'\t'+
                                               'plus'+'\t'+
                                                str(j1-klen1-1)+'\t'+
                                                s0+'\t'+
                                                ups1+'\t'+
                                                pep1+'\t'+
                                                dwn1+'\t'+
                                               '%.3f'%score1+
                                                ('\t'+s1+('\t'+s3)*bool(AnyR2File1))*MatchByMatchWholeSequence1+
                                                Delimiter1)
                                if (j1>=StartCheckingAnti1) and PreFilter1(a0) and (j1+BlockStartToCoreStart3<len(s1)):                                
                                    s0 = s1[j1-StartCheckingAnti1:j1+BlockStartToCoreStart3+1]
                                    score100 = FullValueA1(s0)
                                    if score100>minMatchScore100:
                                        score2 = score100/100
                                        pep2 = translate1(vantisense(s0))
                                        if enforceCore1 and not(re.match(CoreRE1,pep2[BlockToCoreOffset1:])): continue
                                        SM3 = min(j1+BlockStartToCoreStart3+1+MaxFlankingLength1,len(s1))
                                        SM4 = max(0,j1-StartCheckingAnti1-MaxFlankingLength1)
                                        ups2 = translate1(vantisense(s1[j1+BlockStartToCoreStart3+1:SM3-(SM3-(j1+1))%3]))
                                        dwn2 = translate1(vantisense(s1[SM4:j1-StartCheckingAnti1]))
                                        pro2 = ups2.lower()+pep2+dwn2.lower()
                                        KeepRead1 += '__Match_Candidate_from_'+F1n+'_antisense_@base_'+str(j1+3)+'_Peptide='+pro2+'_Score='+'%.3f'%score2
                                        if not(pep2 in BlockD1):
                                            BlockD1[pep2] = BlockInstance1(pep2,score2,ups2,dwn2,(i1+1,-j1-1))
                                        else:
                                            BlockD1[pep2].add(ups2,dwn2,(i1+1,-j1-1))
                                        MaxScore1 = max(score2,MaxScore1)
                                        if MatchByMatchOutput1:
                                            MatchByMatchFile1.write(F1n+'\t'+
                                               (F2n+'\t')*bool(AnyR2File1)+
                                               os.path.dirname(os.path.abspath(F1))+'\t'+
                                               accu1[0].strip().lstrip('>')+'\t'+
                                               'R1'+'\t'+
                                               'minus'+'\t'+
                                                str(j1+3)+'\t'+
                                                vantisense(s0)+'\t'+
                                                ups2+'\t'+
                                                pep2+'\t'+
                                                dwn2+'\t'+
                                               '%.3f'%score2+
                                                ('\t'+s1+('\t'+s3)*bool(AnyR2File1))*MatchByMatchWholeSequence1+
                                                Delimiter1)
                        if len(s3)>=BlockLength3:                
                            for j3,c3 in enumerate(s3):
                                cd3 = BaseL1[ord(c3)]
                                v3 = ((v3&mask1)<<2)+cd3
                                a3 = (a3>>2)+((3-cd3)<<ksam1)
                                if (j3>=StartCheckingSense1) and PreFilter1(v3) and (j3+CoreEndToBlockEnd3<len(s3)):
                                    s0 = s3[j3-StartCheckingSense1:j3+CoreEndToBlockEnd3+1]
                                    score100 = FullValue1(s0)
                                    if score100>minMatchScore100:
                                        score3 = score100/100
                                        pep3 = translate1(s3[j3-StartCheckingSense1:j3+CoreEndToBlockEnd3+1])
                                        if enforceCore1 and not(re.match(CoreRE1,pep3[BlockToCoreOffset1:])): continue
                                        SM3 = max((j3+1)%3,j3-StartCheckingSense1-MaxFlankingLength1)
                                        SM4 = min(len(s3),j3+CoreEndToBlockEnd3+1+MaxFlankingLength1)
                                        ups3 = translate1(s3[SM3:j3-StartCheckingSense1])
                                        dwn3 = translate1(s3[j3+CoreEndToBlockEnd3+1:SM4])
                                        pro3 = ups3.lower()+pep3+dwn3.lower()
                                        KeepRead3 += '__Match_Candidate_from_'+F2n+'_sense_@base_'+str(j3-klen1-1)+'_Peptide='+pro3+'_Score='+'%.3f'%score3
                                        if not(pep3 in BlockD1):
                                            BlockD1[pep3] = BlockInstance1(pep3,score3,ups3,dwn3,(-i1-1,j3+1))
                                        else:
                                            BlockD1[pep3].add(ups3,dwn3,(-i1-1,j3+1))
                                        MaxScore1 = max(score3,MaxScore1)
                                        if MatchByMatchOutput1:
                                            MatchByMatchFile1.write(F1n+'\t'+
                                               F2n+'\t'+
                                               os.path.dirname(os.path.abspath(F1))+'\t'+
                                               accu2[0].strip().lstrip('>')+'\t'+
                                               'R2'+'\t'+
                                               'plus'+'\t'+
                                                str(j3-klen1-1)+'\t'+
                                                s0+'\t'+
                                                ups3+'\t'+
                                                pep3+'\t'+
                                                dwn3+'\t'+
                                               '%.3f'%score3+
                                                ('\t'+s1+'\t'+s3)*MatchByMatchWholeSequence1+
                                                Delimiter1)
                                if (j3>=StartCheckingAnti1) and  PreFilter1(a3) and (j3+BlockStartToCoreStart3<len(s3)):
                                    s0 = s3[j3-StartCheckingAnti1:j3+BlockStartToCoreStart3+1]
                                    score100 = FullValueA1(s0)
                                    if score100>minMatchScore100:
                                        score4 = score100/100
                                        pep4 = translate1(vantisense(s0))
                                        if enforceCore1 and not(re.match(CoreRE1,pep4[BlockToCoreOffset1:])): continue
                                        SM3 = min(j3+BlockStartToCoreStart3+1+MaxFlankingLength1,len(s3))
                                        SM4 = max(0,j3-StartCheckingAnti1-MaxFlankingLength1)
                                        ups4 = translate1(vantisense(s3[j3+BlockStartToCoreStart3+1:SM3-(SM3-(j3+1))%3]))
                                        dwn4 = translate1(vantisense(s3[:SM4]))
                                        pro4 = ups4.lower()+pep4+dwn4.lower()
                                        KeepRead3 += '__Match_Candidate_from_'+F2n+'_antisense_@base_'+str(j3+3)+'_Peptide='+pro4+'_Score='+'%.3f'%score4
                                        if not(pep4 in BlockD1):
                                            BlockD1[pep4] = BlockInstance1(pep4,score4,ups4,dwn4,(-i1-1,-j3-1))
                                        else:
                                            BlockD1[pep4].add(ups4,dwn4,(-i1-1,-j3-1))
                                        MaxScore1 = max(score4,MaxScore1)
                                        if MatchByMatchOutput1:
                                            MatchByMatchFile1.write(F1n+'\t'+
                                               F2n+'\t'+
                                               os.path.dirname(os.path.abspath(F1))+'\t'+
                                               accu2[0].strip().lstrip('>')+'\t'+
                                               'R2'+'\t'+
                                               'minus'+'\t'+
                                                str(j3+3)+'\t'+
                                                vantisense(s0)+'\t'+
                                                ups4+'\t'+
                                                pep4+'\t'+
                                                dwn4+'\t'+
                                               '%.3f'%score4+
                                                ('\t'+s1+'\t'+s3)*MatchByMatchWholeSequence1+
                                                Delimiter1)
                        if KeepRead1 or KeepRead3:
                            NR1 += 1; nr1+=1
                            if ProvideReads1:
                                accu1 = [x.strip() for x in accu1]; accu2 = [x.strip() for x in accu2]
                                accu1[0]+=KeepRead1; accu2[0]+=KeepRead3
                                OutFile1.write(Delimiter1.join(accu1)+Delimiter1)
                                if AnyR2File1:
                                    OutFile2.write(Delimiter1.join(accu2)+Delimiter1)
                    L1strip1 = L1.strip()
                    accu1 = [L1strip1]
                    accu2 = [L2.strip()]
                    if L1strip1!='>@':
                        TR1 += 1
                        tr1+=1
                    if DataCycle1==0:
                        accu1.append([])
                        accu2.append([])
                    KeepRead1 = '' ## default is to toss
                    KeepRead3 = '' ## default is to toss
            elif (not(Interleave1) or (tr0-1)%Threads1==(myThread1-1)) and tr0>StartRead1:
                if DataCycle1==0:
                    accu1[1].append(L1.strip())
                    accu2[1].append(L2.strip())
                else:
                    accu1.append(L1.strip())
                    accu2.append(L2.strip())
            if EndChar1 and NowChar1>EndChar1:
                break
            NowChar1 += len(L1)

        F1o.close()
        if F2n:
            F2o.close()
        MaxScoreChaser1 = bool(MaxScore1)*('. MaxScore=%.3f'%MaxScore1)
        LogNote1('ScanRabbit File Finished', F1n,F2n*bool(F2n),'. Caught',nr1,'match of',tr1,'candidates',MaxScoreChaser1)
        if SummaryOutput1:
            for b1 in BlockD1:
                BlockD1[b1].countme()
            BList1 = sorted(list(BlockD1.keys()), key=lambda x: (BlockD1[x].score,BlockD1[x].count), reverse=True)
            if MatchSpeciesToReport1==0:
                DetailedList1 = [BlockD1[b1].report() for b1 in BList1]
                SpecList1 = [(BlockD1[b1].score,BlockD1[b1].count) for x in BList1]
            else:
                BList2 = sorted(list(BlockD1.keys()), key=lambda x: (BlockD1[x].count,BlockD1[x].score), reverse=True)
                BSet1 = set(BList1[:MatchSpeciesToReport1])
                BSet2 = set(BList2[:MatchSpeciesToReport1])
                BSetAll1 = sorted(list(BSet1.union(BSet2)), key=lambda x: (BlockD1[x].score,BlockD1[x].count), reverse=True)
                for b1 in BSetAll1:
                    BlockD1[b1].assemble()
                DetailedList1 = [BlockD1[b1].report() for b1 in BSetAll1]
                SpecList1 = [(BlockD1[b1].score,BlockD1[b1].count) for b1 in BSetAll1]
            if len(accu1)>1:
                LastL1len = len(''.join(accu1[1]))
            else:
                LastL1len = 0
            if len(accu2)>1:
                LastL2len = len(''.join(accu2[1]))
            else:
                LastL2len = 0
            if AbbreviateSummary1:
                SummaryFile1.write(F1n+'\t'+
                                   (F2n+'\t')*bool(AnyR2File1)+
                                   os.path.dirname(os.path.abspath(F1))+'\t'+
                                   str(nr1)+'\t'+
                                   '%.3f'%MaxScore1+'\t'+
                                   str(tr1)+'\t'+
                                   str(LastL1len)+'\t'+
                                   (str(LastL2len)+'\t')*bool(AnyR2File1)+
                                   ' ; '.join(DetailedList1)+Delimiter1)
            else:
                SummaryFile1.write(F1n+'\t'+
                                   (F2n+'\t')*bool(AnyR2File1)+
                                   os.path.dirname(os.path.abspath(F1))+'\t'+
                                   str(nr1)+'\t'+
                                   '%.3f'%MaxScore1+'\t'+
                                   str(tr1)+'\t'+
                                   str(LastL1len)+'\t'+
                                   (str(LastL2len)+'\t')*bool(AnyR2File1)+Delimiter1)
                for Spec1,Detail1 in zip(SpecList1,DetailedList1):
                    SummaryFile1.write(F1n+'\t'+
                                   (F2n+'\t')*bool(AnyR2File1)+
                                   os.path.dirname(os.path.abspath(F1))+'\t'+
                                   str(Spec1[1])+'\t'+
                                   '%.3f'%Spec1[0]+'\t'+
                                   str(tr1)+'\t'+
                                   str(LastL1len)+'\t'+
                                   (str(LastL2len)+'\t')*bool(AnyR2File1)+
                                   Detail1+Delimiter1)
        if SummaryOutput1:
            SummaryFile1.flush()
    except:
        LogNote1('!!!!!! Read failure or disk error with files',F1n,F2n,'will continue with remaining files')

LogNote1('ScanRabbit run finished. Caught',NR1,'match of',TR1,'total candidates')
if ProvideReads1:
    OutFile1.close()
    if AnyR2File1:
        OutFile2.close()
        LogNote1('Reads written to:',OutFileName1,',',OutFileName2)
    else:
        LogNote1('Reads written to:',OutFileName1)
if SummaryOutput1:
    SummaryFile1.close()
    LogNote1('FileByFile summary written to:',SummaryFileName1)
if MatchByMatchOutput1:
    MatchByMatchFile1.close()
    LogNote1('MatchByMatch summary written to:',MatchByMatchFileName1)
try:
    pypyreminder1('after')
except:
    pass
