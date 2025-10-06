import pysam
from argparse import ArgumentParser
import sys
import subprocess

parser = ArgumentParser()
#~ parser.add_argument('--indexes', dest='inFile', required=True,
                    #~ help='Path to indexes file, one per line.')
parser.add_argument('--in', dest='inPrefix', required=True)
o=parser.parse_args()

#~ outFile = open("summary.csv",'wba')

#~ indexes = []
#~ indexFile = open(o.inFile, 'rb')
#~ for line in indexFile:
    #~ indexes.append(line.strip())
#~ indexFile.close()
#~ outFile.write("Index,Raw reads,SSCS reads,Mapped SSCS,DCS reads,Mapped DCS,Peak Family Size,Max Family Size,SSCS On Target,DCS On Target,DCS Mean Depth,DCS Max Depth,DCS Uncovered Target,Nucleotides Sequenced,Mutations,,A>T,A>C,A>G,T>A,T>C,T>G,C>A,C>T,C>G,G>A,G>T,G>C\n")

index=o.inPrefix
#~ sys.stderr.write("Index %s" % index)
sys.stderr.write("Reading config")
# Get the run ID from the config file
#~ runID=""
c=0
C=0.1
d=100
#~ configFile = open("%s_config.sh" % (index), 'rb')
#~ for line in configFile:
	#~ if "RUN_ID=" in line:
		#~ runID = line.strip().split('=')[1].strip('"')
	#~ elif "minClonal=" in line:
		#~ c=line.strip().split('=')[1].split()[0]
	#~ elif "maxClonal=" in line:
		#~ C=line.strip().split('=')[1].split()[0]
	#~ elif "minDepth=" in line:
		#~ d=line.strip().split('=')[1].split()[0]
#~ configFile.close()

sys.stderr.write("Getting read counts")
# get read counts
#rawReads = pysam.flagstat("%s.sort.bam" % (index)).split('\n')[0].split()[0]
#sscsFlagstat=pysam.flagstat("%s_mem.sscs.sort.bam" % (index)).split('\n')
#sscsReads=sscsFlagstat[0].split()[0]
#mappedSscs=sscsFlagstat[4].split()[0]
rawReads=0
sscsReads=0
mappedSscs=0

dcsFlagstat=pysam.flagstat("%s.sort.bam" % (index)).split('\n')
dcsReads=dcsFlagstat[0].split()[0]
mappedDcs=dcsFlagstat[4].split()[0]

#~ sys.stderr.write("Processing Tagstats")
#~ # get tagstats numbers
#~ tagstatsFile = open("%s.tagstats.txt" % (index), 'rb')
lastProportion=1
peakProportion = 0
peakSize = 1
maxSize=0
#~ for line in tagstatsFile:
	#~ if float(line.split()[2]) <= lastProportion:
		#~ lastProportion = float(line.split()[2])
	#~ elif float(line.split()[2]) >= peakProportion:
		#~ lastProportion = 0
		#~ peakSize = line.split()[0]
		#~ peakProportion = float(line.split()[2])
	#~ maxSize = line.split()[0]
#~ tagstatsFile.close()

# get hs_metrics numbers
#~ hsSscsFile = open("%s.sscs.filt.no_overlap.hs_metrics.txt" % (index), 'rb')
#~ lineCtr=0
#~ line = hsSscsFile.readline()
#~ while lineCtr < 7:
	#~ line = hsSscsFile.readline()
	#~ lineCtr += 1
#~ sscsOnTarget=line.split()[18]
#~ hsSscsFile.close()
#~ sscsOnTarget = 0
#~ sys.stderr.write("Processing hs_metrics")
#hsDcsFile = open("%s.hs_metrics.txt" % (index), 'rb')
hsDcsFile = open("%s.hs_metrics.txt" % (index), 'r')
lineCtr=0
line = hsDcsFile.readline()
while lineCtr < 7:
	line = hsDcsFile.readline()
	lineCtr += 1
rawReads=line.split('\t')[25]
#dcsOnTarget=line.split()[18]
dcsOnTarget=line.split('\t')[3]
#dcsMeanDepth=line.split()[22]
dcsMeanDepth=line.split('\t')[33]
#dcsMaxDepth=line.split()[24]
dcsMaxDepth=line.split('\t')[35]
#dcsUncovered=line.split()[28]
dcsUncovered=line.split('\t')[43]
hsDcsFile.close()

sys.stderr.write("Processing countmuts")
# get countmuts data
#cmFile = open("%s.region.c%s-%s.d%s.countmuts.txt" % (index,c,C,d), 'rb')
cmFile = open("%s.region.c%s-%s.d%s.countmuts.txt" % (index,c,C,d), 'r')

AtoT=""
AtoC=""
AtoG=""
TtoA=""
TtoC=""
TtoG=""
CtoA=""
CtoT=""
CtoG=""
GtoA=""
GtoT=""
GtoC=""
totalNt=""
totalMuts=""
ins=""
dels=""

for line in cmFile:
	if "A to T" in line:
		AtoT=line.split()[3]
	elif "A to C" in line:
		AtoC=line.split()[3]
	elif "A to G" in line:
		AtoG=line.split()[3]
	elif "T to A" in line:
		TtoA=line.split()[3]
	elif "T to C" in line:
		TtoC=line.split()[3]
	elif "T to G" in line:
		TtoG=line.split()[3]
	elif "C to A" in line:
		CtoA=line.split()[3]
	elif "C to T" in line:
		CtoT=line.split()[3]
	elif "C to G" in line:
		CtoG=line.split()[3]
	elif "G to A" in line:
		GtoA=line.split()[3]
	elif "G to T" in line:
		GtoT=line.split()[3]
	elif "G to C" in line:
		GtoC=line.split()[3]
#	elif "Total nucleotides sequenced:" in line:
#		totalNt = line.strip().split()[3]
	elif "Total point mutations:" in line:
		totalMuts = line.strip().split()[3]
	elif "Total insertion events:" in line:
		ins=line.strip().split()[3]
	elif "Total deletion events:" in line:
		dels=line.strip().split()[3]

pileup_count = subprocess.run(['wc', '-l', f"{index}.region.pileup"], capture_output=True, text=True)
totalNt = int(pileup_count.stdout.strip().split()[0])

cmFile.close()
sys.stdout.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (index,rawReads,dcsOnTarget,dcsMeanDepth,dcsMaxDepth,dcsUncovered,totalNt,totalMuts,AtoT,AtoC,AtoG,TtoA,TtoC,TtoG,CtoA,CtoT,CtoG,GtoA,GtoT,GtoC,ins,dels))
            
#~ outFile.close()
        
