import sys
import os
import pysam
import random
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--input', dest='in_bam', required=True,
                    help='Path to input bam file.')
parser.add_argument('--prefix', dest='out', required=True, 
                    help='Prefix for outputs')
parser.add_argument('--thresholds', dest='thresh', required=True)
o=parser.parse_args()

in_bam_file = pysam.AlignmentFile(o.in_bam, "rb", check_sq=False)

preThresh = sorted([float(x) for x in o.thresh.split(',')])
thresholds = [0]
for x in preThresh:
    thresholds.append(thresholds[-1] + x)
thresholds = sorted(thresholds, reverse=True)
sys.stderr.write("Thresholds: ")
sys.stderr.write(" ".join([str(x) for x in thresholds]))
sys.stderr.write("\n")
out_bam = [pysam.AlignmentFile("%s_s%s.bam" % (o.out, x), 'wb', template=in_bam_file) for x in range(len(thresholds))]


processed = 0
written = [0 for x in thresholds]


for line in in_bam_file.fetch(until_eof=True):
    randTest = random.random()
    for i in range(len(thresholds)):
        if randTest >= thresholds[i]:
            out_bam[i].write(line)
            written[i] += 1
            break
    
    processed += 1
    if processed%100000 == 0:
        sys.stderr.write("%s/%s written\n" % (str(written), processed))
in_bam_file.close()
for i in range(len(thresholds)):
    out_bam[i].close()

