
set -e
set -o pipefail
set -u
set -x

#~ cd 3062583A1_mergeALL_20180129
#~ for x in 1 2 3 4; do
#~ python ../subsampleBam3.py \
#~ --input 3062583A1_mergeALL_20180129.sort.bam \
#~ --prefix 3062583A1_mergeALL_20180129_sub1_${x}_ \
#~ --thresholds "0.096992174,0.215118154,0.256288767"

#~ bash ../processSubsampled.sh 3062583A1_mergeALL_20180129_sub1_${x}_
#~ done
#~ cd ../3070030A1_mergeALL_20180129
#~ for x in 1 2 3 4; do
#~ python ../subsampleBam3.py \
#~ --input 3070030A1_mergeALL_20180129.sort.bam \
#~ --prefix 3070030A1_mergeALL_20180129_sub1_${x}_ \
#~ --thresholds "0.056515262,0.195148173,0.27091011"

#~ bash ../processSubsampled.sh 3070030A1_mergeALL_20180129_sub1_${x}_
#~ done
cd 3070052A1_All
for x in 1 2 3 4; do

python ../subsampleBam3.py \
--input 3070052A1_All.sort.bam \
--prefix 3070052A1_All_sub1_${x}_ \
--thresholds "0.074912688,0.110356413,0.253707493"

bash ../processSubsampled.sh 3070052A1_All_sub1_${x}_
done
cd ../3071007A1_All
for x in 1 2 3 4; do
python ../subsampleBam3.py \
--input 3071007A1_All.sort.bam \
--prefix 3071007A1_All_sub1_${x}_ \
--thresholds "0.078444361,0.140899578,0.319195217"

bash ../processSubsampled.sh 3071007A1_All_sub1_${x}_
done
cd ../
