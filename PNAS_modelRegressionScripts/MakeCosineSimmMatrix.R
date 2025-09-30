cos.sim <- function(A,B) 
{
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}   
myMCS3 = read_csv("/Users/loeblabm11/Desktop/CRC_mcs/mcs3.csv")
samples = myMCS3$Sample
outMatrix = data.frame(matrix(ncol = length(samples), nrow = length(samples)), row.names = samples)
colnames(outMatrix) <- samples
outMatrix
for (samp1 in samples) {
  for (samp2 in samples) {
    outMatrix[samp1, samp2] = cos.sim(myMCS3[myMCS3$Sample == samp1, 2:97], myMCS3[myMCS3$Sample == samp2, 2:97])
  }
}
outMatrix

samples = row.names(scores)
outMatrix = data.frame(matrix(ncol = length(samples), nrow = length(samples)), row.names = samples)
colnames(outMatrix) <- samples
outMatrix
for (samp1 in samples) {
  for (samp2 in samples) {
    outMatrix[samp1, samp2] = cos.sim(scores[samp1,], scores[samp2, ])
  }
}
outMatrix
