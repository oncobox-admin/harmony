library(preprocessCore)
library(HARMONY)

prefix = "/home/R/harmony/"

IFN = paste(prefix, "P0.txt", sep = "")
P0 = read.table(IFN, header = TRUE)

IFN = paste(prefix, "Input.txt", sep = "")
Input = read.table(IFN, header = TRUE)

N = ncol(Input)
CN0 = colnames (Input)

CN = CN0[-1]

SYMBOL = as.vector(Input[,1])

set.seed(1234)

for ( i in 1:(N-1)) {
  print (i)
   
  VECTOR = Input[,(i+1)]

  sample = cbind(SYMBOL,VECTOR)

  P = merge(P0,sample,by = "SYMBOL")
 
  RN = as.vector(P[,1])
  
  P = P[,-1]
  
  FN = paste(prefix, "add.txt", sep = "")
  write.table (P,FN)
  P = read.table(FN)
  
  P = log(P)

  P = as.matrix(P)

  P = normalize.quantiles(P) 
  
  write.table(P,FN,row.names = RN)
  P = read.table(FN)

  P = as.data.frame(P)

  result = harmony.afx.static(P, gene.cluster = "skmeans", assay.cluster = "hclust", iterations=3)

  Q = ncol(result)

  harmonized = result[,Q]

  genes <- intersect(rownames(XPN.AFX.RDA), RN)

  if (i == 1) output <- harmonized
  if (i > 1)  output <- cbind(output,harmonized)
 
}

output = exp(output)
SYMBOL1 = rownames(output)
output = cbind(SYMBOL1,output)
CN = c("SYMBOL", CN)
OFN = paste(prefix, "Output.txt",  sep = "")  
write.table(output, OFN, col.names = CN, row.names = FALSE, sep = "\t")
