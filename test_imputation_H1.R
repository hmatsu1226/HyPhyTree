library(ape)
library(RSpectra)

args <- commandArgs(trailingOnly = T)
L <- as.numeric(args[1])

tree <- read.tree(file="data_imputation/tree_1.txt")

N <- 100 + 1 #Leaves num
r <- 101 #outgroup idx

myfn_H1 <- function(para){
  for(i in 1:length(na_pairs)){
    D2[na_pairs[[i]][1], na_pairs[[i]][2]] <- para[i]
    D2[na_pairs[[i]][2], na_pairs[[i]][1]] <- para[i]
  }
  
  H <- -cosh(D2)

  return(-eigs_sym(H, 2, opts = list(retvec = FALSE), which="SA")$values[1])
}

#
idx1 <- read.table(paste("data_imputation/idx1_L",L,".txt",sep=""))[,1]
idx2 <- read.table(paste("data_imputation/idx2_L",L,".txt",sep=""))[,1]

D <- matrix(NA, N, N)
D[idx1, idx1] <- cophenetic(tree)[idx1, idx1]
D[idx2, idx2] <- cophenetic(tree)[idx2, idx2]

na_pairs <- list()
for(i in 1:(N-1)){
  for(j in (i+1):N){
    if(is.na(D[i,j])){
      na_pairs <- c(na_pairs, list(c(i,j)))
    }
  }
}

D2 <- D
#initialization
tmp <- c()
for(i in 1:length(na_pairs)){
  D2[na_pairs[[i]][1], na_pairs[[i]][2]] <- max(abs(D[na_pairs[[i]][1],]-D[na_pairs[[i]][2],]), na.rm=TRUE)
  D2[na_pairs[[i]][2], na_pairs[[i]][1]] <- D2[na_pairs[[i]][1], na_pairs[[i]][2]]
  tmp <- c(tmp, D2[na_pairs[[i]][1], na_pairs[[i]][2]])
}

optpara <- optim(tmp, myfn_H1, method="L-BFGS-B", lower=0, upper=2*max(D, na.rm = TRUE), control=list(trace=1, maxit=2000))
print(optpara$value)
for(i in 1:length(na_pairs)){
  D2[na_pairs[[i]][1], na_pairs[[i]][2]] <- optpara$par[i]
  D2[na_pairs[[i]][2], na_pairs[[i]][1]] <- optpara$par[i]
}

write.table(D2, paste("out_imputation/optD_L",L,"_H1.txt",sep=""), sep="\t", col.names=F, row.names=F)
