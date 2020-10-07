library(hydra)
library(ape)
library(MASS)

Ntrees <- 100
Nleaves <- 100 + 1
Ms <- c(5,10,20,30)

#outgroup node idx
o <- Nleaves

Poincare_dist <- function(u,v){
	if(sum(u**2)>=1 || sum(v**2)>=1){
		return(Inf)
	}else{
		return(acosh(1+2*(sum((u-v)**2))/(1-sum(u**2))/(1-sum(v**2))))
	}
}

objfn_H1 <- function(tmpz){
	ret <- 0
	for(i in 1:length(sub_leaves)){
		a <- Poincare_dist(Z1[sub_leaves[i],], tmpz)
		ret <- ret + (a - Dall[sub_leaves[i],n])**2
	}
	return(ret)
}

objfn_H2 <- function(tmpz){
	ret <- 0
	for(i in 1:length(sub_leaves)){
		a <- Poincare_dist(Z2[sub_leaves[i],], tmpz)
		ret <- ret + (log(cosh(a)) - Dall[sub_leaves[i],n])**2
	}
	return(ret)
}

objfn_MDS <- function(tmpz){
	ret <- 0
	for(i in 1:length(sub_leaves)){
		a <- sqrt(sum((Zmds[sub_leaves[i],]-tmpz)**2))
		ret <- ret + (a - Dall[sub_leaves[i],n])**2
	}
	return(ret)
}

for(nt in 1:Ntrees){
	tree <- read.tree(file=paste("data/tree_",nt,".txt",sep=""))
	tree$mrca <- mrca(tree)
	Dall <- dist.nodes(tree)
	Dleaves <- cophenetic(tree)

	#partial tree
	sub_leaves <- c(o,sample(1:(Nleaves-1))[1:20])
	sub_inter <- c()
	for(i in 1:(length(sub_leaves)-1)){
		for(j in (i+1):length(sub_leaves)){
			sub_inter <- c(sub_inter, tree$mrca[sub_leaves[i],sub_leaves[j]])
		}
	}
	sub_inter <- unique(sub_inter)

	Nall <- dim(Dall)[1]
	X1 <- Dleaves
	X2 <- acosh(exp(Dleaves))

	MSE_H1 <- rep(0, length(Ms))
	MSE_H2 <- rep(0, length(Ms))
	MSE_MDS <- rep(0, length(Ms))

	for(M in Ms){
		#general hyperbolic embeddings
		X1.hydra <- hydraPlus(X1, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z1 <- X1.hydra$r * X1.hydra$directional

		Z1inter <- matrix(0,length(sub_inter),M)
		for(i in 1:length(sub_inter)){
			n <- sub_inter[i]
			tmp <- optim(rep(0,M),objfn_H1,gr=NULL,method="BFGS")
			Z1inter[i,] <- tmp$par
		}

		#our hyperbolic embeddings
		X2.hydra <- hydraPlus(X2, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z2 <- X2.hydra$r * X2.hydra$directional

		Z2inter <- matrix(0,length(sub_inter),M)
		for(i in 1:length(sub_inter)){
			n <- sub_inter[i]
			tmp <- optim(rep(0,M),objfn_H2,gr=NULL,method="BFGS")
			Z2inter[i,] <- tmp$par
		}

		#Euclidean embeddings
		X1.mds <- sammon(X1, k=M)
		Zmds <- X1.mds$points

		Zmdsinter <- matrix(0,length(sub_inter),M)
		for(i in 1:length(sub_inter)){
			n <- sub_inter[i]
			tmp <- optim(rep(0,M),objfn_MDS,gr=NULL,method="BFGS")
			Zmdsinter[i,] <- tmp$par
		}

		#MSE
		mse_H1 <- 0
		mse_H2 <- 0
		mse_MDS <- 0
		for(i in 1:Nleaves){
			if(length(which(sub_leaves == i)) == 1){
				next
			}
			for(j in 1:length(sub_inter)){
				mse_H1 <- mse_H1 + (Dall[i,sub_inter[j]] - Poincare_dist(Z1[i,], Z1inter[j,]))**2
				mse_H2 <- mse_H2 + (Dall[i,sub_inter[j]] - log(cosh(Poincare_dist(Z2[i,], Z2inter[j,]))))**2
				mse_MDS <- mse_MDS + (Dall[i,sub_inter[j]] - sqrt(sum((Zmds[i,]-Zmdsinter[j,])**2)))**2
			}
		}

		MSE_H1[which(Ms==M)] <- mse_H1 / ((Nleaves-length(sub_leaves))*length(sub_inter))
		MSE_H2[which(Ms==M)] <- mse_H2 / ((Nleaves-length(sub_leaves))*length(sub_inter))
		MSE_MDS[which(Ms==M)] <- mse_MDS / ((Nleaves-length(sub_leaves))*length(sub_inter))
	}

	output <- rbind(MSE_MDS, MSE_H1, MSE_H2)
	rownames(output) <- c("MDS","H1","H2")
	colnames(output) <- Ms
	write.csv(output, paste("out/MSE_Dpart_",nt,".csv",sep=""))
}


#merge result
output <- matrix(0, 3, length(Ms))
rownames(output) <- c("MDS","H1","H2")
colnames(output) <- Ms
for(i in 1:Ntrees){
	tmp <- read.csv(paste("out/MSE_Dpart_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/MSE_Dpart_all.csv")
