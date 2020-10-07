library(hydra)
library(ape)
library(MASS)

Ntrees <- 100
Nleaves <- 100 + 1
Ms <- c(5,10,20,30)

Poincare_dist <- function(u,v){
	if(sum(u**2)>=1 || sum(v**2)>=1){
		return(Inf)
	}else{
		return(acosh(1+2*(sum((u-v)**2))/(1-sum(u**2))/(1-sum(v**2))))
	}
}

objfn_H1 <- function(tmpz){
	ret <- 0
	for(i in 1:Nleaves){
		a <- Poincare_dist(Z1[i,], tmpz)
		ret <- ret + (a - Dall[i,Nleaves+n])**2
	}
	return(ret)
}

objfn_H2 <- function(tmpz){
	ret <- 0
	for(i in 1:Nleaves){
		a <- Poincare_dist(Z2[i,], tmpz)
		ret <- ret + (log(cosh(a)) - Dall[i,Nleaves+n])**2
	}
	return(ret)
}

objfn_MDS <- function(tmpz){
	ret <- 0
	for(i in 1:Nleaves){
		a <- sqrt(sum((X1.mds$points[i,]-tmpz)**2))
		ret <- ret + (a - Dall[i,Nleaves+n])**2
	}
	return(ret)
}

for(nt in 1:Ntrees){
	tree <- read.tree(file=paste("data/tree_",nt,".txt",sep=""))
	tree$mrca <- mrca(tree)
	Dall <- dist.nodes(tree)
	Dleaves <- cophenetic(tree)

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

		Z1inter <- matrix(0,Nall-Nleaves,M)
		for(n in 1:(Nall-Nleaves)){
			tmp <- optim(rep(0,M),objfn_H1,gr=NULL,method="BFGS")
			Z1inter[n,] <- tmp$par
		}

		Z1inter_dist <- matrix(0, Nall, Nall)
		for(i in 1:Nall){
			for(j in 1:Nall){
				if(i > Nleaves){
					tmpz1 <- Z1inter[i-Nleaves,]
				}else{
					tmpz1 <- Z1[i,]
				}

				if(j > Nleaves){
					tmpz2 <- Z1inter[j-Nleaves,]
				}else{
					tmpz2 <- Z1[j,]
				}

				Z1inter_dist[i,j] <- Poincare_dist(tmpz1, tmpz2)
			}
		}

		#our hyperbolic embeddings
		X2.hydra <- hydraPlus(X2, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z2 <- X2.hydra$r * X2.hydra$directional

		Z2inter <- matrix(0,Nall-Nleaves,M)
		for(n in 1:(Nall-Nleaves)){
			tmp <- optim(rep(0,M),objfn_H2,gr=NULL,method="BFGS")
			Z2inter[n,] <- tmp$par
		}

		Z2inter_dist <- matrix(0, Nall, Nall)
		for(i in 1:Nall){
			for(j in 1:Nall){
				if(i > Nleaves){
					tmpz1 <- Z2inter[i-Nleaves,]
				}else{
					tmpz1 <- Z2[i,]
				}

				if(j > Nleaves){
					tmpz2 <- Z2inter[j-Nleaves,]
				}else{
					tmpz2 <- Z2[j,]
				}

				Z2inter_dist[i,j] <- Poincare_dist(tmpz1, tmpz2)
			}
		}

		#Euclidean embeddings
		X1.mds <- sammon(X1, k=M)

		Zmdsinter <- matrix(0,Nall-Nleaves,M)
		for(n in 1:(Nall-Nleaves)){
			tmp <- optim(rep(0,M),objfn_MDS,gr=NULL,method="BFGS")
			Zmdsinter[n,] <- tmp$par
		}

		Zmdsinter_dist <- matrix(0, Nall, Nall)
		for(i in 1:Nall){
			for(j in 1:Nall){
				if(i > Nleaves){
					tmpz1 <- Zmdsinter[i-Nleaves,]
				}else{
					tmpz1 <- X1.mds$points[i,]
				}

				if(j > Nleaves){
					tmpz2 <- Zmdsinter[j-Nleaves,]
				}else{
					tmpz2 <- X1.mds$points[j,]
				}

				Zmdsinter_dist[i,j] <- sqrt(sum((tmpz1-tmpz2)**2))
			}
		}

		#MSE
		MSE_H1[which(Ms==M)] <- sum((c(Dall[upper.tri(Dall)])-c(Z1inter_dist[upper.tri(Dall)]))**2) / choose(Nall,2)
		MSE_H2[which(Ms==M)] <- sum((c(Dall[upper.tri(Dall)])-c(log(cosh(Z2inter_dist[upper.tri(Dall)]))))**2) / choose(Nall,2)
		MSE_MDS[which(Ms==M)] <- sum((c(Dall[upper.tri(Dall)])-c(Zmdsinter_dist[upper.tri(Dall)]))**2) / choose(Nall,2)
	}

	output <- rbind(MSE_MDS, MSE_H1, MSE_H2)
	rownames(output) <- c("MDS","H1","H2")
	colnames(output) <- Ms
	write.csv(output, paste("out/MSE_Dint_",nt,".csv",sep=""))
}


#merge result
output <- matrix(0, 3, length(Ms))
rownames(output) <- c("MDS","H1","H2")
colnames(output) <- Ms
for(i in 1:Ntrees){
	tmp <- read.csv(paste("out/MSE_Dint_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/MSE_Dint_all.csv")
