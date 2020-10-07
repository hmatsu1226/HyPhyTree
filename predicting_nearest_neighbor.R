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

	#our hyperbolic embeddings
	output <- c()
	for(M in Ms){
		X2.hydra <- hydraPlus(X2, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z2 <- X2.hydra$r * X2.hydra$directional
		X2.hydra$dist <- hydra:::get.distance(X2.hydra$r, X2.hydra$directional, X2.hydra$curvature)

		Z2inter <- matrix(0,length(sub_inter),M)
		for(i in 1:length(sub_inter)){
			n <- sub_inter[i]
			tmp <- optim(rep(0,M),objfn_H2,gr=NULL,method="BFGS")
			Z2inter[i,] <- tmp$par
		}

		n1s <- setdiff(1:Nleaves, sub_leaves)
		res_rank <- rep(1, length(n1s))
		for(n1 in n1s){
			tmp <- order(Dall[n1,c(sub_leaves, sub_inter)])[1]
			nearest_point <- c(sub_leaves, sub_inter)[tmp]
			nearest_dist <- Inf
			if(nearest_point <= Nleaves){
				nearest_dist <- X2.hydra$dist[n1,nearest_point]
			}else{
				nearest_dist <- Poincare_dist(Z2[n1,], Z2inter[which(sub_inter==nearest_point),])
			}

			for(i in 1:length(sub_leaves)){
				if(sub_leaves[i] == nearest_point){
					next
				}
				if(X2.hydra$dist[n1,sub_leaves[i]] < nearest_dist){
					res_rank[which(n1s==n1)] <- res_rank[which(n1s==n1)] + 1
				}
			}

			for(i in 1:length(sub_inter)){
				if(sub_inter[i] == nearest_point){
					next
				}
				if(Poincare_dist(Z2[n1,], Z2inter[i,]) < nearest_dist){
					res_rank[which(n1s==n1)] <- res_rank[which(n1s==n1)] + 1
				}
			}
		}
	
		output <- rbind(output, c(M, mean(res_rank)))
		print(paste(M, sep="  "))
		print(mean(res_rank))
	}

	write.csv(output, paste("out/Result_NN_H2_",nt,".csv",sep=""))


	#general hyperbolic embeddings
	output <- c()
	for(M in Ms){
		X1.hydra <- hydraPlus(X1, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z1 <- X1.hydra$r * X1.hydra$directional
		X1.hydra$dist <- hydra:::get.distance(X1.hydra$r, X1.hydra$directional, X1.hydra$curvature)

		Z1inter <- matrix(0,length(sub_inter),M)
		for(i in 1:length(sub_inter)){
			n <- sub_inter[i]
			tmp <- optim(rep(0,M),objfn_H1,gr=NULL,method="BFGS")
			Z1inter[i,] <- tmp$par
		}

		n1s <- setdiff(1:Nleaves, sub_leaves)
		res_rank <- rep(1, length(n1s))
		for(n1 in n1s){
			tmp <- order(Dall[n1,c(sub_leaves, sub_inter)])[1]
			nearest_point <- c(sub_leaves, sub_inter)[tmp]
			nearest_dist <- Inf
			if(nearest_point <= Nleaves){
				nearest_dist <- X1.hydra$dist[n1,nearest_point]
			}else{
				nearest_dist <- Poincare_dist(Z1[n1,], Z1inter[which(sub_inter==nearest_point),])
			}

			for(i in 1:length(sub_leaves)){
				if(sub_leaves[i] == nearest_point){
					next
				}
				if(X1.hydra$dist[n1,sub_leaves[i]] < nearest_dist){
					res_rank[which(n1s==n1)] <- res_rank[which(n1s==n1)] + 1
				}
			}

			for(i in 1:length(sub_inter)){
				if(sub_inter[i] == nearest_point){
					next
				}
				if(Poincare_dist(Z1[n1,], Z1inter[i,]) < nearest_dist){
					res_rank[which(n1s==n1)] <- res_rank[which(n1s==n1)] + 1
				}
			}
		}
	
		output <- rbind(output, c(M, mean(res_rank)))
		print(paste(M, sep="  "))
		print(mean(res_rank))
	}

	write.csv(output, paste("out/Result_NN_H1_",nt,".csv",sep=""))


	#Euclidean embeddings
	output <- c()
	for(M in Ms){
		X1.mds <- sammon(X1, k=M)
		Zmds <- X1.mds$points
		X1.mds$dist <- matrix(0, Nleaves, Nleaves)
		for(i in 1:Nleaves){
			for(j in 1:Nleaves){
				X1.mds$dist[i,j] <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[j,])**2))
			}
		}

		Zmdsinter <- matrix(0,length(sub_inter),M)
		for(i in 1:length(sub_inter)){
			n <- sub_inter[i]
			tmp <- optim(rep(0,M),objfn_MDS,gr=NULL,method="BFGS")
			Zmdsinter[i,] <- tmp$par
		}

		n1s <- setdiff(1:Nleaves, sub_leaves)
		res_rank <- rep(1, length(n1s))
		for(n1 in n1s){
			tmp <- order(Dall[n1,c(sub_leaves, sub_inter)])[1]
			nearest_point <- c(sub_leaves, sub_inter)[tmp]
			nearest_dist <- Inf
			if(nearest_point <= Nleaves){
				nearest_dist <- X1.mds$dist[n1,nearest_point]
			}else{
				nearest_dist <- sqrt(sum((Zmds[n1,]-Zmdsinter[which(sub_inter==nearest_point),])**2))
			}

			for(i in 1:length(sub_leaves)){
				if(sub_leaves[i] == nearest_point){
					next
				}
				if(X1.mds$dist[n1,sub_leaves[i]] < nearest_dist){
					res_rank[which(n1s==n1)] <- res_rank[which(n1s==n1)] + 1
				}
			}

			for(i in 1:length(sub_inter)){
				if(sub_inter[i] == nearest_point){
					next
				}
				if(sqrt(sum((Zmds[n1,]-Zmdsinter[i,])**2)) < nearest_dist){
					res_rank[which(n1s==n1)] <- res_rank[which(n1s==n1)] + 1
				}
			}			
		}
	
		output <- rbind(output, c(M, mean(res_rank)))
		print(paste(M, sep="  "))
		print(mean(res_rank))
	}

	write.csv(output, paste("out/Result_NN_MDS_",nt,".csv",sep=""))
}



#merge results
output <- read.csv("out/Result_NN_MDS_1.csv", row.names=1)
for(i in 2:Ntrees){
	tmp <- read.csv(paste("out/Result_NN_MDS_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/Result_NN_MDS_all.csv")

output <- read.csv("out/Result_NN_H1_1.csv", row.names=1)
for(i in 2:Ntrees){
	tmp <- read.csv(paste("out/Result_NN_H1_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/Result_NN_H1_all.csv")

output <- read.csv("out/Result_NN_H2_1.csv", row.names=1)
for(i in 2:Ntrees){
	tmp <- read.csv(paste("out/Result_NN_H2_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/Result_NN_H2_all.csv")

