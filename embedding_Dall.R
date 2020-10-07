library(hydra)
library(ape)
library(MASS)

Ntrees <- 100
Nleaves <- 100 + 1
Ms <- c(5,10,20,30)

for(nt in 1:Ntrees){
	tree <- read.tree(file=paste("data/tree_",nt,".txt",sep=""))
	tree$mrca <- mrca(tree)
	D <- dist.nodes(tree)

	Nall <- dim(D)[1]
	X1 <- D
	X2 <- acosh(exp(D))

	MSE_H1 <- rep(0, length(Ms))
	MSE_H2 <- rep(0, length(Ms))
	MSE_MDS <- rep(0, length(Ms))

	for(M in Ms){
		#general hyperbolic embeddings
		X1.hydra <- hydraPlus(X1, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z1 <- X1.hydra$r * X1.hydra$directional
		X1.hydra$dist <- hydra:::get.distance(X1.hydra$r, X1.hydra$directional, X1.hydra$curvature)

		#our hyperbolic embeddings
		X2.hydra <- hydraPlus(X2, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z2 <- X2.hydra$r * X2.hydra$directional
		X2.hydra$dist <- hydra:::get.distance(X2.hydra$r, X2.hydra$directional, X2.hydra$curvature)

		#Euclidean embeddings
		X1.mds <- sammon(X1, k=M)
		X1.mds$dist <- matrix(0, Nall, Nall)
		for(i in 1:Nall){
			for(j in 1:Nall){
				X1.mds$dist[i,j] <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[j,])**2))
			}
		}

		#MSE
		MSE_H1[which(Ms==M)] <- sum((c(D[upper.tri(D)])-c(X1.hydra$dist[upper.tri(D)]))**2) / choose(Nall,2)
		MSE_H2[which(Ms==M)] <- sum((c(D[upper.tri(D)])-c(log(cosh(X2.hydra$dist[upper.tri(D)]))))**2) / choose(Nall,2)
		MSE_MDS[which(Ms==M)] <- sum((c(D[upper.tri(D)])-c(X1.mds$dist[upper.tri(D)]))**2) / choose(Nall,2)
	}

	output <- rbind(MSE_MDS, MSE_H1, MSE_H2)
	rownames(output) <- c("MDS","H1","H2")
	colnames(output) <- Ms
	write.csv(output, paste("out/MSE_Dall_",nt,".csv",sep=""))
}


#merge result
output <- matrix(0, 3, length(Ms))
rownames(output) <- c("MDS","H1","H2")
colnames(output) <- Ms
for(i in 1:Ntrees){
	tmp <- read.csv(paste("out/MSE_Dall_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/MSE_Dall_all.csv")
