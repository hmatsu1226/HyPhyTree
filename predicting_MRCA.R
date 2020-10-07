library(hydra)
library(ape)
library(MASS)
library(ROCR)

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

for(nt in 1:Ntrees){
	tree <- read.tree(file=paste("data/tree_",nt,".txt",sep=""))
	tree$mrca <- mrca(tree)
	D <- dist.nodes(tree)

	Nall <- dim(D)[1]
	X1 <- D
	X2 <- acosh(exp(D))

	AUC_H1 <- rep(0, length(Ms))
	AUC_H2 <- rep(0, length(Ms))
	AUC_MDS <- rep(0, length(Ms))

	for(M in Ms){
		#general hyperbolic embeddings
		X1.hydra <- hydraPlus(X1, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z1 <- X1.hydra$r * X1.hydra$directional

		#our hyperbolic embeddings
		X2.hydra <- hydraPlus(X2, dim=M, curvature=1, alpha=1, equi.adj=0, control=list(return.dist=1, isotropic.adj=FALSE))
		Z2 <- X2.hydra$r * X2.hydra$directional

		#Euclidean embeddings
		X1.mds <- sammon(X1, k=M)

		itenum <- 1000
		angle_H1_pos <- rep(0,itenum)
		angle_H1_neg <- rep(0,itenum)
		angle_H2_pos <- rep(0,itenum)
		angle_H2_neg <- rep(0,itenum)
		angle_MDS_pos <- rep(0,itenum)
		angle_MDS_neg <- rep(0,itenum)
		for(ite in 1:itenum){
			tmp <- sample(1:(Nleaves-1))
			i <- tmp[1]
			j <- tmp[2]
			k <- tree$mrca[i,j]

			#H1
			a <- Poincare_dist(Z1[i,], Z1[k,])
			b <- Poincare_dist(Z1[j,], Z1[k,])
			c <- Poincare_dist(Z1[i,], Z1[j,])
			tmp1 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z1[i,], Z1[k,])
			b <- Poincare_dist(Z1[o,], Z1[k,])
			c <- Poincare_dist(Z1[i,], Z1[o,])
			tmp2 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z1[o,], Z1[k,])
			b <- Poincare_dist(Z1[j,], Z1[k,])
			c <- Poincare_dist(Z1[o,], Z1[j,])
			tmp3 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			angle_H1_pos[ite] <- min(tmp1,tmp2,tmp3)

			#H2
			a <- Poincare_dist(Z2[i,], Z2[k,])
			b <- Poincare_dist(Z2[j,], Z2[k,])
			c <- Poincare_dist(Z2[i,], Z2[j,])
			tmp1 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z2[i,], Z2[k,])
			b <- Poincare_dist(Z2[o,], Z2[k,])
			c <- Poincare_dist(Z2[i,], Z2[o,])
			tmp2 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z2[o,], Z2[k,])
			b <- Poincare_dist(Z2[j,], Z2[k,])
			c <- Poincare_dist(Z2[o,], Z2[j,])
			tmp3 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			angle_H2_pos[ite] <- min(tmp1,tmp2,tmp3)

			#MDS
			a <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[k,])**2))
			b <- sqrt(sum((X1.mds$points[j,]-X1.mds$points[k,])**2))
			c <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[j,])**2))
			tmp1 <- acos((a**2+b**2-c**2)/(2*a*b))

			a <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[k,])**2))
			b <- sqrt(sum((X1.mds$points[o,]-X1.mds$points[k,])**2))
			c <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[o,])**2))
			tmp2 <- acos((a**2+b**2-c**2)/(2*a*b))

			a <- sqrt(sum((X1.mds$points[o,]-X1.mds$points[k,])**2))
			b <- sqrt(sum((X1.mds$points[j,]-X1.mds$points[k,])**2))
			c <- sqrt(sum((X1.mds$points[o,]-X1.mds$points[j,])**2))
			tmp3 <- acos((a**2+b**2-c**2)/(2*a*b))

			angle_MDS_pos[ite] <- min(tmp1,tmp2,tmp3)


			#negative set
			while(1){
				if(runif(1) < 0.5){
					tmp <- sample(1:(Nleaves-1))
					if(i == tmp[1]){
						next
					}

					tmp <- tree$mrca[i,tmp[1]]
					if(k != tmp){
						k <- tmp
						break
					}
				}else{
					tmp <- sample(1:(Nleaves-1))
					if(j == tmp[1]){
						next
					}

					tmp <- tree$mrca[j,tmp[1]]
					if(k != tmp){
						k <- tmp
						break
					}
				}
			}

			#H1
			a <- Poincare_dist(Z1[i,], Z1[k,])
			b <- Poincare_dist(Z1[j,], Z1[k,])
			c <- Poincare_dist(Z1[i,], Z1[j,])
			tmp1 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z1[i,], Z1[k,])
			b <- Poincare_dist(Z1[o,], Z1[k,])
			c <- Poincare_dist(Z1[i,], Z1[o,])
			tmp2 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z1[o,], Z1[k,])
			b <- Poincare_dist(Z1[j,], Z1[k,])
			c <- Poincare_dist(Z1[o,], Z1[j,])
			tmp3 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			angle_H1_neg[ite] <- min(tmp1,tmp2,tmp3)

			#H2
			a <- Poincare_dist(Z2[i,], Z2[k,])
			b <- Poincare_dist(Z2[j,], Z2[k,])
			c <- Poincare_dist(Z2[i,], Z2[j,])
			tmp1 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z2[i,], Z2[k,])
			b <- Poincare_dist(Z2[o,], Z2[k,])
			c <- Poincare_dist(Z2[i,], Z2[o,])
			tmp2 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			a <- Poincare_dist(Z2[o,], Z2[k,])
			b <- Poincare_dist(Z2[j,], Z2[k,])
			c <- Poincare_dist(Z2[o,], Z2[j,])
			tmp3 <- acos((cosh(a)*cosh(b)-cosh(c))/(sinh(a)*sinh(b)))

			angle_H2_neg[ite] <- min(tmp1,tmp2,tmp3)

			#MDS
			a <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[k,])**2))
			b <- sqrt(sum((X1.mds$points[j,]-X1.mds$points[k,])**2))
			c <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[j,])**2))
			tmp1 <- acos((a**2+b**2-c**2)/(2*a*b))

			a <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[k,])**2))
			b <- sqrt(sum((X1.mds$points[o,]-X1.mds$points[k,])**2))
			c <- sqrt(sum((X1.mds$points[i,]-X1.mds$points[o,])**2))
			tmp2 <- acos((a**2+b**2-c**2)/(2*a*b))

			a <- sqrt(sum((X1.mds$points[o,]-X1.mds$points[k,])**2))
			b <- sqrt(sum((X1.mds$points[j,]-X1.mds$points[k,])**2))
			c <- sqrt(sum((X1.mds$points[o,]-X1.mds$points[j,])**2))
			tmp3 <- acos((a**2+b**2-c**2)/(2*a*b))

			angle_MDS_neg[ite] <- min(tmp1,tmp2,tmp3)
		}

		#ROC
		rocdata_H1 <- matrix(0,itenum*2,2)
		rocdata_H2 <- matrix(0,itenum*2,2)
		rocdata_MDS <- matrix(0,itenum*2,2)
		for(ite in 1:itenum){
			rocdata_H1[ite*2-1,1] <- angle_H1_pos[ite]
			rocdata_H1[ite*2-1,2] <- 1
			rocdata_H1[ite*2,1] <- angle_H1_neg[ite]
			rocdata_H1[ite*2,2] <- 0


			rocdata_H2[ite*2-1,1] <- angle_H2_pos[ite]
			rocdata_H2[ite*2-1,2] <- 1
			rocdata_H2[ite*2,1] <- angle_H2_neg[ite]
			rocdata_H2[ite*2,2] <- 0

			rocdata_MDS[ite*2-1,1] <- angle_MDS_pos[ite]
			rocdata_MDS[ite*2-1,2] <- 1
			rocdata_MDS[ite*2,1] <- angle_MDS_neg[ite]
			rocdata_MDS[ite*2,2] <- 0
		}

		pred_H1 <- prediction(rocdata_H1[,1], rocdata_H1[,2])
		perf_H1 <- performance(pred_H1, "tpr", "fpr")
		auc_H1.tmp <- performance(pred_H1,"auc")
		AUC_H1[which(Ms==M)] <- as.numeric(auc_H1.tmp@y.values)

		pred_H2 <- prediction(rocdata_H2[,1], rocdata_H2[,2])
		perf_H2 <- performance(pred_H2, "tpr", "fpr")
		auc_H2.tmp <- performance(pred_H2,"auc")
		AUC_H2[which(Ms==M)] <- as.numeric(auc_H2.tmp@y.values)

		pred_MDS <- prediction(rocdata_MDS[,1], rocdata_MDS[,2])
		perf_MDS <- performance(pred_MDS, "tpr", "fpr")
		auc_MDS.tmp <- performance(pred_MDS,"auc")
		AUC_MDS[which(Ms==M)] <- as.numeric(auc_MDS.tmp@y.values)
		
	}

	output <- rbind(AUC_MDS, AUC_H1, AUC_H2)
	rownames(output) <- c("MDS","H1","H2")
	colnames(output) <- Ms
	write.csv(output, paste("out/AUC_MRCA_",nt,".csv",sep=""))
}


#merge result
output <- matrix(0, 3, length(Ms))
rownames(output) <- c("MDS","H1","H2")
colnames(output) <- Ms
for(i in 1:Ntrees){
	tmp <- read.csv(paste("out/AUC_MRCA_",i,".csv",sep=""), row.names=1)
	output <- output + tmp
}
output <- output/Ntrees
write.csv(output, "out/AUC_MRCA_all.csv")
