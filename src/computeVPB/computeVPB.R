library(TDAvec)

# creating dataset
N <- 100 
set.seed(123)
X <- TDA::circleUnif(N,r=1) + rnorm(2*N,mean = 0,sd = 0.2)
write.csv(X, "./unitCircle.csv", row.names = FALSE)

# compute a persistence diagram using the Rips filtration built on top of X
D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 
DD <- D
# switch from the birth-death to the birth-persistence coordinates
D[,3] <- D[,3] - D[,2] 
colnames(D)[3] <- "Persistence"
write.csv(D, "./PD.csv", row.names = FALSE)

# construct one-dimensional grid of scale values
ySeqH0 <- unique(quantile(D[D[,1]==0,3],probs = seq(0,1,by=0.2))) 
tau <- 0.3 # parameter in [0,1] which controls the size of blocks around each point of the diagram 
# compute VPB for homological dimension H_0
vpb.0 <- computeVPB(D,homDim = 0,xSeq=NA,ySeqH0,tau)
write.csv(vpb.0, "vpb0.csv", row.names = FALSE)

xSeqH1 <- unique(quantile(D[D[,1]==1,2],probs = seq(0,1,by=0.2)))
ySeqH1 <- unique(quantile(D[D[,1]==1,3],probs = seq(0,1,by=0.2)))
# compute VPB for homological dimension H_1
vpb.1 <- computeVPB(D,homDim = 1,xSeqH1,ySeqH1,tau)
write.csv(vpb.1, "vpb1.csv", row.names = FALSE)
