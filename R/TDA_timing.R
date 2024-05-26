library(TDAvec)

tau <- 0.3 # parameter in [0,1] which controls the size of blocks around each point of the diagram 

totTimes <- c()
genTimes <- c()
vectTimes <- c()
times <- data.frame()
nSim <- 50
for(N in c(10, 50, 100, 200, 250, 300, 350, 400)) {
  cat("N=", N, "\n")
  for(i in 1:nSim) {
    set.seed(nSim*n + i)
    cat(".")
    start_totTime <- Sys.time()
    start_genTime <- Sys.time()
    X <- TDA::circleUnif(N, r = 1) + rnorm(2*N,mean = 0, sd = 0.2)
    D <- TDA::ripsDiag(X,maxdimension = 1,maxscale = 2)$diagram 
    DD <- D
    # switch from the birth-death to the birth-persistence coordinates
    D[,3] <- D[,3] - D[,2] 
    colnames(D)[3] <- "Persistence"
    
    # construct one-dimensional grid of scale values
    ySeqH0 <- unique(quantile(D[D[,1] == 0,3],probs = seq(0,1, by = 0.2))) 
    # compute VPB for homological dimension H_0
    xSeqH1 <- unique(quantile(D[D[,1] == 1,2],probs = seq(0,1,by = 0.2)))
    ySeqH1 <- unique(quantile(D[D[,1] == 1,3],probs = seq(0,1,by = 0.2)))
    genTime_ <- as.numeric(difftime(Sys.time(), start_genTime, units = "secs"))
    
    start_vectTime <- Sys.time()
    vpb0 <- computeVPB(D,homDim = 0,xSeq = NA, ySeqH0, tau)
    vpb1 <- computeVPB(D,homDim = 1,xSeqH1,ySeqH1,tau) 
    vectTime_ <- as.numeric(difftime(Sys.time(), start_vectTime, units = "secs"))
  
    totTime_ <- as.numeric(difftime(Sys.time(), start_totTime, units = "secs"))
    
    times <- rbind(times, data.frame(N = N, totTime = totTime_, genTime = genTime_, vectTime = vectTime_))
  }
  cat("\n")
  write.csv(times, "./times.csv", row.names = FALSE)
}


suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

times %>% tidyr::pivot_longer(cols = -N, names_to = "time", values_to = "duration") %>% group_by(N, time) %>%
  summarize(mean = mean(duration), sd=sd(duration)) %>%
  ggplot(aes(x=N, y=mean, color = time)) + geom_point() + geom_line() + 
  scale_y_log10()
