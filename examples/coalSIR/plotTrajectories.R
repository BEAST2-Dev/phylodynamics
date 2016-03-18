require(pea)

tree <- read.tree('VolzSIRgamma_truth.tree')
treeHeight <- max(getNodeHeights(tree)$heights)
nLeaves <- length(getLeaves(tree))

getTrajectories <- function(df) {
    nTraj <- dim(df)[1]
    trajectories <- list()

    for (idx in 1:nTraj) {
        
        Istr <- strsplit(df$trajI[[idx]],",")[[1]]
        tstr <- strsplit(df$trajTime[[idx]],",")[[1]]
        I <- rep(0, length(Istr))
        time <- rep(0, length(tstr))
        
        for (tidx in 1:length(Istr)) {
            I[tidx] <- as.numeric(Istr[tidx])
            time[tidx] <- as.numeric(tstr[tidx])
        }
        trajectories[[idx]] <- list()
        trajectories[[idx]]$I <- I
        trajectories[[idx]]$t <- time

    }

    return (trajectories)
}

df <- read.table('VolzSIRgamma_saveTrajectories.log', header=T, as.is=T)
t0 <- df$volz.origin[[1]]
trajectories <- getTrajectories(df)

Iend <- c()
Imax <- 0
for (idx in 1:length(trajectories)) {
    n <- length(trajectories[[idx]]$I)
    Iend[idx] <- trajectories[[idx]]$I[n]
    Imax <- max(Imax, max(trajectories[[idx]]$I))
}



# Generate figure
pdf('TreeAndTraj.pdf', width=7, height=5)
par(mar=c(3,3,1,1))
par(mgp=c(2,0.5,0))

layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(1,.5))

plot(tree, show.tip.label=F, y.lim=c(-nLeaves,nLeaves), x.lim=c(0,ceiling(treeHeight)))
title(xlab='Time since MRCA')

box()
nTraj <- min(length(trajectories), 500)
for (idx in 1:nTraj) {
    lines(trajectories[[idx]]$t - (t0-treeHeight), trajectories[[idx]]$I*nLeaves/Imax - nLeaves, col=rgb(0,0,0,.1))
}
lines(c(treeHeight, treeHeight), c(-nLeaves,nLeaves), col='red', lty=2)
axis(side=1)

#hist(Iend, breaks=seq(-.5,1000.5,by=1), xlim=c(0,5.5), prob=T,
hist(Iend, prob=T,
     xlab="Final count", ylab="Probability", main="",
     col="green")

dev.off()
