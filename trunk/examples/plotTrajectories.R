require(ape)

tree <- read.tree('VolzSIRgamma_truth.tree')
treeHeight <- 44.4188604629 # Used phylostat. Why is this non-trivial to obtain in APE!?


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

deltaT <- c()
for (idx in 1:length(trajectories)) {
    deltaT[idx] <- trajectories[[idx]]$t[1]-(t0-treeHeight)
}


# Generate figure
pdf('TreeAndTraj.pdf', width=7, height=5)
par(mar=c(3,3,1,1))
par(mgp=c(2,0.5,0))

plot(tree, show.tip.label=F, y.lim=c(-1000,1000), x.lim=c(0,44))
title(xlab='Time since MRCA')

box()
for (idx in 1:length(trajectories)) {
    lines(trajectories[[idx]]$t - (t0-treeHeight), 2*(trajectories[[idx]]$I - 500), col=rgb(0,0,0,.2))
}
lines(c(treeHeight, treeHeight), c(-1000,1000), col='red', lty=2)
axis(side=1)

dev.off()
