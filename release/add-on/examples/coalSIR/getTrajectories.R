getTrajectories <- function(df) {
    nTraj <- dim(df)[1]
    trajectories <- list()

    for (idx in 1:nTraj) {
        
        Istr <- strsplit(df$trajI[[idx]],",")[[1]]
        Sstr <- strsplit(df$trajS[[idx]],",")[[1]]
        tstr <- strsplit(df$trajTime[[idx]],",")[[1]]
        I <- rep(0, length(Istr))
        S <- rep(0, length(Istr))
        time <- rep(0, length(tstr))
        
        for (tidx in 1:length(Istr)) {
            I[tidx] <- as.numeric(Istr[tidx])
            S[tidx] <- as.numeric(Sstr[tidx])
            time[tidx] <- as.numeric(tstr[tidx])
        }
        trajectories[[idx]] <- list()
        trajectories[[idx]]$I <- I
        trajectories[[idx]]$S <- S
        trajectories[[idx]]$t <- time

    }

    return (trajectories)
}

df <- read.table("stochCoalTreeMCMC.traj", header=T,
                 colClasses=c("numeric","numeric","numeric","character","character", "character"))
traj <- getTrajectories(df)

for (i in seq(1:length(traj))) {
    if (i==1) {
        plot(traj[[1]]$t, traj[[1]]$I/(2*0.00075*traj[[1]]$S), 'l', ylim=c(0,1000))
    } else {
        lines(traj[[i]]$t, traj[[i]]$I/(2*0.00075*traj[[i]]$S), col=rgb(0,0,0,0.2))
    }
}
