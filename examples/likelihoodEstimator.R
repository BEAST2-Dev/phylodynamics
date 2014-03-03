# Monte Carlo estimate of SIR parameter likelihoods given tree.
# This script requires PEA (github.com/tgvaughan/PEA).

require(pea)

simSIRTraj <- function(beta, gamma, S0, T) {

    S <- S0
    I <- 1
    t <- 0

    idx <- 1
    while (TRUE) {

        ainfect <- beta*S[idx]*I[idx]
        aremove <- gamma*I[idx]
        atot <- ainfect + aremove

        if (atot==0) {
            tnext <- Inf
        } else {
            tnext <- t[idx] + rexp(1, atot)
        }

        if (tnext>T)
            break

        t[idx+1] <- tnext

        if (runif(1,min=0,max=atot)<ainfect) {
            S[idx+1] <- S[idx] - 1
            I[idx+1] <- I[idx] + 1
        } else {
            S[idx+1] <- S[idx]
            I[idx+1] <- I[idx] - 1
        }

        idx <- idx + 1
    }

    res <- list()
    res$S <- S
    res$I <- I
    res$t <- t

    return(res)
}

simSIRTrajCond <- function(beta, gamma, S0, T, k) {

    while (TRUE) {
        res <- simSIRTraj(beta, gamma, S0, T)
        if (tail(res$I,n=1)>=k)
            return(res)
    }
}

simSIRTrajTL <- function(beta, gamma, S0, T, steps=1000) {

    dt <- T/(steps-1)
    
    t <- seq(0, T, length.out=steps)
    S <- S0
    I <- 1
    
    for (tidx in 2:steps) {
        ainfect <- beta*S[tidx-1]*I[tidx-1]
        aremove <- gamma*I[tidx-1]

        ninfect <- rpois(1, ainfect*dt)
        nremove <- rpois(1, aremove*dt)

        S[tidx] <- S[tidx-1] - ninfect
        I[tidx] <- I[tidx-1] + ninfect - nremove

        if (S[tidx]<0)
            S[tidx] <- 0 # hack to avoid -ve pops due to finite time-step errors
        
        if (I[tidx]<0)
            I[tidx] <- 0 # hack to avoid -ve pops due to finite time-step errors
    }

    res <- list()
    res$t <- t
    res$S <- S
    res$I <- I

    return(res)
}

simSIRTrajCondTL <- function(beta, gamma, S0, T, k, steps=1000) {

    while (TRUE) {
        res <- simSIRTrajTL(beta, gamma, S0, T, steps)
        if (tail(res$I,n=1)>=k)
            return(res)
    }
}


# Likelihood estimate

getCoalescentTreeDensity <- function(tree, beta, gamma, S0, origin, Ntraj=1000) {

    logDensity <- 0
    logDensitySq <- 0

    treeEvents <- getNodeHeights(tree)

    for (trajIdx in 1:Ntraj) {
        cat(paste("beta:",beta,"gamma:",gamma,"S0:",S0,"origin:",origin,"Trajectory",trajIdx,"of",Ntraj,"\n"))

        thisLogDensity <- 0

        traj <- simSIRTrajCond(beta, gamma, S0, origin, 1)

        Svec <- rev(traj$S)
        Ivec <- rev(traj$I)
        tvec <- origin - rev(traj$t)
        
        trajIdx <- 1
        t <- 0

        for (idx in 1:length(treeEvents$heights)) {
            while (treeEvents$heights[idx]>tvec[trajIdx]) {
                rate <- choose(treeEvents$lineages[idx],2)*2*beta*Svec[trajIdx]/Ivec[trajIdx]
                thisLogDensity <- thisLogDensity + -(tvec[trajIdx]-t)*rate
                t <- tvec[trajIdx]
                trajIdx <- trajIdx + 1
            }

            rate <- choose(treeEvents$lineages[idx],2)*2*beta*Svec[trajIdx]/Ivec[trajIdx]
            thisLogDensity <- thisLogDensity + -(treeEvents$heights[idx]-t)*rate
            thisLogDensity <- thisLogDensity + log(2*beta*Svec[trajIdx]/Ivec[trajIdx])

            t <- treeEvents$heights[idx]
        }

        logDensity <- logDensity + thisLogDensity
        logDensitySq <- logDensitySq + thisLogDensity^2
    }

    res <- list()
    res$mean <- logDensity/Ntraj
    res$sd <- sqrt(logDensitySq/Ntraj - res$mean^2)
    res$SEM <- res$sd/sqrt(Ntraj)

    return (res)
}


# Estimate likelihood from test tree 

gamma <- 0.3
beta <- 0.00075
S0 <- 999
origin <- 47.1834399026

tree <- read.tree('VolzSIRgamma_truth.tree')

gammaVec <- seq(.1,.5,by=.05)
ll <- rep(0, length(gammaVec))
llsd <- rep(0, length(gammaVec))
llSEM <- rep(0, length(gammaVec))

Ntraj <- 100

for (i in 1:length(gammaVec)) {
    res <- getCoalescentTreeDensity(tree, beta, gammaVec[i], S0, origin, Ntraj)
    ll[i] <- res$mean
    llsd[i] <- res$sd
    llSEM[i] <- res$SEM
}


# Create figure
pdf('gammaLikelihoodFromR.pdf', width=7, height=5)

plot(gammaVec, ll, 'o',
     xlab=expression(gamma),
     ylab='Log likelihood',
     main=paste('Log likelihoods from simulated tree (',Ntraj,' trajectories)',sep=''))
lines(c(0.3,0.3), c(-1e10,1e10), lty=2, col='blue', lwd=2)
lines(gammaVec, ll+llSEM*2, lty=2)
lines(gammaVec, ll-llSEM*2, lty=2)
legend('bottomleft', inset=.05, c('Truth', '+/- 2*SEM'), lty=2, lwd=c(2,1), col=c('blue', 'black'))

dev.off()
