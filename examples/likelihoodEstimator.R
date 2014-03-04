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


# Likelihood estimate

getCoalescentTreeDensity <- function(tree, beta, gamma, S0, origin, Ntraj=1000) {

    logDensity <- rep(0,Ntraj)

    treeEvents <- getNodeHeights(tree)

    trajIdx <- 0
    goodTrajIdx <- 0
    while (goodTrajIdx<Ntraj) {
        trajIdx <- trajIdx + 1
        cat(paste("beta:",beta,"gamma:",gamma,"S0:",S0,"origin:",origin,"Trajectory",goodTrajIdx,"of",Ntraj,"\n"))

        thisLogDensity <- 0

        traj <- simSIRTraj(beta, gamma, S0, origin)

        # Check for trajectories shorter than tree
        if (traj$I[length(traj$I)]==0) {
            logDensity[trajIdx] <- -Inf
            next
        } else {
            goodTrajIdx <- goodTrajIdx + 1
        }

        Svec <- rev(traj$S)
        Ivec <- rev(traj$I)
        tvec <- origin - rev(traj$t)
        
        tidx <- 1
        t <- 0

        for (idx in 2:length(treeEvents$heights)) {
            while (treeEvents$heights[idx]>tvec[tidx]) {
                
                rate <- 2*beta*Svec[tidx]/Ivec[tidx]

                # Waiting time contribution
                thisLogDensity <- thisLogDensity + -(tvec[tidx]-t)*choose(treeEvents$lineages[idx-1],2)*rate
                
                t <- tvec[tidx]
                tidx <- tidx + 1
            }

            rate <- 2*beta*Svec[tidx]/Ivec[tidx]

            # Waiting time contribution
            thisLogDensity <- thisLogDensity +
                -(treeEvents$heights[idx]-t)*choose(treeEvents$lineages[idx-1],2)*rate

            # Event time contribution (only if coalescence)
            if (treeEvents$lineages[idx]<treeEvents$lineages[idx-1])
                thisLogDensity <- thisLogDensity + log(rate)

            t <- treeEvents$heights[idx]
        }

        logDensity[trajIdx] <- thisLogDensity
    }

    # Calculate log of mean of densities
    
    maxLogDensity <- max(logDensity)

    logDensityShifted <- logDensity - maxLogDensity
    scaledDensities <- exp(logDensityShifted)
    meanScaledDensity <- mean(scaledDensities)
    quantiles <- quantile(scaledDensities, probs=c(0.025, 0.975))

    res <- list()
    res$mean <- log(meanScaledDensity) + maxLogDensity
    res$lower <- log(quantiles[1]) + maxLogDensity
    res$upper <- log(quantiles[1]) + maxLogDensity    
    res$logDensity <- logDensity

    return (res)
}


# Estimate likelihood from test tree 

gamma <- 0.3
beta <- 0.00075
S0 <- 999
origin <- 12.7808530307

tree <- read.tree('VolzSIRgamma_truth.tree')

gammaVec <- seq(.1,.7,by=.05)
llmean <- rep(0, length(gammaVec))
llres <- list()

Ntraj <- 1000

for (i in 1:length(gammaVec)) {
    res <- getCoalescentTreeDensity(tree, beta, gammaVec[i], S0, origin, Ntraj)
    llmean[i] <- res$mean
    llres[[i]] <- res
}

# Load in Java code results for same tree:
df10000 <- read.table('alex_likelihood10000.txt')
df1000 <- read.table('alex_likelihood1000.txt')
df100 <- read.table('alex_likelihood100.txt')

# Create figure
pdf('gammaLikelihoodFromR.pdf', width=7, height=5)

plot(gammaVec, llmean, 'o',
     xlab=expression(gamma),
     ylab='Log likelihood',
     main=paste('Log likelihoods from simulated tree (',Ntraj,' full trajectories)',sep=''),
     col='blue')
lines(df10000[[1]], df10000[[3]], 'o', col='red')
#lines(df1000[[1]], df1000[[3]], 'o', col='red')
#lines(df100[[1]], df100[[3]], 'o', col='red')
lines(c(0.3,0.3), c(-1e10,1e10), lty=2, col='grey', lwd=2)
#lines(gammaVec, ll+llSEM*2, lty=2)
#lines(gammaVec, ll-llSEM*2, lty=2)
legend('bottomright', inset=.05, c('R','Java (10000)','Truth'), lty=c(1,1,2), pch=c(1,1,NA), lwd=c(1,1,2), col=c('blue','red','grey'))

dev.off()
