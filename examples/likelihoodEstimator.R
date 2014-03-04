# Monte Carlo estimate of SIR parameter likelihoods given tree.
# This script requires PEA (github.com/tgvaughan/PEA).

require(pea)

# Integrate the SIR ODEs over the specified domain
deterministicSIR <- function(beta, gamma, S0, T, steps, maxIter=3) {

    S <- S0
    I <- 1
    t <- 0

    dt <- T/(steps-1)

    for (tidx in 2:steps) {

        Sp <- S[tidx-1]
        Ip <- I[tidx-1]

        for (iter in 1:maxIter) {
            dSdt <- -beta*Sp*Ip
            dIdt <- beta*Sp*Ip - gamma*Ip

            Sp <- S[tidx-1] + 0.5*dt*dSdt
            Ip <- I[tidx-1] + 0.5*dt*dIdt
        }

        S[tidx] <- 2*Sp - S[tidx-1]
        I[tidx] <- 2*Ip - I[tidx-1]
        t[tidx] <- tidx*dt

    }

    res <- list()
    res$S <- S
    res$I <- I
    res$t <- t

    return(res)
}


# Likelihood estimate from the deterministic solution
getDeterministicCoalescentTreeDensity <- function(tree, beta, gamma, S0, origin, steps=1000) {

    treeEvents <- getNodeHeights(tree)

    logDensity <- 0

    traj <- deterministicSIR(beta, gamma, S0, origin, steps)

    Svec <- rev(traj$S)
    Ivec <- rev(traj$I)
    tvec <- origin - rev(traj$t)
        
    tidx <- 1
    t <- 0

    for (idx in 2:length(treeEvents$heights)) {
        while (treeEvents$heights[idx]>tvec[tidx]) {
                
            rate <- 2*beta*Svec[tidx]/Ivec[tidx]

            # Waiting time contribution
            logDensity <- logDensity + -(tvec[tidx]-t)*choose(treeEvents$lineages[idx-1],2)*rate
                
            t <- tvec[tidx]
            tidx <- tidx + 1
        }

        rate <- 2*beta*Svec[tidx]/Ivec[tidx]

        # Waiting time contribution
        logDensity <- logDensity +
            -(treeEvents$heights[idx]-t)*choose(treeEvents$lineages[idx-1],2)*rate
        
        # Event time contribution (only if coalescence)
        if (treeEvents$lineages[idx]<treeEvents$lineages[idx-1])
            logDensity <- logDensity + log(rate)
        
        t <- treeEvents$heights[idx]
    }

    return (logDensity)
}


# Stochastic simulation of an SIR trajectory
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


# Stochastic likelihood estimate using a number of simulated SIR epidemics
getCoalescentTreeDensity <- function(tree, beta, gamma, S0, origin, Ntraj=1000) {

    logDensity <- rep(0,Ntraj)

    treeEvents <- getNodeHeights(tree)

    trajIdx <- 0
    goodTrajIdx <- 0
    while (goodTrajIdx<Ntraj) {
        trajIdx <- trajIdx + 1
        if (trajIdx%%100 == 0)
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



###### MAIN ######

gamma <- 0.3
beta <- 0.00075
S0 <- 999
origin <- 12.7808530307

tree <- read.tree('VolzSIRgamma_truth.tree')

# Estimate STOCHASTIC coalescent likelihoods for different gammas

gammaVec <- seq(.1,.7,by=.05)
llmean <- rep(0, length(gammaVec))
llres <- list()

Ntraj <- 1000

for (i in 1:length(gammaVec)) {
    res <- getCoalescentTreeDensity(tree, beta, gammaVec[i], S0, origin, Ntraj)
    llmean[i] <- res$mean
    llres[[i]] <- res
}

# Estimate DETERMINISTIC coalescent likelihoods for different gammas

gammaVecDet <- seq(0.1,0.35,by=0.01)
lldet <- rep(0, length(gammaVecDet))
for (i in 1:length(gammaVecDet)) {
    lldet[i] <- getDeterministicCoalescentTreeDensity(tree, beta, gammaVecDet[i], S0, origin)
}

# Load in Java code results for same tree:
df10000 <- read.table('alex_likelihood10000.txt')
df1000 <- read.table('alex_likelihood1000.txt')
df100 <- read.table('alex_likelihood100.txt')

# Create figure
pdf('gammaLikelihoodFromR.pdf', width=7, height=5)

plot(gammaVec, llmean, 'o', ylim=c(-440,-400),
     xlab=expression(gamma),
     ylab='Log likelihood',
     main=paste('Log likelihoods from simulated tree (',Ntraj,' full trajectories)',sep=''),
     col='blue')
lines(df10000[[1]], df10000[[3]], 'o', col='red')
lines(gammaVecDet, lldet, 'o', col='purple')
             
lines(c(0.3,0.3), c(-1e10,1e10), lty=2, col='grey', lwd=2)
legend('bottomright', inset=.05, c('R','Java (10000)','R (det.)', 'Truth'), lty=c(1,1,1,2), pch=c(1,1,1,NA), lwd=c(1,1,1,2), col=c('blue','red','purple','grey'))

dev.off()
