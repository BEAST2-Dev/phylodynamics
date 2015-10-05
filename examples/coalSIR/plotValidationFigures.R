# Plot tree height sample densities obtained via simulation under the
# stochastic and deterministic coalescent models against densities
# obtained via MCMC as implemented by the phylodynamics package.
#
# Each of the four BEAST scripts in this directory must be run before
# this script is executed.

removeBurnin <- function(df, burninFrac=0.1) {
    N <- dim(df)[1]
    return(df[-(1:ceiling(burninFrac*N)),])
}

# Deterministic coalescent
dfdetsim <- read.table('detCoalTreeSim.treestat', sep='\t', header=T, colClasses="numeric")
dfdetmcmc <- removeBurnin(read.table('detCoalTreeMCMC.treestat', sep='\t', header=T, colClasses="numeric"))

pdf('detCoalFigure.pdf', width=3, height=8)
par(mfcol=c(3,1))
par(mgp=c(2,0.5,0), mar=c(5,3,2,1))

hdetsim <- hist(dfdetsim$Tree.Height, breaks=seq(19,31,length.out=25), plot=F)
hdetmcmc <- hist(dfdetmcmc$Tree.Height, breaks=seq(19,31,length.out=25), plot=F)

plot(hdetmcmc$mids, hdetmcmc$density, 'o',
     xlim=c(24,31), ylim=c(0,0.4),lwd=2, col='blue',
     xlab='Transmission tree height',
     ylab='Relative frequency',
     main='(a)')
lines(hdetsim$mids, hdetsim$density, 'o', lwd=2, col='red')

legend('topleft', inset=.05, c('MCMC', 'Direct simulation'),
       lty=1, pch=1, lwd=2, col=c('blue','red'),
       cex=0.8)

hdetsim <- hist(dfdetsim$Tree.Length, breaks=seq(180,300,by=10), plot=F)
hdetmcmc <- hist(dfdetmcmc$Tree.Length, breaks=seq(180,300,by=10), plot=F)

plot(hdetmcmc$mids, hdetmcmc$density, 'o',
     lwd=2, col='blue',
     xlab='Transmission tree length',
     ylab='Relative frequency',
     main='(b)')
lines(hdetsim$mids, hdetsim$density, 'o', lwd=2, col='red')

hdetsim <- hist(dfdetsim$Cherry.count, breaks=seq(2.5,10.5,by=1), plot=F)
hdetmcmc <- hist(dfdetmcmc$Cherry.count, breaks=seq(2.5,10.5,by=1), plot=F)

plot(hdetmcmc$mids, hdetmcmc$density, 'o',
     lwd=2, col='blue',
     xlab='Cherry count',
     ylab='Relative frequency',
     main='(c)')
lines(hdetsim$mids, hdetsim$density, 'o', lwd=2, col='red')

dev.off()

# Stochastic coalescent
dfstochsim <- read.table('stochCoalTreeSim.treestat', sep='\t', header=T, colClasses="numeric")
dfstochmcmc <- removeBurnin(read.table('stochCoalTreeMCMC.treestat', sep='\t', header=T, colClasses="numeric"))

pdf('stochCoalFigure.pdf', width=3, height=8)
par(mfcol=c(3,1))
par(mgp=c(2,0.5,0), mar=c(5,3,2,1))

hstochsim <- hist(dfstochsim$Tree.Height, breaks=seq(19,31,length.out=25), plot=F)
hstochmcmc <- hist(dfstochmcmc$Tree.Height, breaks=seq(19,31,length.out=25), plot=F)

plot(hstochmcmc$mids, hstochmcmc$density, 'o',
     xlim=c(19,31), ylim=c(0,0.3), lwd=2, col='blue',
     xlab='Transmission tree height',
     ylab='Relative frequency',
     main='(a)')
lines(hstochsim$mids, hstochsim$density, 'o', lwd=2, col='red')

legend('topleft', inset=.05, c('MCMC', 'Direct simulation'), lty=1, pch=1, lwd=2, col=c('blue','red'))

hstochsim <- hist(dfstochsim$Tree.Length, breaks=seq(50,400,by=25), plot=F)
hstochmcmc <- hist(dfstochmcmc$Tree.Length, breaks=seq(50,400,by=25), plot=F)

plot(hstochmcmc$mids, hstochmcmc$density, 'o',
     lwd=2, col='blue',
     xlab='Transmission tree length',
     ylab='Relative frequency',
     main='(b)')
lines(hstochsim$mids, hstochsim$density, 'o', lwd=2, col='red')

hstochsim <- hist(dfstochsim$Cherry.count, breaks=seq(2.5,10.5,by=1), plot=F)
hstochmcmc <- hist(dfstochmcmc$Cherry.count, breaks=seq(2.5,10.5,by=1), plot=F)

plot(hstochmcmc$mids, hstochmcmc$density, 'o',
     lwd=2, col='blue',
     xlab='Cherry count',
     ylab='Relative frequency',
     main='(c)')
lines(hstochsim$mids, hstochsim$density, 'o', lwd=2, col='red')

dev.off()
