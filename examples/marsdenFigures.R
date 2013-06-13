##### Cross-sections

# vary_beta_0 results
df <- read.table('vary_beta_0.txt', header=T)
dfx <- read.table('exactlikbeta_0.txt', header=T)

pdf('vary_beta_0.pdf', onefile=F, width=6, height=5)
plot(df$beta, df$volzSIR-max(df$volzSIR), 'o', col='blue', ylim=c(-300,0),
     xlab=expression(beta),
     ylab='likelihood',
     main='25 taxon tree likelihoods versus beta')
lines(df$beta, df$volz2009-max(df$volz2009[is.finite(df$volz2009)]), 'o', col='red')
lines(df$beta, df$BDSIR-max(df$BDSIR[is.finite(df$BDSIR)]), 'o', col='black')
lines(dfx$beta, dfx$X25-max(dfx$X25), 'o', col='purple')
legend('bottomright', inset=.05,
       c('VolzSIR', 'Volz2009', 'BDSIR','Exact'),
       lty=1, pch=1, col=c('blue','red','black','purple'))
dev.off()

# vary_gamma_0 results
df <- read.table('vary_gamma_0.txt', header=T)
dfx <- read.table('exactlikgamma_0.txt', header=T)

pdf('vary_gamma_0.pdf', onefile=F, width=6, height=5)
plot(df$gamma, df$volzSIR-max(df$volzSIR), 'o', col='blue', ylim=c(-100,0),
     xlab=expression(gamma),
     ylab='likelihood',
     main='25 taxon tree likelihoods versus gamma')
lines(df$gamma, df$volz2009-max(df$volz2009), 'o', col='red')
lines(df$gamma, df$BDSIR-max(df$BDSIR), 'o', col='black')
lines(dfx$gamma, dfx$X25-max(dfx$X25), 'o', col='purple')
legend('bottom', inset=.05,
       c('VolzSIR', 'Volz2009', 'BDSIR','Exact'),
       lty=1, pch=1, col=c('blue','red','black','purple'))
dev.off()

# vary_S0_0 results
df <- read.table('vary_S0_0.txt', header=T)

pdf('vary_S0_0.pdf', onefile=F, width=6, height=5)
plot(df$S0, df$volzSIR-max(df$volzSIR), 'o', col='blue', ylim=c(-150,0),
     xlab=expression(S[0]),
     ylab='likelihood',
     main='25 taxon tree likelihoods versus S0')
lines(df$S0, df$volz2009-max(df$volz2009), 'o', col='red')
lines(df$S0, df$BDSIR-max(df$BDSIR[is.finite(df$BDSIR)]), 'o', col='black')
legend('bottomright', inset=.05,
       c('VolzSIR', 'Volz2009', 'BDSIR','Exact'),
       lty=1, pch=1, col=c('blue','red','black','purple'))
dev.off()


##### Contour plots

beta <- read.table('coord_beta.txt')[[1]]
gamma <- read.table('coord_gamma.txt')[[1]]
S0 <- read.table('coord_S0.txt')[[1]]

# fix_beta_0_VolzSIR and fix_beta_0_VolzSIR2009

z1 <- as.matrix(read.table('fix_beta_0_VolzSIR.txt'))
z2 <- as.matrix(read.table('fix_beta_0_VolzSIR2009.txt'))

pdf('fix_beta_0.pdf', onefile=F, width=6, height=5)
contour(gamma, S0, z1-max(z1[is.finite(z1)]), levels=seq(-10,0,by=.5), col='blue',
        xlab=expression(gamma), ylab=expression(S[0]))
contour(gamma, S0, z2-max(z2[is.finite(z2)]), levels=seq(-10,0,by=.5), add=T, col='red')

# Truth
lines(c(0,1e5),c(1000,1000))
lines(c(2,2),c(0,1e5))

# ML VolzSIR:
zmax <- max(z1[is.finite(z1)])
mlidx <- which(z1==zmax)
ridx <- (mlidx%%100) + 1
cidx <- floor(mlidx/100) + 1
maxgam <- gamma[ridx]
maxS0 <- S0[cidx]
lines(c(0,1e5), c(maxS0,maxS0), col='blue')
lines(c(maxgam,maxgam), c(0,1e5), col='blue')

# ML VolzSIR2009:
zmax <- max(z2[is.finite(z2)])
mlidx <- which(z2==zmax)
ridx <- (mlidx%%100) + 1
cidx <- floor(mlidx/100) + 1
maxgam <- gamma[ridx]
maxS0 <- S0[cidx]
lines(c(0,1e5), c(maxS0,maxS0), col='red')
lines(c(maxgam,maxgam), c(0,1e5), col='red')

# Legend
legend('bottomright', inset=.05, lty=1,
       c('Volz2012', 'Volz2009', 'Truth'),
       col=c('blue','red', 'black'))

dev.off()


# fix_gamma_0_VolzSIR and fix_gamma_0_VolzSIR2009

z1 <- as.matrix(read.table('fix_gamma_0_VolzSIR.txt'))
z2 <- as.matrix(read.table('fix_gamma_0_VolzSIR2009.txt'))

pdf('fix_gamma_0.pdf', onefile=F, width=6, height=5)
contour(S0, beta, z1-max(z1[is.finite(z1)]), levels=seq(-10,0,by=.5), col='blue',
        xlab=expression(S[0]), ylab=expression(beta))
contour(S0, beta, z2-max(z2[is.finite(z2)]), levels=seq(-10,0,by=.5), add=T, col='red')

# Truth
lines(c(0,1e5),c(0.005,0.005))
lines(c(1000,1000),c(0,1e5))

# ML VolzSIR:
zmax <- max(z1[is.finite(z1)])
mlidx <- which(z1==zmax)
ridx <- (mlidx%%100) + 1
cidx <- floor(mlidx/100) + 1
maxS0 <- S0[ridx]
maxbeta <- beta[cidx]
lines(c(maxS0,maxS0), c(0,1e5), col='blue')
lines(c(0,1e5), c(maxbeta,maxbeta), col='blue')

# ML VolzSIR2009:
zmax <- max(z2[is.finite(z2)])
mlidx <- which(z2==zmax)
ridx <- (mlidx%%100) + 1
cidx <- floor(mlidx/100) + 1
maxS0 <- S0[ridx]
maxbeta <- beta[cidx]
lines(c(maxS0,maxS0), c(0, 1e5), col='red')
lines(c(0,1e5), c(maxbeta,maxbeta), col='red')

# Legend
legend('topright', inset=.05, lty=1,
       c('Volz2012', 'Volz2009', 'Truth'),
       col=c('blue','red', 'black'))

dev.off()

# fix_S0_0_VolzSIR and fix_S0_0_VolzSIR2009

z1 <- as.matrix(read.table('fix_S0_0_VolzSIR.txt'))
z2 <- as.matrix(read.table('fix_S0_0_VolzSIR2009.txt'))

pdf('fix_S0_0.pdf', onefile=F, width=6, height=5)
contour(beta, gamma, z1-max(z1[is.finite(z1)]), levels=seq(-10,0,by=.5), col='blue',
        xlab=expression(beta), ylab=expression(gamma))
contour(beta, gamma, z2-max(z2[is.finite(z2)]), levels=seq(-10,0,by=.5), add=T, col='red')

# Truth
lines(c(0.005,0.005),c(0,1e5))
lines(c(0,1e5),c(2,2))

# ML VolzSIR:
zmax <- max(z1[is.finite(z1)])
mlidx <- which(z1==zmax)
ridx <- (mlidx%%100) + 1
cidx <- floor(mlidx/100) + 1
maxbeta <- beta[ridx]
maxgamma <- gamma[cidx]
lines(c(maxbeta,maxbeta), c(0,1e5), col='blue')
lines(c(0,1e5), c(maxgamma,maxgamma), col='blue')

# ML VolzSIR2009:
zmax <- max(z2[is.finite(z2)])
mlidx <- which(z2==zmax)
ridx <- (mlidx%%100) + 1
cidx <- floor(mlidx/100) + 1
maxbeta <- beta[ridx]
maxgamma <- gamma[cidx]
lines(c(maxbeta,maxbeta), c(0, 1e5), col='red')
lines(c(0,1e5), c(maxgamma,maxgamma), col='red')

# Legend
legend('bottomleft', inset=.05, lty=1,
       c('Volz2012', 'Volz2009', 'Truth'),
       col=c('blue','red', 'black'))

dev.off()
