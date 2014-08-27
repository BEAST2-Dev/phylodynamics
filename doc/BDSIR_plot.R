library(s20x)
library(boa)
# input : a matrix M and a column ascii name 
# output : the numeric index of the column 
colnameindex = function(M , colname0) 
{ 
colsnames = names(M[1,]); 
theindex = which(colsnames==colname0); 
return(theindex); 
}

loglist = read.table("logs.list", as.is=TRUE, header=F) 

finalSampleTime = 1999.42 # time of most recent sample
dim =length(loglist[,1])

intervals=101


for(i in 1:dim){

	log = read.table(loglist[i,], header=TRUE)
	samples=length(log[,1])
	burnin = round(.1 * samples)
	ds0 = rep(NA,samples-burnin+1)
	dr0 = rep(NA,samples-burnin+1)

	R_perInterval=matrix(data = NA, nrow = dim, ncol = intervals)
	R_perInterval_HPD=matrix(data = NA, nrow = 2, ncol = intervals)
	allS = matrix(NA, ncol=intervals, nrow=samples-burnin+1)

	attach(log)

    S0 = get(names(log)[which(regexpr("S0\\w+", names(log))>0)])

	allS[,1]=S0[burnin:samples]

	
		finalS = matrix(NA,nrow=(samples-burnin+1),ncol=5)
		finalI = matrix(NA,nrow=(samples-burnin+1),ncol=5)
		finalR = matrix(NA,nrow=(samples-burnin+1),ncol=5)

		
		#if (min(S0) < minS) minS = min(S0)
        N= median(S0[burnin:samples])
		T = get(names(log)[which(regexpr("origin\\w+", names(log))>0)])

		if (plot==1){
			plot(1, xlim=c(1983, finalSampleTime), ylim = c(0, N), col="white", xlab='time', ylab='# individuals in compartment', main="SIR trajectories")
			legend("topleft", inset=.01, c("S", "I", "R"), lty=c(1,1,1), col=c('#76bf8a', '#ffa217', '#4f7eff'), bty="n")
		}
	
		
		
# calculations for plotting trajectories:	
		
	
	  S = I = inc = R = t = rep(NA,intervals)
	  S[1] = median(S0[burnin:samples])
	  I[1] = 1
	  inc[1] = 1
	  R[1] = 0
	  t[1] = finalSampleTime - median(T[burnin:samples])
	  step = median(T[burnin:samples])/intervals
	  
	  for (j in 2:intervals){
			dS = get(names(log)[which(grepl("dS",names(log)))[j-1]])
			dR = get(names(log)[which(grepl("dR",names(log)))[j-1]])
			S[j] =S[j-1] - median(dS[burnin:samples])
			R[j] = R[j-1] + median(dR[burnin:samples])
			I[j] = I[j-1] + median(dS[burnin:samples]) - median(dR[burnin:samples])
			inc[j] = median(dS[burnin:samples])
			t[j] = t[j-1] + step
	  }
	  
	layout20x(3,1)  
	  
	# plot SIR trajs
	plot(1, main='', xlim=c(t[1], finalSampleTime), ylim = c(0, N), col="white", ylab='', xlab=paste("(a)",sep=''))
	legend("topleft", inset=.01, c("S", "I", "R"), lty=c(1,1,1), col=c('#76bf8a', '#ffa217', '#4f7eff'), bty="n")
	lines(t,S, col="#76bf8a")
	  lines(t,R, col="#4f7eff")
	  lines(t,I, col="#ffa217")
			
	# plot prevalence (=I) and incidence 
	plot(1, xlim=c(t[1], finalSampleTime), ylim = c(0, max(inc,I)), col="white", ylab='', xlab=paste("(b)",sep=''))
	legend("topleft", inset=.05, c("Prevalence","Incidence"), lty=c(1,1,1), col=c('#ffa217','red'), bty="n")
	lines(t,I, col="#ffa217")
	lines(t,inc, col="red")


# calculations for plotting R over time:	

	R0 = get(names(log)[which(grepl("R0",names(log)))])
	birth = get(names(log)[which(grepl("birth",names(log)))])
	becomeUninfectiousRate = get(names(log)[which(grepl("becomeUninfectiousRate",names(log)))])
	
	for (j in 1:intervals){
		  if (j>1) allS[,j]=allS[,j-1]- get(names(log)[which(grepl("dS",names(log)))[j-1]])[(burnin):samples]
	  
		  R_perInterval[i,j] = median(birth[(burnin):samples] * allS[,j] / becomeUninfectiousRate[(burnin):samples])
		  R_perInterval_HPD[,j] = boa.hpd((R0[(burnin):samples] * (becomeUninfectiousRate[(burnin):samples]) / allS[,1] * allS[,j] / (becomeUninfectiousRate[(burnin):samples]))[which(!is.na(allS[,j]))],0.05)[1:2]
	}		

	#plot median R
	plot(t,R_perInterval[i,], ylim=c(0, max(R_perInterval_HPD,na.rm=T)), type='l', ylab='', xlab='')		
	lines(t,R_perInterval_HPD[1,], lty=3)
	lines(t,R_perInterval_HPD[2,], lty=3)
	abline(1,0,col='grey')

			
	dev.copy2pdf(file=paste("SIR_Plot_Re_",i,".pdf",sep=''),width=9,height=10)		
	  


}

