
#Functions needed to apply Catch-MSY of S Martell and R Froese
# FISH and FISHERIES 2012 (14)4:504-514.

.schaefer	<- function(theta)
{
	with(as.list(theta), {  ## for all combinations of ri & ki
		bt=vector()
		ell = 0  ## initialize ell
		for (j in startbt)
		{
			if(ell == 0) 
			{
				bt[1]=j*k*exp(rnorm(1,0, sigR))  ## set biomass in first year
				for(i in 1:nyr) ## for all years in the time series
				{
					xt=rnorm(1,0, sigR)
					bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
				}
		
				#Bernoulli likelihood, assign 0 or 1 to each combination of r and k
				ell = 0
				if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
				ell = 1
			}	
		}
		return(list(ell=ell))
		
		
	})
}


sraMSY	<-function(theta, N)
{
	#This function conducts the stock reduction
	#analysis for N trials
	#args:
	#	theta - a list object containing:
	#		r (lower and upper bounds for r)
	#		k (lower and upper bounds for k)
	#		lambda (limits for current depletion)
	
	
	with(as.list(theta), 
	{
		ri = exp(runif(N, log(r[1]), log(r[2])))  ## get N values between r[1] and r[2], assign to ri
		ki = exp(runif(N, log(k[1]), log(k[2])))  ## get N values between k[1] and k[2], assing to ki
		itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2], sigR=sigR) ## assign ri, ki, and final biomass range to itheta
		M = apply(itheta,1,.schaefer) ## call Schaefer function with parameters in itheta
		i=1:N
		## prototype objective function
		get.ell=function(i) M[[i]]$ell
		ell = sapply(i, get.ell) 
		return(list(r=ri,k=ki, ell=ell))	
	})
}
