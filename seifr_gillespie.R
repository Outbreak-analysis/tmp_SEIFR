###
###
###   SEIFR MODEL WITH GILLESPIE IMPLEMENTATION
###   (USING adaptivetau R PACKAGE)
###   
###   Generate several Monte-Carlo iterations
###   and store the results in a dataframe.
###
###


library(ggplot2);theme_set(theme_bw())
library(plyr)
library(adaptivetau)  


lappend <- function (lst, ...){
	lst <- c(lst, list(...))
	return(lst)
}

trans.SEIFR <- function(nE,nI,nF){
	###
	### Generate the list of transitions for SEIFR 
	###
	
	# infection:
	z <- list(c(S=-1, E1=1))
	
	# transition through E_k compartments
	for(i in 2:nE){
		tmp <- c(-1,1)
		names(tmp) <- c(paste0("E",i-1),paste0("E",i))
		z <- lappend(z,tmp)
	}
	
	# infectiousness triggered:
	tmp <- c(-1,1)
	names(tmp) <-  c(paste0("E",nE),"I1")
	z <- lappend(z,tmp)
	
	# transition through I_k compartments
	for(i in 2:nI){
		tmp <- c(-1,1)
		names(tmp) <- c(paste0("I",i-1),paste0("I",i))
		z <- lappend(z,tmp)
	}
	
	# Recovery
	tmp <- c(-1,1)
	names(tmp) <-  c(paste0("I",nI),"R")
	z <- lappend(z,tmp)
	
	# Death
	tmp <- c(-1,1)
	names(tmp) <-  c(paste0("I",nI),"F1")
	z <- lappend(z,tmp)
	
	# transition through F_k compartments
	for(i in 2:nF){
		tmp <- c(-1,1)
		names(tmp) <- c(paste0("F",i-1),paste0("F",i))
		z <- lappend(z,tmp)
	}
	
	# Burial
	tmp <- c(-1,1)
	names(tmp) <-  c(paste0("F",nF),"B")
	z <- lappend(z,tmp)
	
	return(z)
}

trans.rate.SEIFR <- function(x, prm, t){
	###
	### Generate transition rates for SEIFR 
	###
	###   * WARNING * : same order as transitions!
	
	# unpack all parameters:
	b_IS <- prm$beta_IS
	b_FS <- prm$beta_FS
	sigma <- prm$sigma
	gam <- prm$gamma
	nE <- prm$nE
	nI <- prm$nI
	nF <- prm$nF
	delta <- prm$delta
	f <- prm$f
	
	# change in contact rate (behaviour, iterventions,...)
	ccr <- prm$chgcontact
	t.ccr <- prm$t.chgcontact
	
	# total population
	N <- sum(x)
	
	# prevalence infectious alive
	idx.I <- which(grepl(pattern = "I",x = names(x)))
	Itot <- sum(x[idx.I])
	# prevalence infectious dead in funerals
	idx.F <- which(grepl(pattern = "F",x = names(x)))
	Ftot <- sum(x[idx.F])
	
	# infection rate
	inf.rate <- (b_IS*Itot+b_FS*Ftot)*x["S"]/N*(x["S"]>=1)
	
	if(!is.null(t.ccr)){
		if(t>=t.ccr) inf.rate <- inf.rate*exp(ccr*t)
	}
	
	# transition through E_k
	trans.E <- numeric(nE)
	for(i in 1:nE){
		trans.E[i] <- sigma*nE * x[paste0("E",i)]*(x[paste0("E",i)]>=1)
	}
	
	# transition through I_k
	# (because transition to either F or R, 
	#  last compartment treated separately [hence 'nI-1'])
	trans.I <- numeric(nI-1)
	for(i in 1:(nI-1)){
		trans.I[i] <- sigma*nI * x[paste0("I",i)]*(x[paste0("I",i)]>=1)
	}
	
	
	# Recovery (I_nI)
	recov.rate <- (1-delta)*gam*nI*x[paste0("I",nI)]*(x[paste0("I",nI)]>=1)
	
	# Death (I_nI)
	death.rate <- delta*gam*nI*x[paste0("I",nI)]*(x[paste0("I",nI)]>=1)
	
	# transition through F_k compartments
	# (last element is for Burial)
	trans.F <- numeric(nF)
	for(i in 1:nF){
		trans.F[i] <- f*nF * x[paste0("F",i)]*(x[paste0("F",i)]>=1)
	}
	
	
	# Gather all transition rates that match
	# transition events defined elsewhere:
	return(c(inf.rate,
			 trans.E,
			 trans.I, 
			 recov.rate,
			 death.rate,
			 trans.F))
}

SEIFR.sim <-function(beta_IS, # contact rate I -> S
					 beta_FS, # contact rate F -> S
					 DOL.days,  # avg duration of latency
					 DOI.days,  # avg duration of infectiousness
					 funeral.days, # avg duration of funerals
					 horizon,  # horizon of the simulation
					 n.MC,   # Monte carlo iterations
					 delta , # proportion infected that will die
					 pop.size, # effective population size
					 I.init,  # initial number of infectious individuals
					 nE,   # number of (artificial) compartments for E
					 nI,   # number of (artificial) compartments for I
					 nF,   # number of (artificial) compartments for F
					 chgcontact = NULL,  # temporal change in the contact rates (intervention,behaviour,...)
					 t.chgcontact = NULL,  # time when the change starts
					 seed = 1234,
					 do.adaptivetau = TRUE, # FALSE = true Gillespie algo ; TRUE=approximation
					 epsilon = 0.05, # larger=faster but less accurate
					 time.bucket = 1,  # aggregation of incidence in time units (Gillespie events happens any time)
					 remove.fizzles = FALSE  # remove Monte carlo iterations that are fizzles
){
	###
	### RUN SEIFR SIMULATIONS
	###
	
	set.seed(seed)
	
	gamma <- 1/DOL.days
	sigma <- 1/DOI.days
	f <- 1/funeral.days
	
	params <- list(beta_IS=beta_IS,
				   beta_FS=beta_FS,
				   sigma=sigma,
				   gamma=gamma, 
				   nE=nE, nI=nI, nF=nF,
				   delta=delta,
				   f=f,
				   chgcontact = chgcontact,
				   t.chgcontact = t.chgcontact)
	
	# Initial values
	I0 <- I.init
	x0 <- c(pop.size-I0,
			rep(0,nE), 
			I0,
			rep(0,nI-1),
			rep(0,nF), 
			0, # <-- R (removed: alive & immune after infection)
			0  # <-- B (buried: dead after infection)
	)
	names(x0) <- c("S",
				   paste0("E",1:nE),
				   paste0("I",1:nI),
				   paste0("F",1:nF),
				   "R","B")
	
	# Monte Carlo iterations:
	
	for(mc in 1:n.MC){
		print(paste("MC",mc,"/",n.MC))
		
		if (!do.adaptivetau){
			res <- ssa.exact(init.values = x0,
							 transitions = trans.SEIFR(nE,nI,nF), 
							 rateFunc = trans.rate.SEIFR, 
							 params = params,
							 tf=horizon)
		}
		if(do.adaptivetau){
			res <- ssa.adaptivetau(init.values = x0,
								   transitions = trans.SEIFR(nE,nI,nF), 
								   rateFunc = trans.rate.SEIFR, 
								   params = params,
								   tl.params = list(epsilon=epsilon),  
								   tf=horizon)
		}
		
		idx.E <- grepl("E",colnames(res))
		idx.I <- grepl("I",colnames(res))
		idx.F <- grepl("F",colnames(res))
		idx.EIF <- as.logical(idx.E+idx.I+idx.F)
		idx.S <- which(grepl("S",colnames(res)))
		idx.B <- which(grepl("B",colnames(res)))
		
		tmp <- data.frame(t = res[,1], 
						  inc = c(I0,-diff(res[,idx.S])),
						  prev = apply(res[,idx.EIF],1,sum),
						  prevI = apply(res[,idx.I],1,sum),
						  cumburied = res[,idx.B],
						  mc = rep(mc,nrow(res)))
		
		tmp$cuminc <- cumsum(tmp$inc)
		tmp$buried <- c(0,diff(tmp$cumburied))
		
		### Aggregate in time buckets
		tmp$tb <- ceiling(tmp$t/time.bucket)
		tmp$tb[tmp$tb==0] <- 1
		
		# Store all results in dataframe
		if(mc==1) all.sim <- tmp
		if(mc>1) all.sim <- rbind(all.sim, tmp)
	}
	
	### Remove fizzles:
	if(remove.fizzles){
		# identify the fizzles
		thres <- 0.1
		sim.nofizz = ddply(all.sim,.variables = c("mc"),
						   summarize,
						   i.max = max(prev))
		imax.all <- max(sim.nofizz$i.max)
		mc.nofizz <- which(sim.nofizz$i.max>thres*imax.all)
		# Probability of fizzles (not used)
		p.fizz <- 1-length(mc.nofizz)/n.MC
		# Filter out the fizzles:
		all.sim.nofizz <- all.sim[all.sim$mc %in% mc.nofizz,]
		all.sim <- all.sim.nofizz
	}
	### Incidence only at time buckets ('tb')
	inc.tb <- ddply(all.sim,c("tb","mc"),summarize,
					inc = sum(inc),
					cuminc = max(cuminc),
					buried = sum(buried),
					cumburied = max(cumburied))
	
	return(inc.tb)
}


lag.fct <- function(x,lag.mean,lag.var, seed=1234){
	
	set.seed(seed)
	
	x.lag <- rep(0,times=2*length(x))
	x.lag[1:length(x)] <- x
	
	for(t in 1:length(x)){
		
# 		print(paste("time",t))
# 		print(length(x.lag))
		
		
		if(x[t]>0){
			# If incidence positive, calculate lags for each case at that date:
			lag <- round(rlnorm(n = x[t],
								meanlog = log(lag.mean),
								sdlog = lag.var), 
						 digits = 0)
			# print(lag) # DEBug
			for(k in 1:length(lag)){
				# Apply the lag for each case
				if(lag[k]>0){
					# Move incidence, one by one, to lagged reporting date
					x.lag[t] <- x.lag[t]-1
					x.lag[t+lag[k]] <- x.lag[t+lag[k]] + 1
				}
			}
		}
	}
	return(x.lag)
}


reporting.filter <- function(sim,
							 report.inc.prob, 
							 report.inc.lag.mean, 
							 report.inc.lag.var, 
							 report.bur.prob, 
							 report.bur.lag.mean,
							 report.bur.lag.var,
							 seed = 1234,
							 do.plot = FALSE){
	###
	### Add a reporting "layer" to the simulated epidemic
	###
	
	set.seed(seed)
	n.mc <- max(sim$mc)
	
	for(m in 1:n.mc) {
		tmp <- subset(sim,mc==m)
		inc <- tmp$inc
		bur <- tmp$buried
		
		# Reduce actual numbers to reported ones
		inc.rep <- rbinom(n = length(inc),size = inc, prob = report.inc.prob)
		bur.rep <- rbinom(n = length(bur),size = bur, prob = report.bur.prob)
		
		# Introduce lag
		inc.rep.lag <- lag.fct(x = inc.rep, 
							   lag.mean = report.inc.lag.mean, 
							   lag.var = report.inc.lag.var)
		bur.rep.lag <- lag.fct(x = bur.rep, 
							   lag.mean = report.bur.lag.mean, 
							   lag.var = report.bur.lag.var)
		
		# Format data frame same as input 'sim'
		sim.tmp <- data.frame(tb=1:length(inc.rep.lag), 
							  inc=inc.rep.lag, 
							  buried=bur.rep.lag)
		sim.tmp$mc <- m
		sim.tmp$cuminc <- cumsum(sim.tmp$inc)
		sim.tmp$cumburied <- cumsum(sim.tmp$buried)
		# stack all data frames together
		if(m==1) sim2 <- sim.tmp
		if(m>1) sim2 <- rbind(sim2,sim.tmp)
		
		if (do.plot){
			plot(inc,typ="o")
			lines(inc.rep,col="orange")
			lines(inc.rep.lag[1:length(inc)],col="red",lwd=3)
			
			lines(bur,col="navyblue",typ="o")
			lines(bur.rep,col="lightblue")
			lines(bur.rep.lag,col="blue",lwd=3)
			abline(h=0, lty=2)		
		}
	}
	return(sim2)
}



