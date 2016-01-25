###
###   SEIFR fit using ABC method from "EasyABC" package
###

library(EasyABC)

source("seifr_gillespie.R")



abc.stats <- function(sim,
					  prm, # input parameters for calculating the summary stats
					  stat.type, # type of summary stats desired
					  do.plot=FALSE){
	###
	### Summary statistics from epidemic data
	###
	
	# Unpack input parameters:
	t.first <- prm[["first.time"]]
	t.last <- prm[["last.time"]]
	t.rng <- t.first:t.last
	
	# Retrieve data:
	t <- sim$tb[t.rng]
	inc <- sim$inc[t.rng]
	bur <- sim$buried[t.rng]
	
	# Poisson regressions on incidence and burials:
	rg.inc <- try(glm(inc ~ t,family = "poisson")[["coefficients"]]	,silent = TRUE)
	rg.bur <- try(glm(bur ~ t,family = "poisson")[["coefficients"]]	,silent = TRUE)
	
	if (class(rg.inc)=="try-error") rg.inc <- 9E9
	if (class(rg.bur)=="try-error") rg.bur <- 9E9
	
	# Maximum and its timing
	inc.mx <- max(inc)
	inc.tmx <- t[which.max(inc)]
	bur.mx <- max(bur)
	bur.tmx <- t[which.max(bur)]
	
	if(do.plot){
		plot(t,inc)
		lines(t,exp(rg.inc[1]+rg.inc[2]*t))
		plot(t,bur)
		lines(t,exp(rg.bur[1]+rg.bur[2]*t))
	}
	
	# Return the stats types requested,
	res <- vector()
	
	if(stat.type[["inc.poisson.reg"]]) res <- c(res,rg.inc=rg.inc)
	if(stat.type[["bur.poisson.reg"]]) res <- c(res,rg.bur=rg.bur)
	if(stat.type[["inc.max"]]) res <- c(res,c(inc.mx=inc.mx,inc.tmx=inc.tmx))
	if(stat.type[["bur.max"]]) res <- c(res,c(bur.mx=bur.mx,bur.tmx=bur.tmx))
	
	return(res)
}


fit.abc <- function(prm.fit, 
					prm.fixed, 
					obs.data,
					prm.stats,
					stat.type,
					priors,
					horizon,  
					n.ABC,
					tol.ABC
){
	# Summary stats of observed data:
	sum.stat.obs <- abc.stats(sim = obs.data,
							  prm = prm.stats,
							  stat.type = stat.type)
	
	prm.fit.vec <- rapply(prm.fit, c)
	str(prm.fit.vec)
	
	wrap.abc <- function(x=prm.fit.vec){
		# Merge all model parameters:
		all.prm <-  as.list(c(prm.fit.vec,prm.fixed))
		# Simulate:
		sim <- SEIFR.sim(model.prm = all.prm,
						 horizon = horizon,
						 n.MC = 1,
						 do.adaptivetau = TRUE,
						 epsilon = 0.05)
		# Calculate summary stats:
		sum.stat <- abc.stats(sim = sim, 
							  prm = prm.stats,
							  stat.type = stat.type)
		return(sum.stat)
	}
	
	# Use of 'EasyABC' package:
	posterior <- ABC_rejection(model = wrap.abc, 
							   prior = priors, 
							   nb_simul = n.ABC,
							   summary_stat_target = sum.stat.obs,
							   tol = tol.ABC,
							   progress_bar = TRUE)
	return(posterior)
}
