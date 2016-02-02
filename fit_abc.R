###
###   SEIFR fit using ABC method from "EasyABC" package
###

library(EasyABC)
library(parallel)

source("seifr_gillespie.R")

abc.stats <- function(sim,
					  prm, # input parameters for calculating the summary stats
					  stat.type, # type of summary stats desired
					  do.plot = FALSE) {
	###
	### Summary statistics from epidemic data
	###
	
	res <- vector()
	
	# Unpack input parameters:
	t.first <- prm[["first.time"]]
	t.last <- prm[["last.time"]]
	t.rng <- t.first:t.last
	
	# Retrieve data:
	t <- sim$tb[t.rng]
	inc <- sim$inc[t.rng]
	bur <- sim$buried[t.rng]
	
	# Poisson regressions on incidence and burials:
	if (stat.type[["inc.poisson.reg"]]){
		rg.inc <-
			try(glm(inc ~ t,family = "poisson")[["coefficients"]]	,silent = TRUE)
		if (class(rg.inc) == "try-error") rg.inc <- 9E9
		res <- c(res,rg.inc = rg.inc)
	}
	if (stat.type[["bur.poisson.reg"]]){
		rg.bur <-
			try(glm(bur ~ t,family = "poisson")[["coefficients"]]	,silent = TRUE)	
		if (class(rg.bur) == "try-error") rg.bur <- 9E9
		res <- c(res,rg.bur = rg.bur)
	}
	# Maximum and its timing
	if (stat.type[["inc.max"]]){
		inc.mx <- max(inc)
		inc.tmx <- t[which.max(inc)]
		res <- c(res,c(inc.mx = inc.mx,inc.tmx = inc.tmx))
	}
	if (stat.type[["bur.max"]]){
		bur.mx <- max(bur)
		bur.tmx <- t[which.max(bur)]
		res <- c(res,c(bur.mx = bur.mx,bur.tmx = bur.tmx))
	}
	
	if (do.plot) {
		if (stat.type[["inc.poisson.reg"]]) {
			plot(t,inc)
			lines(t,exp(rg.inc[1] + rg.inc[2] * t))
		}
		if (stat.type[["bur.poisson.reg"]]) {
			plot(t,bur)
			lines(t,exp(rg.bur[1] + rg.bur[2] * t))
		}
	}
	return(res)
}


get.clustername <- function() {
	res <- NA
	objlist <- ls(envir = .GlobalEnv)
	# print(objlist)
	for (i in 1:length(objlist)) {
		tmp <- get(objlist[i])
		print(class(tmp))
		if (any(grepl("cluster",class(tmp)))) {
			res <- objlist[i]
			break
		}
	}
	return(res)
}



fit.abc <- function(prm.fit,
					prm.fixed,
					obs.data,
					prm.stats,
					stat.type,
					priors,
					horizon,
					n.MC,  # number of Monte Carlo iterations per sampled parameter set
					n.ABC, # number of sampled parameter sets
					tol.ABC,
					tau.espilon,
					multi.core = 0,
					do.plot = FALSE) {
	# Summary stats of observed data:
	sum.stat.obs <- NULL
	if (!is.null(obs.data)) {
		sum.stat.obs <- abc.stats(sim = obs.data,
								  prm = prm.stats,
								  stat.type = stat.type,
								  do.plot = do.plot)
	}
	prm.fit.vec <- rapply(prm.fit, c)
	
	wrap.abc <- function(x) {
		# Initialize random generator:
		set.seed(x[1])
		# Get rid of the seed
		# ("EasyABC" forces the seed to be 
		# the first parameter of the function)
		x <- x[2:length(x)]
		
		# Merge all model parameters:
		all.prm <-  as.list(c(x, prm.fixed))
		names(all.prm)[1:length(x)] <- names(prm.fit)

		simul.prm <- list(
			horizon = horizon,
			n.MC = n.MC,
			do.adaptivetau = TRUE,
			epsilon = tau.espilon
		)
		# Simulate:
		sim <- SEIFR.sim(model.prm = all.prm,
						 simul.prm = simul.prm)
		# Take average of all MC iterations:
		sim.avg <- ddply(sim, .variables = "tb", summarize,
						 inc = mean(inc),
						 buried = mean(buried))
		# Calculate summary stats:
		sum.stat <- abc.stats(sim = sim.avg,
							  prm = prm.stats,
							  stat.type = stat.type, 
							  do.plot = do.plot)
		# DEBUG:
# 		print("currvalue:"); print(x)
# 		print("currStat:"); print(sum.stat)
		return(sum.stat)
	}
	
	# Use of 'EasyABC' package:
	maxcores <- detectCores()
	if (multi.core > 0) nc <- min(multi.core,maxcores)
	if (multi.core <= 0) nc <- max(1,maxcores - multi.core)
	
	posterior <- ABC_rejection(
		model = wrap.abc,
		use_seed = TRUE, #(nc > 1),		
		# seed_count = 
		n_cluster = nc,
		prior = priors,
		nb_simul = n.ABC,
		summary_stat_target = sum.stat.obs,
		tol = tol.ABC,
		progress_bar = TRUE
	)
	return(posterior)
}
