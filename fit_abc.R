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
					tau.epsilon,
					multi.core = 0,
					do.plot = FALSE) {
	###
	### ABC fit using 'EasyABC' package
	###
	
	# Integrity checks:
	if (length(prm.fit)!=length(priors)) stop("prm.fit and priors must have the same length")
	if (floor(tol.ABC*n.ABC)<=2) stop("ABC tolerance or n.ABC is too low")
	
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
			epsilon = tau.epsilon
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


fit.abc.seifr <- function(obs.data,
						  prm.fit.file,
						  prm.model.file,
						  prm.fit,
						  priors,
						  stat.type) {
	
	# Fixed parameters (will NOT be calibrated):
	# (take the original model params and remove the ones that will be)
	model.prm <- calc.beta(read.param(prm.model.file))
	prm.fixed  <- model.prm
	prm.fixed[which(names(prm.fixed) %in% names(prm.fit))] <- NULL  
	
	# Read calibration parameters:
	pf <- read.csv(prm.fit.file,header=FALSE)
	horizon <- pf[pf[,1]=="horizon",2]  
	n.MC <- pf[pf[,1]=="n.MC",2]  
	n.ABC <- pf[pf[,1]=="n.ABC",2]  
	tol.ABC <- pf[pf[,1]=="tol.ABC",2]  
	first.time <- pf[pf[,1]=="first.t.stat",2] 
	tau.epsilon <- pf[pf[,1]=="tau.epsilon",2] 
	multi.core <- pf[pf[,1]=="multi.core",2] 
	
	# Time range where the 
	# summary stats are aplied:
	prm.stats <- list(first.time = first.time, 
					  last.time = max(obs.data$tb))
	
	# Calibration with ABC:
	t1 <- as.numeric(Sys.time())
	
	post.abc <- fit.abc(prm.fit, 
						prm.fixed, 
						obs.data,
						prm.stats,
						stat.type,
						priors,
						horizon,  
						n.MC,
						n.ABC,
						tol.ABC,
						tau.epsilon,
						multi.core,
						do.plot = FALSE
	)
	
	t2 <- as.numeric(Sys.time())
	message(paste0("ABC fit done in ",round((t2-t1)/60,2)," minutes."))
	
	save.image("fit_abc_seifr.RData")
	return(post.abc)
}


post.stat <- function(pp,i,j,CI) {
	q.lo <- (1-CI)/2
	q.hi <- 1-q.lo
	x.m <- mean(pp[,i])
	x.lo <- quantile(pp[,i],probs=q.lo)
	x.hi <- quantile(pp[,i],probs=q.hi)
	y.m <- mean(pp[,j])
	y.lo <- quantile(pp[,j],probs=q.lo)
	y.hi <- quantile(pp[,j],probs=q.hi)
	return(list(x.m=x.m,x.lo=x.lo,x.hi=x.hi,y.m=y.m,y.lo=y.lo,y.hi=y.hi))
}

plot.abcfit <- function(post.abc,prm.fit,priors,true.values=NULL,location=NULL){
	###
	### Plot posterior distributions
	
	pp <- as.data.frame(post.abc$param,row.names=c(1:nrow(post.abc$param)))
	colnames(pp) <- names(prm.fit)
	np <- ncol(pp)
	ns <- nrow(pp)
	pr <- as.data.frame(priors)
	
	snp <- ceiling(sqrt(np))
	par(mfrow=c(snp,snp))
	for (i in 1:np) {
		for (j in 1:np) {
			if(j>i){
				# print(paste(i,j))
				title <- "ABC fit results"
				if(!is.null(location)) title <- paste(title,"\n location",location)
				plot(pp[,i],pp[,j],
					 main = title,
					 cex.main = 0.67,
					 xlab = names(prm.fit)[i],
					 ylab = names(prm.fit)[j],
					 xlim = as.numeric(as.character(pr[-1,i])), 
					 ylim = as.numeric(as.character(pr[-1,j])))
				grid()
				z <- post.stat(pp,i,j,CI = 0.90)
				segments(x0 = z[["x.lo"]], x1=z[["x.hi"]], y0=z[["y.m"]], y1=z[["y.m"]])
				segments(x0 = z[["x.m"]], x1=z[["x.m"]], y0=z[["y.lo"]], y1=z[["y.hi"]])
				
				if(!is.null(true.values)){
					points(true.values[[i]],true.values[[j]],pch=5,col="red",cex=2,lwd=4)
				}
				
			}
		}
	}
	
	par(mfrow=c(snp,snp))
	for (i in 1:np) {
		title <- "ABC fit results"
		if(!is.null(location)) title <- paste(title,"\n location",location)
		hist(pp[,i],
			 main=title,
			 cex.main = 0.67,
			 breaks=15,
			 xlim = as.numeric(as.character(pr[-1,i])), 
			 xlab = names(prm.fit)[i], 
			 col="grey",border="darkgrey")
		
		if(!is.null(true.values))
			abline(v=true.values[[i]],col="red",lwd=3)
	}
}





