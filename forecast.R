library(snowfall)
library(parallel)
source("seifr_gillespie.R")

forecast.one.postsample <- function(obs.data,
									param,
									prm.fit,
									prm.fixed,
									simul.fcast.prm){
	###
	### FORECAST BASED ON ONE FITTED PARAMETER VECTOR ('param')
	###
	### Returns a data frame of n.MC simulations run with
	### the unique parameter set provided
	
	
	# Check if forecasting horizon is long enough
	if (max(obs.data$tb)>=simul.fcast.prm[["horizon"]]) 
		stop("Forecasting horizon is not beyond observed data!")
	
	print(paste("Each forecast made with",
				simul.fcast.prm[["n.MC"]],
				"MC iterations for a given posterior sample parameter"))
	
	
	# Reconstruct model parameters:
	prm.fit.tmp <- as.list(param)
	names(prm.fit.tmp) <- names(prm.fit)
	model.prm <- c(prm.fit.tmp, prm.fixed)
	
	# Simulate
	sim.post <- SEIFR.sim(model.prm = model.prm, 
						  simul.prm = simul.fcast.prm)
	return(sim.post)
}

forecast.fullreport <- function(obs.data,
								post.abc,
								prm.fit,
								prm.model.file,
								prm.fcast.file,
								multi.core=0) {
	###
	### FORECAST ASSUMING FULL REPORTING
	### AND SAMPLING FROM POSTERiORS
	###
	t1 <- Sys.time()
	# Unpack forecasting parameters:
	pfcst <- read.csv(prm.fcast.file,header=FALSE)
	simul.fcast.prm <- list()
	simul.fcast.prm[["horizon"]] <- pfcst[pfcst[,1]=="horizon",2]  
	simul.fcast.prm[["n.MC"]] <- pfcst[pfcst[,1]=="n.MC",2]  
	simul.fcast.prm[["epsilon"]] <- pfcst[pfcst[,1]=="tau.epsilon",2]  
	simul.fcast.prm[["do.adaptivetau"]] <- TRUE
	simul.fcast.prm[["time.bucket"]] <- 1
	simul.fcast.prm[["remove.fizzles"]] <- 0
	
	n.prm <- nrow(post.abc$param)
	message(paste("Number of posterior samples:",n.prm))
	message(paste("MC iterations for each sample:",simul.fcast.prm[["n.MC"]]))
	message(paste("==>",n.prm*simul.fcast.prm[["n.MC"]],"simulations in total"))
	
	# Determine the fixed parameters
	model.prm <- calc.beta(read.param(prm.model.file))
	prm.fixed  <- model.prm
	prm.fixed[which(names(prm.fixed) %in% names(prm.fit))] <- NULL  
	
	# Parallel setup:
	maxcores <- detectCores()
	if (multi.core > 0) nc <- min(multi.core,maxcores)
	if (multi.core <= 0) nc <- max(1,maxcores - multi.core)
	sfInit(parallel = (nc>1), cpu = nc)
	sfLibrary(adaptivetau)
	sfLibrary(plyr)
	
	snow.wrap <- function(i){
		x <- forecast.one.postsample(obs.data,
									 param = post.abc$param[i,],
									 prm.fit,
									 prm.fixed,
									 simul.fcast.prm)
		x$prmset <- i
		return(x)
	}
	
	# Run all simulations for each posterior parameter set:
	idx.apply <- 1:n.prm
	sfExportAll()
	res <- sfSapply(idx.apply, snow.wrap, simplify = FALSE)
	sfStop()
	# Merge all results in one data frame
	df <- do.call("rbind",res)
	
	# Summarize across all
	# Monte Carlo iterations and
	# posterior samples:
	df.ensemble <- ddply(df,c("tb"),summarize, 
						 inc.m = mean(inc),
						 inc.md = median(inc),
						 inc.lo = quantile(inc, probs=0.1),
						 inc.hi = quantile(inc, probs=0.9),
						 buried.m = mean(buried),
						 buried.md = median(buried),
						 buried.lo = quantile(buried, probs=0.1),
						 buried.hi = quantile(buried, probs=0.9)
	)
	
	df.peak <- ddply(df,c("mc"),summarize, 
					 inc.pk=max(inc),
					 t.inc.pk=tb[which.max(inc)],
					 buried.pk=max(buried),
					 t.buried.pk=tb[which.max(buried)])
	
	t2 <- Sys.time()
	message(paste("\n Forecasting done in",round(as.numeric(t2-t1)/60,2),"minutes"))
	
	return(list(df.ensemble = df.ensemble,
				df.peak = df.peak))
}

ME <- function(target,y,relative) {
	# Mean error
	# target : benchmark
	n <- length(target)
	stopifnot(n == length(y))
	res <- sum(y-target)/n
	if(relative) res <- sum(y/target-1)/n
	return(res)
}

MAE <- function(target,y,relative) {
	# Mean absolute error
	n <- length(target)
	stopifnot(n == length(y))
	res <- sum(abs(target-y))/n
	if(relative) res <- sum(abs(y/target-1))/n
	return(res)
}

MQE <- function(target, y.lo,y.hi,relative) {
	# Mean absolute error
	n <- length(target)
	stopifnot(n == length(y.lo))
	res <- sum(max(0,target-y.hi))/n + sum(max(0,y.lo - target))/n
	if(relative) res <- sum(max(0,target-y.hi)/target)/n + sum(max(0,y.lo - target)/target)/n
	return(res)
}


assess.fcast.ens <- function(var, 
							 x,
							 obs.data, 
							 obs.data.full,
							 relative) {
	### ASSESS THE FORECASTING QUALITY
	### WHEN FULL DATA SET IS KNOWN
	
	# Crop to relevant data
	# (full data set may be larger)
	obs.data.full2 <- obs.data.full[1:nrow(x),]
	
	last.obs <- max(obs.data$tb)
	
	xx <- subset(x, tb > last.obs)
	target <- subset(obs.data.full2, tb > last.obs)
	
	var.m <- paste0(var,".m")
	var.lo <- paste0(var,".lo")
	var.hi <- paste0(var,".hi")
	
	s <- MAE(target = target[,var], y = xx[,var.m], relative) 
	s <- s + MQE(target = target[,var], y.lo = xx[,var.lo],y.hi = xx[,var.hi],relative)
	
	return(c(
		b = ME(target = target[,var], y = xx[,var.m], relative),
		s = s
	))

}

assess.fcast <- function(x,
						 obs.data, 
						 obs.data.full,
						 relative){
	
	inc <- assess.fcast.ens (var="inc", x, obs.data, obs.data.full, relative) 
	bur <- assess.fcast.ens (var="buried", x, obs.data, obs.data.full, relative) 
	
	return(list(inc = inc, bur = bur ))
}



plot.ens <- function(var, 
					 x,
					 obs.data, 
					 obs.data.full=NULL) {
	
	var.m <- paste0(var,".m")
	var.lo <- paste0(var,".lo")
	var.hi <- paste0(var,".hi")
	
	plot(x$tb, 
		 x[,var.m],
		 ylim = range(obs.data,x[,var.lo],x[,var.hi]),
		 xlab="time",ylab="",main=paste("Time series",var))
	lines(x$tb,x[,var.lo])
	lines(x$tb,x[,var.hi])
	grid()
	
	if(!is.null(obs.data.full)){
		points(obs.data.full$tb,obs.data.full[,var],
			   pch=3,col="red")
	}
	points(obs.data$tb,obs.data[,var],
		   pch=3,col="red",lwd=5)
	
	abline(v=obs.data$tb[length(obs.data$tb)],lty=2)
}

plot.peak <- function(var, 
					  X,
					  obs.data.full=NULL,
					  do.hist = FALSE){
	
	x <- X[["df.peak"]] # peak calculations based on RAW simulation
	xe <- X[["df.ensemble"]] # peak calculations based on ENSEMBLE of simulations
	
	var.pk <- paste0(var,".pk")
	var.t.pk <- paste0("t.",var,".pk")
	
	qlo <- 0.1
	qhi <- 0.9
	
	P <- x[,var.pk]
	TP <- x[,var.t.pk]
	P.m <- mean(P)
	P.lo <- quantile(P,probs=qlo)
	P.hi <- quantile(P,probs=qhi)
	TP.m <- mean(TP)
	TP.lo <- quantile(TP,probs=qlo)
	TP.hi <- quantile(TP,probs=qhi)
	
	var.m <- paste0(var,".m")
	var.lo <- paste0(var,".lo")
	var.hi <- paste0(var,".hi")
	
	PE.m <- max(xe[,var.m])
	PE.lo <- max(xe[,var.lo])
	PE.hi <- max(xe[,var.hi])
	
	TPE.m <- which.max(xe[,var.m])
	TPE.lo <- which.max(xe[,var.lo])
	TPE.hi <- which.max(xe[,var.hi])
	
	Prng <- range(P)
	if(!is.null(obs.data.full)) Prng <- range(P,PE.m,PE.lo,PE.hi,max(obs.data.full[,var]))
	TPrng <- range(TP)
	if(!is.null(obs.data.full)) TPrng <- range(TP,TPE.m,TPE.lo,TPE.hi,which.max(obs.data.full[,var]))
	
	# From raw forecast simulations:
	plot(x=TP.m, y=P.m, 
		 lwd=3, pch=15,
		 xlim=TPrng, ylim=Prng, 
		 xlab="Timing peak", ylab="Peak",
		 main=paste("Peak",var))
	polygon(x = c(TP.lo,TP.hi,TP.hi,TP.lo), 
			y = c(P.lo,P.lo,P.hi,P.hi),
			col = rgb(0,0,0,0.1))
	
	# From ensemble 
	# (which does not make much sens, but plot it anyway)
	points(TPE.m,PE.m)
	polygon(x = c(TPE.lo,TPE.hi,TPE.hi,TPE.lo), 
			y = c(PE.lo,PE.lo,PE.hi,PE.hi),
			lty=2, border="darkgrey")
	
	if(!is.null(obs.data.full)) 
		points(which.max(obs.data.full[,var]),
			   max(obs.data.full[,var]),
			   col="red",pch=3,cex=2,lwd=6)
	
	if(do.hist){
		hist(P,
			 xlim=xrng,
			 xlab=var.pk,
			 ylab="",
			 main ="",
			 col="lightgrey",border="grey")
		if(!is.null(obs.data.full)){
			abline(v=max(obs.data.full[,var]),col="red",lwd=2)
		}
		abline(v=P.m, lty=2,lwd=3)
		abline(v=P.lo, lty=3,lwd=2)
		abline(v=P.hi, lty=3,lwd=2)
		
		hist(TP,
			 xlim=t.xrng,
			 xlab=var.t.pk,
			 ylab="",
			 main ="",
			 col="lightgrey",border="grey")
		if(!is.null(obs.data.full)){
			abline(v=which.max(obs.data.full[,var]),col="red",lwd=2)
		}
		abline(v=TP.m, lty=2,lwd=3)
		abline(v=TP.lo, lty=3,lwd=2)
		abline(v=TP.hi, lty=3,lwd=2)
	}
}


plot.forecast <- function(x, obs.data, obs.data.full=NULL) {
	par(mfrow=c(2,2))
	plot.ens(var="inc", x[["df.ensemble"]], obs.data, obs.data.full)
	plot.ens(var="buried", x[["df.ensemble"]], obs.data, obs.data.full)
	
	plot.peak(var="inc", x, obs.data.full)
	plot.peak(var="buried", x, obs.data.full)
}


