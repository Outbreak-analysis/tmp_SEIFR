library(snowfall)

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
	if (max(obs.data$tb)>=simul.prm[["horizon"]]) 
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
								prm.fixed,
								simul.fcast.prm,
								multi.core=0) {
	###
	### FORECAST ASSUMING FULL REPORTING
	### AND SAMPLING FROM POSTERiORS
	###
	n.prm <- nrow(post.abc$param)
	
	message(paste("Number of posterior samples:",n.prm))
	message(paste("MC iterations for each sample:",simul.fcast.prm[["n.MC"]]))
	message(paste("==>",n.prm*simul.fcast.prm[["n.MC"]],"simulations in total"))
	
	library(parallel)
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
	
	df2 <- ddply(df,c("tb"),summarize, 
				 inc.m = mean(inc),
				 inc.md = median(inc),
				 inc.lo = quantile(inc,probs=0.10),
				 inc.hi = quantile(inc,probs=0.90),
				 buried.m = mean(buried),
				 buried.md = median(buried),
				 buried.lo = quantile(buried,probs=0.10),
				 buried.hi = quantile(buried,probs=0.90)
				 )
	return(df2)
}



plot.tmp <- function(var, x,obs.data, obs.data.full=NULL) {
	
	var.m <- paste0(var,".m")
	var.lo <- paste0(var,".lo")
	var.hi <- paste0(var,".hi")
	
	plot(x$tb, 
		 x[,var.m],
		 ylim = range(obs.data,x[,var.lo],x[,var.hi]),
		 xlab="time",ylab="",main=var)
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

plot.forecast <- function(x, obs.data, obs.data.full=NULL) {
	
	plot.tmp(var="inc", x,obs.data, obs.data.full)
	plot.tmp(var="buried", x,obs.data, obs.data.full)
	
}


