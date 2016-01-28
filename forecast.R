source("seifr_gillespie.R")


simul.fcast.prm = simul.prm
simul.fcast.prm[["horizon"]] <- max(obs.data$tb)+10

forecast.after.fit <- function(obs.data,
							   post.abc,
							   prm.fit,
							   prm.fixed,
							   simul.fcast.prm){

	# Check forecasting horizon is long enough
	if (max(obs.data$tb)>=simul.prm[["horizon"]]) 
		stop("Forecasting horizon is not beyond observed data!")
	
	# number of posterior samples
	nps <- nrow(post.abc$param)
	
	sim.post <- list()
	for (i in 1:nps){
		print(i)
		# Reconstruct model parameters:
		prm.fit.i <- as.list(post.abc$param[i,])
		names(prm.fit.i) <- names(prm.fit)
		model.prm.i <- c(prm.fit.i, prm.fixed)
		
		# Simulate
		sim.post[[i]] <- SEIFR.sim(model.prm = model.prm.i, 
								   simul.prm = simul.fcast.prm)
		sim.post[[i]]$key <- paste(i,sim.post[[i]]$mc,sep="-")
		
	}
	
	df <- do.call("rbind", sim.post)
	
	t.last <- max(obs.data$tb)
	last.rng <- 3
	t.keep <- c((t.last-last.rng+1):t.last)
	
	# Calculate 'true' incidence that
	# is consistent with reporting rate 
	# and observed data
	report.rate.inc <- 0.5
	report.rate.bur <- 0.9
	CI <- 0.98
	q <- 1-CI
	obs.data$inc.lo <- obs.data$inc + qnbinom(p=q/2, size = obs.data$inc, prob = report.rate.inc)
	obs.data$inc.hi <- obs.data$inc + qnbinom(p=1-q/2, size = obs.data$inc, prob = report.rate.inc)
	inc.lo <- obs.data$inc.lo[t.keep]
	inc.hi <- obs.data$inc.hi[t.keep]
	
	# Keep only simulations that
	# go through the latest observed data
	# TO DO: ADD BURIAL NUMBERS
	keys <- unique(df$key)
	nk <- length(keys)
	key.keep <- list()
	cnt <- 1
	for(i in 1:nk){
		# Retrieve just the time window requested for accepting simulations:
		tmp <- subset(df,key==keys[i] & tb<=t.last & tb>=t.keep[1] )
		# Check if simulated incidence is consistent with
		# estimated (true) target incidence  (estimated from observed&reporting rate)
		check.inc <- FALSE
		if (nrow(tmp)>length(inc.lo)){
			check.inc <- all(tmp$inc>= inc.lo & tmp$inc<=inc.hi)
			if(length(tmp$inc)!=length(inc.lo)) print(i)
		}
		# If satisfy all conditions, then accept this simulation:
		if (check.inc) {
			key.keep[[cnt]] <- keys[i]
			cnt <- cnt+1
		}
	}

	
	df.keep <- df[df$key %in% key.keep,]
	
	g <- ggplot(df.keep) + geom_line(aes(x=tb, y=inc,dummy=factor(key)),alpha=0.3)
	g <- g + geom_line(aes(x=tb, y=inc*report.rate.inc,dummy=factor(key)),alpha=0.3,colour="red")
	g <- g + geom_step(data=obs.data,aes(x=tb,y=inc),colour="red",size=2)
# 	g <- g + geom_line(data=obs.data,aes(x=tb,y=inc.hi),colour="orange",size=1)
# 	g <- g + geom_line(data=obs.data,aes(x=tb,y=inc.lo),colour="orange",size=2)
	g <- g + geom_point(data=obs.data,aes(x=tb,y=buried),colour="blue",size=2)
	g <- g + geom_segment(data=obs.data,aes(x=tb,xend=tb,yend=inc.hi,y=inc.lo),alpha=0.8,colour="orange",size=2)
	# g <- g + scale_y_log10()
	plot(g)
	
	
	
	
	
}