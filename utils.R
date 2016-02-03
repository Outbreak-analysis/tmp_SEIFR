
read.param <- function(filename){
	x <- read.csv(filename,header=FALSE)
	res <- list()
	for (i in 1:nrow(x)){
		res[[i]] <- x[i,2]
	}
	names(res) <- x[,1]
	return(res)
}

calc.beta <- function(prm){
	# WARNING: quick and dirty to calculate betas (and not correct!)
	beta_IS <- prm[["R0"]]/(1/prm[["DOI.days"]]+1/prm[["funeral.days"]])/2 
	beta_FS <- prm[["rho"]]*beta_IS
	return(c(list(beta_IS=beta_IS,beta_FS=beta_FS),prm)	)
}

simulate.syntdata <- function(prm.model.file,
							  prm.simul.file,
							  prm.report.file,
							  n.data.set = NULL,
							  doplot = TRUE) {
	
	# Read parameters from files:
	simul.prm <- read.param(prm.simul.file)
	model.prm <- calc.beta(read.param(prm.model.file))
	model.prm.true <- model.prm
	
	# overwrites MC iterations read in file
	if (!is.null(n.data.set)) simul.prm[["n.MC"]] <- n.data.set
	
	# Simulate epidemics:
	t1 <- as.numeric(Sys.time())
	sim <- SEIFR.sim(model.prm, simul.prm, seed = 1234)
	t2 <- as.numeric(Sys.time())
	message(paste("Synthetic data sets simulated in",round((t2-t1)/60,2),"minutes"))
	
	# Apply reporting layer:
	report.prm <- read.param(prm.report.file)
	sim2 <- reporting.filter(sim,prm=report.prm)
	
	### Plots just one Monte Carlo iteration
	if (doplot){
		g <- ggplot(sim2) + geom_step(aes(x=tb,y=inc),size=2)
		g <- g + geom_step(aes(x=tb,y=buried),size=2,colour="red") + facet_wrap(~mc)
		g <- g + ggtitle("All synthetic data sets\n incidence and burials")
		plot(g)
		
		g <- ggplot(sim2) + geom_line(aes(x=tb,y=inc+1,colour = factor(mc)),
									  size=1, alpha=0.8)
		g <- g + scale_y_log10() + annotation_logticks(sides="lr")
		g <- g + scale_color_brewer(palette = "Paired")
		g <- g + ggtitle("All synthetic data sets\n log incidence only")
		plot(g)
	}
	return(sim2)
}

