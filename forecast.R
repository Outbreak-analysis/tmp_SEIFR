source("seifr_gillespie.R")

forecast.after.fit <- function(obs.data,
							   post.abc,
							   prm.fit,
							   prm.fixed,
							   simul.prm){
	
	# Check forecasting horizon is long enough
	
	for (i in 1:nrow(post.abc)){
		
		prm.fit.i <- as.list(post.abc$param[i,])
		names(prm.fit.i) <- names(prm.fit)
		model.prm.i <- c(prm.fit.i, prm.fixed)
		sim.i <- SEIFR.sim(model.prm = model.prm.i, 
						   simul.prm = simul.prm)
	}
	
	
}