library(parallel)
source("seifr_gillespie.R")
source("fit_abc.R")
source("utils.R")
source("forecast.R")

pdf.options(width=10,height=10)

# Swicth DC's hack on and off:
# (the hack is about fixing an issue
# in EasyABC package when run on multi-cores)
do.DC_HACK <- TRUE
if(do.DC_HACK){
	source("problem-easyabc/EasyABC-internal-DC_HACK.R")
	source("problem-easyabc/ABC_rejection-DC_HACK.R")
}

# Files where various parameters are defined
prm.model.file <- "prm_model.csv"
prm.simul.file <- "prm_simul.csv"
prm.report.file <- "prm_report.csv"
prm.fit.file <- "prm_fit.csv"
prm.fcast.file <- "prm_forecast.csv"

# Simulate synthetic data:
syntdata <- simulate.syntdata(prm.model.file,
							  prm.simul.file,
							  prm.report.file,
							  n.data.set = 12,
							  doplot = TRUE) 

# Hypothetical last observation time
# of synthetic data:
last.obs <- 19

# Specify which parameters will be calibrated
# and their prior distributions:
prm.fit <- list(beta_IS = 0.999, 
				beta_FS = 0.999,
				DOL.days = 99)  # <- param values do not matter here

priors <-  list(c("unif",0.05,3),  
				c("unif",0.05,3),
				c("unif",1,7))

# Type of summary stats
# inc.poisson.reg : poisson regression on incidence
# bur.poisson.reg : poisson regression on burials
# inc.max : level and timing of max incidence
# bur.max : level and timing of max burials
stat.type <- list(inc.poisson.reg = TRUE,
				  bur.poisson.reg = TRUE,
				  inc.max = TRUE,
				  bur.max = FALSE)

multi.data.fit <- function(syntdata,
						   last.obs,  # (hypothetical) last observation used for synthetic data 
						   prm.fit.file,
						   prm.model.file,
						   prm.fcast.file,
						   prm.fit ,
						   priors,
						   stat.type) {
	###
	###   PERFORM FIT & FORECAST ON SEVERAL SYNTHETIC DATA SETS
	###
	
	tt1 <- as.numeric(Sys.time())
	
	n.data.set <- max(syntdata$mc)
	qf <- list()
	
	for (i in 1:n.data.set) {
		message(paste("   = = = Forecasting data set",i,"/",n.data.set, "   = = = "))
		# Retrieve synthetic data and 
		# crop it to the supposed last observation
		obs.data.full <- subset(syntdata, mc == i)
		obs.data <- subset(obs.data.full, tb < last.obs)
		
		# Fit SEIFR model
		post.abc <- fit.abc.seifr(obs.data,
								  prm.fit.file,
								  prm.model.file,
								  prm.fit ,
								  priors,
								  stat.type)
		# Visualize posteriors:
		model.prm.true <- calc.beta(read.param(prm.model.file))
		true.values <- model.prm.true[names(prm.fit)]
		plot.abcfit(post.abc, prm.fit, priors, true.values)
		
		# Forecast (based on fit done)
		fcast <- forecast.fullreport(obs.data,
									 post.abc,
									 prm.fit,
									 prm.model.file,
									 prm.fcast.file)
		
		plot.forecast(x=fcast, obs.data,obs.data.full)
		
		# Assess forecasting quality
		# (because we know the future of synthetic data)
		qf[[i]] <- assess.fcast(x = fcast[["df.ensemble"]],
								obs.data, 
								obs.data.full,
								relative = FALSE)
	}
	tt2 <- as.numeric(Sys.time())
	message(paste("\n      -- Completed in",round((tt2-tt1)/60,1),"minutes -- \n"))
	return(qf)
}


qf <- multi.data.fit(syntdata,
			   last.obs,  # (hypothetical) last observation used for synthetic data 
			   prm.fit.file,
			   prm.model.file,
			   prm.fcast.file,
			   prm.fit ,
			   priors,
			   stat.type)
s0 <- vector()
b0 <- vector()
for (i in 1:length(qf)) {
	s0[i] <-qf[[i]]$inc['s']
	b0[i] <-qf[[i]]$inc['b']
}
idx <- is.infinite(b0)
b <- b0[!idx] 
s <- s0[!idx]

plot(s,b,pch=16,
	 xlim=range(s,0),
	 ylim=range(b,0))
text(s,b,labels = c(1:length(b)),pos = 3)
abline(h=0)
grid()

save.image("bcktst.RData")


