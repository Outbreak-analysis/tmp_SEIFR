#### 
####  Example file to show how to 
####  call the SEIFR model
####
####  Time unit can be anything, but must be consistent!
####

source("seifr_gillespie.R")
source("fit_abc.R")
source("utils.R")
source("forecast.R")

# Read parameters from files:
simul.prm <- read.param("prm_simul.csv")
model.prm <- calc.beta(read.param("prm_model.csv"))

# Simulate epidemics:
sim <- SEIFR.sim(model.prm, simul.prm, seed = 1234)

# Apply reporting layer:
report.prm <- read.param("prm_report.csv")
sim2 <- reporting.filter(sim,prm=report.prm)

### (optional) Plots just one Monte Carlo iteration
if (FALSE){
	mc.chosen <- 1  # Choose any one
	g <- ggplot(subset(sim2,mc==mc.chosen)) + geom_step(aes(x=tb,y=inc),size=2)
	g <- g + geom_step(aes(x=tb,y=buried),size=2,colour="red")
	plot(g)
}

###
###    FAKE CALIBRATION WITH 'ABC'
###

# We pretend that one simulation data set 'mc'
# is an actual observation until a specified horizon ('hz'):
hz <- 19
obs.data <- subset(sim2,mc==2 & tb<=hz)
# How the data look like:
plot(obs.data$tb, obs.data$inc, typ="s",lwd=3)
lines(obs.data$tb, obs.data$bur, typ="s",col="red",lwd=3)

# Specify which parameters will be calibrated
# and their prior distributions:
prm.fit <- list(beta_IS=0.999, 
				beta_FS=0.999)  # <- param values do not matter here
priors <-  list(c("unif",0.05,3),  
				c("unif",0.05,3))

# Fixed parameters (will NOT be calibrated):
# (take the original model params and remove the ones that will be)
prm.fixed  <- model.prm
prm.fixed[which(names(prm.fixed) %in% names(prm.fit))] <- NULL  

horizon <- hz+20  
n.MC <- 7
n.ABC <- 50
tol.ABC <- 0.20

# Summary statistics definition.
# What kind of summary statistics
# will be used to assess the proximity
# to observed (target) data
#
# Time range where the 
# summary stats are aplied:
prm.stats <- list(first.time = 10, 
				  last.time = hz)
# Type of summary stats
# inc.poisson.reg : poisson regression on incidence
# bur.poisson.reg : poisson regression on burials
# inc.max : level and timing of max incidence
# bur.max : level and timing of max burials

stat.type <- list(inc.poisson.reg = TRUE,
				  bur.poisson.reg = FALSE,
				  inc.max = FALSE,
				  bur.max = FALSE)

# Swicth DC's hack on and off:
# (the hack is about fixing an issue
# in EasyABC package when run on multi-cores)
do.DC_HACK <- TRUE
if(do.DC_HACK){
	library(parallel)
	source("problem-easyabc/EasyABC-internal-DC_HACK.R")
	source("problem-easyabc/ABC_rejection-DC_HACK.R")
}

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
					multi.core = 0 
)

t2 <- as.numeric(Sys.time())
message(paste0("ABC fit done in ",round((t2-t1)/60,2)," minutes."))

save.image("sim.RData")

# Visualize posteriors:
par(mfrow=c(1,3))
beta_IS <- model.prm[["beta_IS"]] 
beta_FS <- model.prm[["beta_FS"]] 
plot(post.abc$param[,1],post.abc$param[,2])
abline(v=mean(post.abc$param[,1]))
abline(h=mean(post.abc$param[,2]))
points(beta_IS, beta_FS,col="red",cex=5,pch=16)
hist(post.abc$param[,1],breaks=12,col="grey")
abline(v=beta_IS,col="red")
hist(post.abc$param[,2],breaks=12,col="grey")
abline(v=beta_FS,col="red")


# Forecast

simul.fcast.prm = simul.prm
simul.fcast.prm[["horizon"]] <- max(obs.data$tb)+10
simul.fcast.prm[["n.MC"]] <- 15

x <- forecast.fullreport(obs.data,
						 post.abc,
						 prm.fit,
						 prm.fixed,
						 simul.fcast.prm)
par(mfrow=c(1,1))
plot.forecast(x, obs.data)

