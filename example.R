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

pdf.options(width=15,height=10)

# Read parameters from files:
simul.prm <- read.param("prm_simul.csv")
model.prm <- calc.beta(read.param("prm_model.csv"))
model.prm.true <- model.prm

# Simulate epidemics:
t1 <- as.numeric(Sys.time())
sim <- SEIFR.sim(model.prm, simul.prm, seed = 1234)
t2 <- as.numeric(Sys.time())
message(paste("first simulations:",round((t2-t1)/60,2),"minutes"))

# Apply reporting layer:
report.prm <- read.param("prm_report.csv")
sim2 <- reporting.filter(sim,prm=report.prm)

### (optional) Plots just one Monte Carlo iteration
if (TRUE){
	mc.chosen <- c(1:9)  # Choose any one
	g <- ggplot(subset(sim2,mc %in% mc.chosen)) + geom_step(aes(x=tb,y=inc),size=2)
	g <- g + geom_step(aes(x=tb,y=buried),size=2,colour="red") + facet_wrap(~mc)
	plot(g)
}

###
###    FAKE CALIBRATION WITH 'ABC'
###

# We pretend that one simulation data set 'mc'
# is an actual observation until a specified horizon ('hz'):
hz <- 19
obs.data.full <- subset(sim2, mc==2)
obs.data <- subset(obs.data.full,tb <= hz)
# How the data look like:
plot(obs.data$tb, obs.data$inc, typ="s",lwd=3)
lines(obs.data$tb, obs.data$bur, typ="s",col="red",lwd=3)

# Specify which parameters will be calibrated
# and their prior distributions:
prm.fit <- list(beta_IS = 0.999, 
				beta_FS = 0.999,
				DOL.days = 99)  # <- param values do not matter here

priors <-  list(c("unif",0.05,3),  
				c("unif",0.05,3),
				c("unif",1,7))

# Fixed parameters (will NOT be calibrated):
# (take the original model params and remove the ones that will be)
prm.fixed  <- model.prm
prm.fixed[which(names(prm.fixed) %in% names(prm.fit))] <- NULL  

# Calibration parameters:
pf <- read.csv("prm_fit.csv",header=FALSE)
horizon <- pf[pf[,1]=="horizon",2]  
n.MC <- pf[pf[,1]=="n.MC",2]  
n.ABC <- pf[pf[,1]=="n.ABC",2]  
tol.ABC <- pf[pf[,1]=="tol.ABC",2]  
first.time <- pf[pf[,1]=="first.t.stat",2] 
tau.epsilon <- pf[pf[,1]=="tau.epsilon",2] 
multi.core <- pf[pf[,1]=="multi.core",2] 

# Summary statistics definition.
# What kind of summary statistics
# will be used to assess the proximity
# to observed (target) data
#
# Time range where the 
# summary stats are aplied:
prm.stats <- list(first.time = first.time, 
				  last.time = max(obs.data$tb))

# Type of summary stats
# inc.poisson.reg : poisson regression on incidence
# bur.poisson.reg : poisson regression on burials
# inc.max : level and timing of max incidence
# bur.max : level and timing of max burials
stat.type <- list(inc.poisson.reg = TRUE,
				  bur.poisson.reg = TRUE,
				  inc.max = TRUE,
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
					tau.epsilon,
					multi.core,
					do.plot = FALSE
)

t2 <- as.numeric(Sys.time())
message(paste0("ABC fit done in ",round((t2-t1)/60,2)," minutes."))

save.image("sim.RData")

# Visualize posteriors:
true.values <- model.prm.true[names(prm.fit)]
plot.abcfit(post.abc,prm.fit,priors,true.values)

###   Forecast   ###



x <- forecast.fullreport(obs.data,
						 post.abc,
						 prm.fit,
						 prm.model.file = "prm_model.csv",
						 prm.fcast.file = "prm_forecast.csv")

plot.forecast(x, obs.data,obs.data.full)

message("--- END ---")
