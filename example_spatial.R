#### 
####  Example file to show how to 
####  call the spatial SEIFR model 
####
####  Time unit can be anything, but must be consistent!
####

source("seifr_gillespie.R")
source("fit_abc.R")
source("utils.R")
source("forecast.R")

pdf.options(width=15,height=10)

# Number of locations:
n.loc <- 2

# Read parameters from files:
simul.prm <- read.param("prm_simul.csv")
# Reporting parameters
report.prm <- read.param("prm_report.csv")
# Model parameter for each location
model.prm <- list()
model.prm.true <- list()
model.prm.filename <- list()

for (i in 1:n.loc) {
	model.prm.filename[[i]] <- paste0("prm_model_loc_",i,".csv")
	model.prm[[i]] <- calc.beta(read.param(model.prm.filename[[i]]))	
	model.prm.true[[i]] <- model.prm[[i]]
}

t1 <- as.numeric(Sys.time())

# Simulate epidemics:
sim <- list()
for (i in 1:n.loc) {
	# simulate epidemic at location #i
	sim.tmp <- SEIFR.sim(model.prm[[i]], simul.prm, seed = 1234)
	# Apply reporting layer:
	sim[[i]] <- reporting.filter(sim.tmp, prm=report.prm)
	sim[[i]]$loc <- i
}
# merge everything in one data frame:
allsim <- do.call("rbind",sim)

t2 <- as.numeric(Sys.time())

message(paste("Synthetic data simulations:",round((t2-t1)/60,2),"minutes"))

g <- ggplot(allsim)+geom_line(aes(x=tb,y=inc,colour=factor(mc)))+facet_wrap(~loc)
g <- g + scale_y_log10() + annotation_logticks(sides = "lr")
plot(g)

###
###    FAKE CALIBRATION WITH 'ABC'
###

# We pretend that one simulation data set 'mc'
# is an actual observation until a specified horizon ('hz'):
hz <- 19

obs.data.full <- list()
obs.data <- list()

for (i in 1:n.loc) {
	obs.data.full[[i]] <- subset(allsim, loc==i & mc==2)
	obs.data[[i]] <- subset(obs.data.full[[i]], tb <= hz)
}

all.obs.data.full <- do.call("rbind",obs.data.full)
all.obs.data <- do.call("rbind",obs.data)

# How the data look like:
# By locations:
g <- ggplot(all.obs.data.full)+geom_step(aes(x=tb,y=inc),colour="grey")
g <- g + geom_step(aes(x=tb,y=buried),colour="pink")
g <- g + geom_step(data = all.obs.data,
				   aes(x=tb,y=inc), size=1.5)
g <- g + geom_step(data = all.obs.data,
				   aes(x=tb,y=buried), size=1.5, colour="red")
g <- g + facet_wrap(~loc,scales = 'free_y')
g <- g + scale_y_log10()
g <- g + ggtitle("Observed data by location")
plot(g)
# Overall (sum of all locations)
obs.sum <- ddply(all.obs.data.full,c("tb"),summarize,
				 inc.sum = sum(inc),
				 bur.sum = sum(buried))
plot(obs.sum$tb,obs.sum$inc.sum,
	 typ="s",
	 main = "Global epidemic \n synthetic data",
	 xlab="time",ylab="count")
lines(obs.sum$tb,obs.sum$bur.sum,col="red")

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
prm.fixed <- list()
for (i in 1:n.loc) {
	prm.fixed[[i]]  <- model.prm[[i]]
	prm.fixed[[i]][which(names(prm.fixed[[i]]) %in% names(prm.fit))] <- NULL  
}

# Calibration parameters:
pf <- read.csv("prm_fit_spatial.csv",header=FALSE)
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
				  last.time = max(all.obs.data$tb))

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

post.abc <- list()
for (i in 1:n.loc) {
	message(paste("fitting loc",i,"..."))
	post.abc[[i]] <- fit.abc(prm.fit, 
							 prm.fixed[[i]], 
							 obs.data[[i]],
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
}

t2 <- as.numeric(Sys.time())
message(paste0("ABC fit done in ",round((t2-t1)/60,2)," minutes."))

save.image("fit_loc.RData")

# Visualize posteriors:
true.values <- list()
for(i in 1:n.loc) {
	true.values[[i]] <- model.prm.true[[i]][names(prm.fit)]
	plot.abcfit(post.abc[[i]], prm.fit, priors, true.values[[i]],
				location = i)
}

###   Forecast   ###

fcast <- list()
for(i in 1:n.loc) {
	message(paste("forecasting loc",i,"..."))
	fcast[[i]] <- forecast.fullreport(obs.data[[i]],
									  post.abc[[i]],
									  prm.fit,
									  prm.model.file = model.prm.filename[[i]],
									  prm.fcast.file = "prm_forecast.csv")
	plot.forecast(fcast[[i]], obs.data[[i]],obs.data.full[[i]])
}

message("--- END ---")
