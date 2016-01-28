#### 
####  Example file to show how to 
####  call the SEIFR model
####
####  Time unit can be anything, but must be consistent!
####

source("seifr_gillespie.R")
source("fit_abc.R")
source("utils.R")

# Read parameters from files:
simul.prm <- read.param("prm_simul.csv")
model.prm <- calc.beta(read.param("prm_model.csv"))

# Simulate epidemics:
sim <- SEIFR.sim(model.prm, simul.prm, seed = 1234)

# Apply reporting layer:
report.prm <- read.param("prm_report.csv")
sim2 <- reporting.filter(sim,prm=report.prm)

### Plots just one Monte Carlo iteration
if (FALSE){
	mc.chosen <- 1  # Choose any one
	g <- ggplot(subset(sim2,mc==mc.chosen)) + geom_step(aes(x=tb,y=inc),size=2)
	g <- g + geom_step(aes(x=tb,y=buried),size=2,colour="red")
	plot(g)
}


###
###    FAKE CALIBRATION WITH 'ABC'
###

# We pretend that the first simulation data set 'mc=1'
# is an actual observation:
hz <- 23
obs.data <- subset(sim2,mc==1 & tb<=hz)

plot(obs.data$tb, obs.data$inc, typ="s",lwd=3)
lines(obs.data$tb, obs.data$bur, typ="s",col="red",lwd=3)

# Specify which parameters will be calibrated
# and their prior distributions:
prm.fit <- list(beta_IS=0.999, beta_FS=0.999)  # <- values do not matter here
priors <-  list(c("unif",0.05,3),
				c("unif",0.05,3))

# Fixed parameters (will NOT be calibrated):
# (take the original model params and remove the ones that will be)
prm.fixed  <- model.prm
prm.fixed[which(names(prm.fixed) %in% names(prm.fit))] <- NULL  

horizon <- hz+20  
n.ABC <- 100
tol.ABC <- 0.20

# Summary statistics definition
prm.stats <- list(first.time = 6, 
				  last.time = hz)

stat.type <- list(inc.poisson.reg=TRUE,
				  bur.poisson.reg=TRUE,
				  inc.max = TRUE,
				  bur.max = FALSE)

# Swicth DC's hack on and off:
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
					n.ABC,
					tol.ABC,
					multi.core = 0 
)

t2 <- as.numeric(Sys.time())
print(paste0("ABC fit done in ",round((t2-t1)/60,2)," minutes."))

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