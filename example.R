#### 
####  Example file to show how to 
####  call the SEIFR model
####
####  Time unit can be anything, but must be consistent!
####

source("seifr_gillespie.R")
source("fit_abc.R")

do.adaptivetau <- TRUE   # FALSE=true Gillespie ; TRUE=approximation of Gillespies algo
epsilon <- 0.06 # default<-0.05 ; larger=faster but less accurate

n.MC <- 5   # Monte carlo iterations
horizon <- 180  

# Epidemiological parameters
DOL.days <- 5   # Mean duration of latency
DOI.days <- 6   # Mean duration of infectiousness
funeral.days <- 1 # Mean duration of funerals
delta <- 0.65  # Death rate following infection

R0 <- 2.5   # Not the correct R0, but shouldn't be far
beta_IS <- R0/(1/DOI.days+1/funeral.days)/2 # <-- quick and dirty (and not correct!)
beta_FS <- 1*beta_IS

print(paste0("beta_IS=",beta_IS))
print(paste0("beta_FS=",beta_FS))

# Population size needs to be very large for 
# convergence stochastic -> deterministic
pop.size <- 1E3

# Number of artificial compartments 
# to mimic a gamma distribution
# (the higher 'nX', the narrower the distribution)
nE <- 3
nI <- 4
nF <- 2

model.prm <- list(beta_IS=beta_IS, # contact rate I -> S
				  beta_FS=beta_FS, # contact rate F -> S
				  DOL.days=DOL.days,  # avg duration of latency
				  DOI.days=DOI.days,  # avg duration of infectiousness
				  funeral.days=funeral.days, # avg duration of funerals
				  delta=delta , # proportion infected that will die
				  pop.size=pop.size, 
				  I.init=1,  # initial number of infectious individuals
				  nE=nE,   # number of (artificial) compartments for E
				  nI=nI,   # number of (artificial) compartments for I
				  nF=nF,   # number of (artificial) compartments for F
				  chgcontact = NULL,  # temporal change in the contact rates (intervention,behaviour,...)
				  t.chgcontact = NULL  # time when the change starts
)


sim <- SEIFR.sim(model.prm,
				 horizon,  # horizon of the simulation
				 n.MC,   # Monte carlo iterations
				 seed = 1234,
				 do.adaptivetau = do.adaptivetau,
				 epsilon = epsilon, # larger=faster but less accurate
				 time.bucket = 1,  # aggregation of incidence in time units (Gillespie events happens any time)
				 remove.fizzles = FALSE  # remove Monte carlo iterations that are fizzles
)


# Apply reporting layer:

report.inc.prob <- 0.6
report.inc.lag.mean <- 4 
report.inc.lag.var <- 0.05 
report.bur.prob <- 0.8 
report.bur.lag.mean <- 3
report.bur.lag.var <- 0.02

sim2 <- reporting.filter(sim,
						 report.inc.prob, 
						 report.inc.lag.mean, 
						 report.inc.lag.var, 
						 report.bur.prob, 
						 report.bur.lag.mean,
						 report.bur.lag.var)



### Plots just one Monte Carlo iteration

mc.chosen <- 1  # Choose any one

g <- ggplot(subset(sim2,mc==mc.chosen))+geom_step(aes(x=tb,y=inc),size=2)
g <- g + geom_step(aes(x=tb,y=buried),size=2,colour="red")
plot(g)


###
###    FAKE CALIBRATION WITH 'ABC'
###

# We pretend that a simulation data set
# is an actual observation:
hz <- 23
obs.data <- subset(sim,mc==1 & tb<=hz)

plot(obs.data$tb, obs.data$inc, typ="s",lwd=3)
lines(obs.data$tb, obs.data$bur, typ="s",col="red",lwd=3)

# Specify which parameters will be calibrated:
prm.fit <- list(beta_IS=0.999, beta_FS=0.999)
priors <-  list(c("unif",0.05,2),
				c("unif",0.05,2))

# Fixed parameters (will NOT be calibrated):
prm.fixed  <- list(DOL.days=DOL.days,  # avg duration of latency
				   DOI.days=DOI.days,  # avg duration of infectiousness
				   funeral.days=funeral.days, # avg duration of funerals
				   delta=delta , # proportion infected that will die
				   pop.size=pop.size, 
				   I.init = 1,  # initial number of infectious individuals
				   nE=nE,   # number of (artificial) compartments for E
				   nI=nI,   # number of (artificial) compartments for I
				   nF=nF,   # number of (artificial) compartments for F
				   chgcontact = NULL,  # temporal change in the contact rates (intervention,behaviour,...)
				   t.chgcontact = NULL  # time when the change starts
)

horizon <- hz+20  
n.ABC <- 1000
tol.ABC <- 0.05

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
					multi.core = TRUE # <-- TRUE does not work!
					)

t2 <- as.numeric(Sys.time())
print(paste0("ABC fit done in ",round((t2-t1)/60,2)," minutes."))

# Visualize posteriors:
par(mfrow=c(1,3))
plot(post.abc$param[,1],post.abc$param[,2])
abline(v=mean(post.abc$param[,1]))
abline(h=mean(post.abc$param[,2]))
points(beta_IS, beta_FS,col="red",cex=5,pch=16)
hist(post.abc$param[,1],breaks=12,col="grey")
abline(v=beta_IS,col="red")
hist(post.abc$param[,2],breaks=12,col="grey")
abline(v=beta_FS,col="red")