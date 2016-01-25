#### 
####  Example file to show how to 
####  call the SEIFR model
####
####  Time unit can be anything, but must be consistent!
####

source("seifr_gillespie.R")

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

# Population size needs to be very large for 
# convergence stochastic -> deterministic
pop.size <- 1E3

# Number of artificial compartments 
# to mimic a gamma distribution
# (the higher 'nX', the narrower the distribution)
nE <- 3
nI <- 4
nF <- 2

sim <- SEIFR.sim(beta_IS, # contact rate I -> S
				 beta_FS, # contact rate F -> S
				 DOL.days,  # avg duration of latency
				 DOI.days,  # avg duration of infectiousness
				 funeral.days, # avg duration of funerals
				 horizon,  # horizon of the simulation
				 n.MC,   # Monte carlo iterations
				 delta , # proportion infected that will die
				 pop.size, 
				 I.init=1,  # initial number of infectious individuals
				 nE,   # number of (artificial) compartments for E
				 nI,   # number of (artificial) compartments for I
				 nF,   # number of (artificial) compartments for F
				 chgcontact = NULL,  # temporal change in the contact rates (intervention,behaviour,...)
				 t.chgcontact = NULL,  # time when the change starts
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
seed  <-  1234

sim2 <- reporting.filter(sim,
						 report.inc.prob, 
						 report.inc.lag.mean, 
						 report.inc.lag.var, 
						 report.bur.prob, 
						 report.bur.lag.mean,
						 report.bur.lag.var,
						 seed)



### Plots just one Monte Carlo iteration

mc.chosen <- n.MC  # Choose any one

g <- ggplot(subset(sim2,mc==mc.chosen))+geom_step(aes(x=tb,y=inc),size=2)
g <- g + geom_step(aes(x=tb,y=buried),size=2,colour="red")
plot(g)

