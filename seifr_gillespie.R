source("seifr_gillespie_FCT.R")

do.adaptivetau <- TRUE
epsilon <- 0.06 # default<-0.05 ; larger=faster but less accurate

n.MC <- 5   # Monte carlo iterations
horizon <- 180

# Epidemiological parameters
DOL.days <- 5
DOI.days <- 6
funeral.days <- 1
delta <- 0.65

R0 <- 2.5
beta_IS <- R0/(1/DOI.days+1/funeral.days)/2 # <-- quick and dirty (and not correct!)
beta_FS <- 1*beta_IS

# Population size needs to be very large for 
# convergence stochastic -> deterministic
pop.size <- 1E3

# Artificial compartments:
nE <- 3
nI <- 4
nF <- 2

sim <- SEIFR.sim(beta_IS, # contact rate
				 beta_FS,
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
				 nF
				 # 							,   # number of (artificial) compartments for F
				 # 							chgcontact = NULL,  # temporal change in the contact rates (intervention,behaviour,...)
				 # 							t.chgcontact = NULL,  # time when the change starts
				 # 							seed = 1234,
				 # 							do.adaptivetau = TRUE,
				 # 							epsilon = 0.05, # larger=faster but less accurate
				 # 							time.bucket = 1,  # aggregation of incidence in time units (Gillespie events happens any time)
				 # 							remove.fizzles = FALSE  # remove Monte carlo iterations that are fizzles
)



### Plots:

# g <- ggplot(all.sim)+geom_line(aes(x=t,y=prev,colour=factor(mc)))
# plot(g)

# g <- ggplot(subset(sim,mc<=5))+geom_step(aes(x=tb,y=cuminc,colour=factor(mc)),size=3)

g <- ggplot(subset(sim,mc==5))+geom_step(aes(x=tb,y=inc),size=3)
g <- g + geom_step(aes(x=tb,y=buried),size=3)
# plot(g)

z <- subset(sim,mc==5)
plot(z$tb,z$inc,type="l",col="red")
lines(z$tb,z$buried)
