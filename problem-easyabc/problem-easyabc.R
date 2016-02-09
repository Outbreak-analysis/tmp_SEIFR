library(EasyABC)


# Define a dummy model that is supposed
# to mimic a stochastic simulation and the calculation 
# of summary statistics (two in this case)
toy.model.1 <- function(x){
	set.seed(x[1]) # so that each core is initialized with a different seed value.
	return(c(x[2]+x[3]+rnorm(1,0,0.1),
			 x[2]*x[3]+rnorm(1,0,0.1)))
}

# prior distributions of parameter to infer:
toy_prior <- list(c("unif",0,1),c("normal",1,2))

# Mimic what we would get
# if we would calculate the
# summary stats on observed data:
x.true <- 0.8
y.true <- 1.7
summ.stat.obs <- c(x.true+y.true, x.true*y.true)


# ABC parameters:
n.ABC <- 100
tol.ABC <- 0.1  # proportion of best samples retained for posterior distribution

# number of cores used:
n_cluster <- 4

# ABC fit 
fit1 <- ABC_rejection(model = toy.model.1,
					  prior = toy_prior,
					  nb_simul = n.ABC,
					  summary_stat_target = summ.stat.obs,
					  tol = tol.ABC,
					  use_seed = TRUE,
					  n_cluster = n_cluster)	
# All works well:
xx<-fit1$param[,1]
yy<-fit1$param[,2]
plot(xx,yy)
points(x.true,y.true,pch=16,cex=3,col="red")



# Now, we want the "toy.model" function to
# wrap a very large function (representing a complicated model)
# and use EasyABC's multi-cores capabilities. 
# But 'EasyABC' forces to have the wrapping function
# that will be called as the argument 'model' of
# 'ABC_rejection' to be self-contained (or coded
# outside R as a binary).

# That may not be convenient, especially if one
# already has existing code that needs to be
# 'plugged' with EsayABC (and you don't want 
# to rewrite everything)

# If we simply call a function defined outside the 
# model wrapper, we get an error saying that the function
# used outside the wrapper is unknown.
# This is because EasyABC does not export functions from
# the global environment before making the clusters 
# for the parallel environments.
#
# To fix this problem, I "re-source" locally the source code
# of 'EasyABC' after tweaking it. The hack consists 
# of including a line that export all functions in
# the clusters that will be used for parallel computing.
# (look for the string "DC HACK"in the locally copied R sources)



huge.fct <- function(z){
	### This represent a function 
	### with 1000s lines of R code
	### ...
	### ...
	
	return(z+1)
}

toy.model.2 <- function(x){
	# We want this wrapping function to
	# be relatively compact, hence we
	# call the huge function (not recoded inside toy.model.2)
	
	set.seed(x[1]) 
	
	# We call the huge function:
	H <- huge.fct(x[4])
	
	return(c(H + x[2]+x[3]+rnorm(1,0,0.1),
			 x[2]*x[3]+rnorm(1,0,0.1)))
}

# just making sure it works in the global environment:
# toy.model.2(c(1234, 1,2,3))

# definig the priors again (there's one more parameter)
toy_prior.2 <- list(c("unif",0,1),
					c("normal",1,2),
					c("unif",1,1))


# ABC fit 
message("starting fit2...")

# Swicth DC's hack on and off:
do.DC_HACK <- TRUE

if(do.DC_HACK){
	library(parallel)
	source("EasyABC-internal-DC_HACK.R")
	source("ABC_rejection-DC_HACK.R")
}


fit2 <- ABC_rejection(model = toy.model.2,
					  prior = toy_prior.2,
					  nb_simul = n.ABC,
					  summary_stat_target = summ.stat.obs,
					  tol = tol.ABC,
					  use_seed = TRUE,
					  n_cluster = 2)	

xx<-fit2$param[,1]
yy<-fit2$param[,2]
plot(xx,yy,main="with huge function and DC HACK")
points(x.true,y.true,pch=16,cex=3,col="red")
message("fit2 SUCCESS!")

