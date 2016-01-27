filename = "prm_simul.csv"

read.param <- function(filename){
	x <- read.csv(filename,header=FALSE)
	res <- list()
	for (i in 1:nrow(x)){
		res[[i]] <- x[i,2]
	}
	names(res) <- x[,1]
	return(res)
}

prm <- read.param("prm_model.csv")

calc.beta <- function(prm){
	# WARNING: quick and dirty to calculate betas (and not correct!)
	beta_IS <- prm[["R0"]]/(1/prm[["DOI.days"]]+1/prm[["funeral.days"]])/2 
	beta_FS <- prm[["rho"]]*beta_IS
	return(c(list(beta_IS=beta_IS,beta_FS=beta_FS),prm)	)
}