#' Function for simulating recombination 
#'
#' @param p_trans A vector of length L-1 specifying the recombination rate between each SNP position
#' @param L a vector of length 1 specifying the number of SNPs to simulate
#' @param convert A vector specifying which snps are under gene conversion
#' @param convert_rate A vector of length 1 specifying the gene conversion rate


make_parents <- function(L){
	if(length(L)!=1){
		stop("The number of loci to simulate needs to be a vector of length 1")
	}
	if(!inherits(L, "numeric")){ 
		stop("The number of loci to simulate needs to be numeric")
	}
	if(L!=as.integer(L)){
		stop("The number of loci to simulate needs to be an integer")
	}
	if(L<2){
		stop("The number of loci to simulate needs to be >= 2")
	}

	# Code the parent genomes (0 or 1)
	p1 <- rep(0, L)
	p2 <- rep(1, L)

	out <- list(p1=p1, p2=p2)
	class(out) <- c("list", "parent.genomes")
	return(out)

}

# Function to make sure values in a vector, a, are between x and y, inclusive:
check_values <- function(a,x,y){
	if(max(a)>y){
		stop(paste("The maximum of 'a' needs to be <=", y, sep=""))
	}
	if(min(a)<x){
		stop(paste("The minimum of 'a' needs to be >=", x,sep=""))
	}
}

# Simulate recombination events given the vector transistion probs and the number of loci
recombine_index <- function(p_trans, L){
	if((length(p_trans)+1)!=L){
		stop(paste("The vector 'p_trans' needs to be of length ", L-1, sep=""))
	}
	if(class(p_trans)!="numeric"){
		stop("The vector 'p_trans' needs to be numeric")
	}

	check_values(a=p_trans, x=0, y=0.5)
	
	# 1==recombination occurs, 0==no recombination
	out <- sapply(p_trans, function(i) rbinom(n=1, size=1, prob=i)

	)
	return(out)
}

# Functions to simulate mutation for a vector of snps
switch_values <- function(i){
	if(i==0){
		return(1)
	}
	if(i==1){
		return(0)
	}
}

mutate <- function(MLG, mu=mu){
	# MLG <- multilocus genotype
	ll <- length(MLG)
	xx <- MLG
	to.mutate <- sample(x=c(0,1), ll, prob=c(1-mu,mu), replace=TRUE)
	if(sum(to.mutate>0)){
		xx[which(to.mutate==1)] <- sapply(xx[which(to.mutate==1)], switch_values)
	}
	return(xx)
}


recombine <- function(parents, r.index, mu.rate){
	if(!inherits(parents, "parent.genomes")){
		stop(paste("Object ", parents, " needs to be of class parent.genomes", sep=""))
	}
	# S-phase and mutation
	chromotids.orig <- list(p1_1=parents$p1, p1_2=parents$p1, p2_1=parents$p2, p2_2=parents$p2)
	chromotids.mutated <- lapply(chromotids.orig, function(i, ...){
		mutate(i, mu=mu.rate)
		})
	# sample one chromotid from each parent to recombine:
	picked.chromotids <- c(sample(x=c(1,2), size=1, prob=c(0.5, 0.5)), sample(x=c(3,4), size=1, prob=c(0.5, 0.5)))
	chromotids <- data.frame(one=chromotids.mutated[[picked.chromotids[1]]], 
					   two=chromotids.mutated[[picked.chromotids[2]]])
	# recombine these two chromotids:				
	l <- (length(chromotids[[1]]))
	for(i in 1:(l-1)){
		if(i %in% which(r.index==1)){
				one <- chromotids[(i+1):l,1]
				two <- chromotids[(i+1):l,2]
				chromotids[(i+1):l,] <- cbind(two, one)
		}
	}	
	chromotids.recombined <- chromotids.mutated
	chromotids.recombined[[picked.chromotids[1]]] <- chromotids[,1]
	chromotids.recombined[[picked.chromotids[2]]] <- chromotids[,2]

	out <- list(parents=parents, r.index=r.index, mu.rate=mu.rate, 
			chromotids.mutated=chromotids.mutated, picked.chromotids=picked.chromotids, 
			chromotids.recombined=chromotids.recombined)
	class(out) <- c("list", "recombine")
	return(out)	
}


# Simulate sequencing of x coverage across simulated data: 
coverage_sim_bionom <- function(a, p_assign, coverage){
	if(!inherits(a, "recombine")){
		stop(paste("Object 'a' needs to be of class 'recombine'.", sep=""))
	}


	rbinom(n=1, prob=p_assign, size=coverage)

	# output has the number of simulated reads for each parent
}


