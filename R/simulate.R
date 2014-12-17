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

# Function to make sure values in a vector,a, are between x and y, inclusive
check_values <- function(a,x,y){
	if(max(a)>=x){
		stop("The maximum of 'a' needs to be <=x")
	}
	if(min(a)<=y){
		stop("The minimum of 'a' needs to be >=y")
	}
}


recombine_index <- function(p_trans, L){
	# Simulate recombination events given the transistion prob and the number of loci
	if((length(p_trans)+1)!=L){
		stop("The vector 'p_trans' needs to be of length L-1")
	}
	if(class(p_trans)!="numeric"){
		stop("The vector 'p_trans' needs to be numeric")
	}

	check_values(a=p_trans, x=0, y=0)
	# 1==recombination occurs, 0==no recombination
	out <- sapply(p_trans, function(i) rbinom(n=1, size=1, prob=i)

	)
	return(out)
}





