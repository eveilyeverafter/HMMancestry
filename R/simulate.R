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

recombine_index <- function(p_trans, L){
	# Simulate recombination events given the transistion prob and the number of loci
	if((length(p_trans)+1)!=L){
		stop("The vector 'p_trans' needs to be of length L-1")
	}
	if(class(p_trans)!="numeric"){
		stop("The vector 'p_trans' needs to be numeric")
	}

	if(max(p_trans)>1){
		stop("The maximum transistion probabilty needs to be <=1")
	}
	if(min(p_trans)<0){
		stop("The minimum transistion probabilty needs to be >0")
	}
	# ^ ^ ^ OR run generic check_prob function to be created

	# 1==recombination occurs, 0==no recombination
	out <- sapply(p_trans, function(i) rbinom(n=1, size=1, prob=i)

	)
	return(out)
}





