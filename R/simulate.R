#' @title Simulate haploid parental genotypes along a chromosome
#' 
#' @description This function creates divergent genotypes between two parents. The first parent's genotypes
#' are denoted by a vector of zeros of length \code{L} while the second parent's genotypes are a vector of ones.
#' 
#' @param L is an integer specifying the number of loci to simulate. L must be greater than 1. 
#' 
#' @return an object of class \code{parent.genomes} giving the parental snp calls at each of the \code{L} loci. 
#' 
#' @references none.
#' 
#' @seealso \code{\link{recombine_index}}, \code{\link{recobmine}}
#' 
#' @author Tyler D. Hether
#' 
#' @export make_parents
#' 
#' @examples
#' 
#' # Make a chromosome of L=50 loci for each parent.
#' make_parents(L=50)

make_parents <- function(snps){
	# if(length(L)!=1){
	# 	stop("The number of loci to simulate needs to be a vector of length 1")
	# }
	# if(!inherits(L, "numeric")){ 
	# 	stop("The number of loci to simulate needs to be numeric")
	# }
	# if(L!=as.integer(L)){
	# 	stop("The number of loci to simulate needs to be an integer")
	# }
	# if(L<2){
	# 	stop("The number of loci to simulate needs to be >= 2")
	# }

	# Code the parent genomes (0 or 1)
	p1 <- rep(0, length(snps))
	p2 <- rep(1, length(snps))

	out <- list(snps=snps, p1=p1, p2=p2)
	class(out) <- c("list", "parent.genomes")
	return(out)

}

#' @title Simulate recombination points along a chromosome.
#' 
#' @description Simulate recombination events given the vector transistion 
#' probs and the number of loci.
#' 
#' @param p.trans a vector of length \code{L-1} specifying the recombination 
#' rate between two snps. Each element in \code{p.trans} must >= 0 (no recombination) 
#' and <= 0.5 (independent assortment).
#' 
#' @return a vector of length \code{L-1} specifying the location of recombiantion points (i.e., between snps).
#' 
#' @seealso \code{\link{make_parents}}, \code{\link{fitEcm}}
#' 
#' @author Tyler D. Hether
#' 
#' @export recombine_index
#' 
#' @examples
#' 
#' set.seed(1234567)
#' l <- 100000 # number of loci to simulate
#' scale <- 0.01 # d = bp*scale
#' snps <- 1:l # floor(seq(from=1, to=1000, length.out=l)) # position of snps (in bps)
#' # snps[25:l] <- 5000+snps[25:l]
#' sum(recombine_index(scale=scale, snps=snps))/l

recombine_index <- function(scale, snps){
	# if(class(p.trans)!="numeric"){
	# 	stop("The vector 'p.trans' needs to be numeric")
	# }

	# check_values(a=p.trans, x=0, y=0.5)
	haldane <- function(d){ (1/2)*(1- exp(-2*d))}
	displacement <- sapply(2:length(snps), function(x){
		return(snps[x]-snps[x-1])
		})
	# 1==recombination occurs, 0==no recombination
	out <- sapply(1:length(displacement), function(i) rbinom(n=1, size=1, prob=haldane(displacement[i]*scale)))

	return(out)
}


#' @title Simulate recombination
#' 
#' @description This function simulates recombination between two parent genomes of class
#' \code{parent.genomes}.
#' 
#' @param parents an object of class \code{parent.genomes} specifying the parental genotypes 
#' at each simulated snp.
#' 
#' @param r.index a vector of length \code{L-1} specifying whether a recombination is to be 
#' simulated (1) or not (0) in between two adjacent snps. 
#' 
#' @param mu.rate a numeric between 0 and 1 (inclusive) specifying the per snp mutation rate. 
#' 
#' @param f.cross a numeric between 0 and 1 (inclusive) giving the frequency of recombination 
#' events that result in crossing over. This is same as 1 minus the frequenc of non-crossovers.
#' 
#' @param f.convert a numeric between 0 and 1 (inclusive) that gives the frequency of gene conversion 
#' during recombination.
#' 
#' @param length.conversion an integer specifying the mean (and variance) of a given gene conversion 
#' tract.
#' 
#' @details (to do)
#' 
#' @return an object of class \code{recombine} that contains the following data:
#' \describe{
#' 	\item{parents}{input data}
#' 	\item{r.index}{input data}
#' 	\item{m.rate}{input data}
#' 	\item{f.cross}{input data}
#' 	\item{f.convert}{input data}
#' 	\item{length.conversion}{input data}
#' 	\item{chromatids.mutated_snps}{a list giving the genotypes of 4 chromatids (two pairs of sister chromatids) following mutation}
#' 	\item{chromatids.recombined}{a list giving the genotypes of 4 chromatids following recombination}
#' }
#' 
#' @references (to do: refernce a crossover vs non-crossover paper)
#' 
#' @seealso \code{\link{parent_genomes}}, \code{\link{recombine_index}}
#' 
#' @author Tyler D. Hether
#' 
#' @export recombine
#' 
#' @examples
#' set.seed(1234567) # For reproducability
#' l <- 50 # number of loci to simulate
#' scale <- 0.01; snps <- 1:l;
#' r <- recombine_index(scale, snps) 
#' p_a <- .999 # probability of correct sequencing assignment (1-sequence error rate)
#' p <- make_parents(l) # make the parents
#' recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=0, length.conversion=10) # recombine parents
#' recomb_sim

recombine <- function(parents, r.index, mu.rate=0, f.cross=0.5, f.convert=0, length.conversion=20){
	if(!inherits(parents, "parent.genomes")){
		stop(paste("Object ", parents, " needs to be of class parent.genomes", sep=""))
	}
	l <- length(parents$p1)

	# S-phase and mutation
	chromatids.orig <- list(p1_1=parents$p1, p1_2=parents$p1, p2_1=parents$p2, p2_2=parents$p2)
	chromatids.mutated_snps <- lapply(chromatids.orig, function(i, ...){
		mutate_snps(i, mu=mu.rate)
		})
	
	# Recombine chromatids:
	chromatids.recombined <- chromatids.mutated_snps
	for(i in 1:(l-1)){
		# Cases where there is a recombination event, pick
		# two of the four chromatids to recombine:
		if(i %in% which(r.index==1)){
			picked.chromatids <- sample(c(1:4), 2, replace=FALSE)
			# Recombine them:
			chromatids <- data.frame(one=chromatids.recombined[[picked.chromatids[1]]], 
					   two=chromatids.recombined[[picked.chromatids[2]]])
			one <- chromatids[(i+1):l,1]
			two <- chromatids[(i+1):l,2]
				# Does recombination result in a crossover (1) or non-crossover (0)?
				to_cross <- cross_vs_noncross(f.cross=f.cross)
				if(to_cross==1){
					chromatids[(i+1):l,] <- cbind(two, one)
				}			
				# Is there gene conversion at this recombination point?
				if(sample(c(1,0),1,prob=c(f.convert,1-f.convert))==1){
					# simulate gene conversion on one of the non-picked chromatids:
					to_convert <- sample(c(1:4)[-picked.chromatids],1) # non-picked chromotid
					length_of_conversion <- rpois(n=1,lambda=length.conversion)
					# if the length of the tract extends pass the chromosome then stop at the end of the chromosome.
					# Otherwise, create a gene conversion tract of length length_of_conversion starting at the 
					# recombination point. 
					if( (i+1+length_of_conversion) > l){
						chromatids.recombined[[to_convert]][c((i+1):l)] <- sapply(chromatids.recombined[[to_convert]][c((i+1):l)], function(j){switch_values(j)})
					} else {
						chromatids.recombined[[to_convert]][c((i+1):(i+1+length_of_conversion))] <- sapply(chromatids.recombined[[to_convert]][c((i+1):(i+1+length_of_conversion))], function(j){switch_values(j)})
					}

				}
			# Save results:
			chromatids.recombined[[picked.chromatids[1]]] <- chromatids[,1]
			chromatids.recombined[[picked.chromatids[2]]] <- chromatids[,2]
		}
	}

	out <- list(parents=parents, r.index=r.index, mu.rate=mu.rate, f.cross=f.cross, f.convert=f.convert, length.conversion=length.conversion,
			chromatids.mutated_snps=chromatids.mutated_snps, chromatids.recombined=chromatids.recombined)
	class(out) <- c("list", "recombine")
	return(out)	
}

#' @title Simulate sequencing across simulated data:
#' 
#' @description This function simulates sequecing across simulated data of class \code{recombine}.
#' 
#' @param simdata an object of class \code{recombine}
#' 
#' @param p.assign a numeric between 0 and 1 (inclusive) that gives the probability of 
#' correct sequencing assignment.
#' 
#' @param coverage the mean (and variance) of sequencing coverage to be simlated.  \code{coverage} is 
#' sampled from a poisson distribution (i.e., lambda=\code{coverage})
#' 
#' @return an object of class \code{snp.recom} providing the simulated sequencing reads for each parent.
#' 
#' @seealso \code{\link{recombine}}, \code{\link{estimate_anc_fwd_back}}
#' 
#' @author Tyler D. Hether 
#' 
#' @export simulate_coverage
#' 
#' @examples
#' set.seed(1234567) # For reproducability
#' l <- 1000 # number of loci to simulate
#' rec <- 0.01 # recombination rate between each snp
#' r <- recombine_index(rep(rec, l-1)) # recombination rate between each snp (vector form)
#' p_a <- .999 # probability of correct sequencing assignment (1-sequence error rate)
#' p <- make_parents(l) # make the parent
#' recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=1, length.conversion=10) # recombine parents
#' sim_reads <- simulate_coverage(simdata=recomb_sim, p.assign=p_a, coverage=1) # simulate sequencing coverage

simulate_coverage <- function(simdata, p.assign, coverage){
	if(!inherits(simdata, "recombine")){
		stop(paste("Object 'simdata' needs to be of class 'recombine'.", sep=""))
	}

	simulated.reads.total <- lapply(simdata$chromatids.recombined, function(i) return(rpois(i, coverage)))
	
	# This block returns the number of simulated reads for each parental type, 0 and 1, for each snp position.
	out <- lapply(1:4, function(i){
		correct_reads <- rbinom(size=simulated.reads.total[[i]], prob=p.assign, n=simulated.reads.total[[i]])
		incorrect_reads <- simulated.reads.total[[i]] - correct_reads 

		p0.assign <- sapply(1:length(simdata$chromatids.recombined[[i]]), function(x){
				if(simdata$chromatids.recombined[[i]][x]==0){
					return(correct_reads[x])
				}
				if(simdata$chromatids.recombined[[i]][x]==1){
					return(incorrect_reads[x])
				}
			})
		p1.assign <- sapply(1:length(simdata$chromatids.recombined[[i]]), function(x){
				if(simdata$chromatids.recombined[[i]][x]==1){
					return(correct_reads[x])
				}
				if(simdata$chromatids.recombined[[i]][x]==0){
					return(incorrect_reads[x])
				}
			})
		return(list(p0.assign=p0.assign, p1.assign=p1.assign, snps=simdata$parents$snps, states_given=simdata$chromatids.recombined[[i]]))

	})
	out[['snps']] <- simdata$parents$snps

	class(out) <- c("list", "snp.recom")
	return(out)
}

#' @title Simulate recombination of a given number of tetrads
#' 
#' @description This is a wrapper function of many other functions that simulates a given number of
#' tetrads each with four haploid spores that are recombinants between two parents.  
#' 
#' @param n.tetrads An integer specifying the number of tetrads to simulate.
#' 
#' @param l an integer describing the number of loci to simulate.
#' 
#' @param rec a vector of length \code{l-1} that contains the the recombation rate 
#' between each snp.
#' 
#' @param p.assign a numeric between 0 and 1 (inclusive) that gives the probability of 
#' correct sequencing assignment.
#' 
#' @param mu.rate a numeric between 0 and 1 (inclusive) specifying the per snp mutation rate. 
#' 
#' @param f.cross a numeric between 0 and 1 (inclusive) giving the frequency of recombination 
#' events that result in crossing over. This is same as 1 minus the frequenc of non-crossovers.
#' 
#' @param f.convert a numeric between 0 and 1 (inclusive) that gives the frequency of gene conversion 
#' during recombination.
#' 
#' @param length.conversion an integer specifying the mean (and variance) of a given gene conversion 
#' tract.
#' 
#' @param coverage the mean (and variance) of sequencing coverage to be simlated.  \code{coverage} is 
#' sampled from a poisson distribution (i.e., lambda=\code{coverage})
#' 
#' @return A list of length n.tetrads of class \code{tetrad}. Each element of 
#' \code{tetrad} is itself a list of class \code{individual.tetrad}, which 
#' contains four lists, one for each spore, and three elements: 
#' \describe{
#' 	\item{p0.assign}{The number of reads that were simulated for parent 0.}
#' 	\item{p1.assign}{The number of reads taht were simulated for parent 1.}
#' 	\item{snps}{The snp id along the simulated chromosome.}
#' } 
#' 
#' @seealso \code{\link{recombine}}, \code{\link{make_parents}}, 
#' \code{\link{simulate_coverage}}, \code{\link{recombine_index}}, \code{\link{id_hotspots}}
#' 
#' @author Tyler D. Hether 
#' 
#' @export sim_tetrad
#' 
#' @examples
#' set.seed(1234567) # For reproducability
#' # simulate a recombination hotspot between the 99th and 100th snp
#' l=200; scale = 0.01; snps=2*(1:l);
#' n.tetrads <- 100 # number of spores to simulate
#' res <- sim_tetrad(n.tetrads=n.tetrads, l=l, scale=scale, snps=snps, p.assign=0.999, 
#'    mu.rate=0, f.cross=0.8, f.convert=0.9, length.conversion=10, coverage=2.5)

sim_tetrad <- function(n.tetrads, l, scale, snps, p.assign, mu.rate, f.cross, f.convert, 
        length.conversion, coverage){

    out <- lapply(1:n.tetrads, function(Z, ...){

        r <- recombine_index(scale, snps)
        p <- make_parents(snps)
        recomb_sim <- recombine(parents=p, r.index=r, mu.rate=mu.rate, f.cross=f.cross, 
                f.convert=f.convert, length.conversion=length.conversion)
        sim_reads <- simulate_coverage(simdata=recomb_sim, p.assign=p.assign, coverage=coverage)
        class(sim_reads) <- list("individual.tetrad")
        return(sim_reads)
    })
    class(out) <- c("list", "tetrad")
    return(out)
}

#' @title Simulate random spores en masse
#' 
#' @description This is a wrapper function of many other functions that simulates a given number of
#' haploids each recombinant between two parents.  
#' 
#' @param n.spores An integer specifying the number of spores to simulate.
#' 
#' @param l an integer describing the number of loci to simulate.
#' 
#' @param rec a vector of length \code{l-1} that contains the the recombation rate 
#' between each snp.
#' 
#' @param p.assign a numeric between 0 and 1 (inclusive) that gives the probability of 
#' correct sequencing assignment.
#' 
#' @param mu.rate a numeric between 0 and 1 (inclusive) specifying the per snp mutation rate. 
#' 
#' @param f.cross a numeric between 0 and 1 (inclusive) giving the frequency of recombination 
#' events that result in crossing over. This is same as 1 minus the frequenc of non-crossovers.
#' 
#' @param f.convert a numeric between 0 and 1 (inclusive) that gives the frequency of gene conversion 
#' during recombination.
#' 
#' @param length.conversion an integer specifying the mean (and variance) of a given gene conversion 
#' tract.
#' 
#' @param coverage the mean (and variance) of sequencing coverage to be simlated.  \code{coverage} is 
#' sampled from a poisson distribution (i.e., lambda=\code{coverage})
#' 
#' @return A list of length n.spores. Each element of the list is of class \code{single.spore}, 
#' which contains three elements: 
#' \describe{
#' 	\item{p0.assign}{The number of reads that were simulated for parent 0.}
#' 	\item{p1.assign}{The number of reads taht were simulated for parent 1.}
#' 	\item{snps}{The snp id along the simulated chromosome.}
#' } 
#' 
#' @seealso \code{\link{recombine}}, \code{\link{make_parents}}, 
#' \code{\link{simulate_coverage}}, \code{\link{recombine_index}}, \code{\link{id_hotspots}}
#' 
#' @author Tyler D. Hether 
#' 
#' @export sim_en_masse
#' 
#' @examples
#' set.seed(1234567) # For reproducability
#' # simulate a recombination hotspot between the 100th and 101st snp
#' l=100; scale = 0.01; snps=1:l
#' n.spores <- 500 # number of spores to simulate
#' spores <- sim_en_masse(n.spores=n.spores, l=l, scale=scale, snps=snps, 
#' p.assign=.999, mu.rate=0.001, f.cross=0.5, 
#'     f.convert=0.5, length.conversion=10, coverage=1)

sim_en_masse <- function(n.spores, l, scale, snps, p.assign, mu.rate, f.cross, f.convert, 
        length.conversion, coverage){

    out <- lapply(1:n.spores, function(Z, ...){

        r <- recombine_index(scale, snps)
        p <- make_parents(l)
        recomb_sim <- recombine(parents=p, r.index=r, mu.rate=mu.rate, f.cross=f.cross, 
                f.convert=f.convert, length.conversion=length.conversion)
        sim_reads <- simulate_coverage(simdata=recomb_sim, p.assign=p.assign, coverage=coverage)
        to.pick <- sample(c(1:4), 1)
        class(sim_reads[[to.pick]]) <- list("single.spore")
        return(sim_reads[[to.pick]])
    })
    class(out) <- c("list", "en.masse")
    return(out)
}

# Minor functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to make sure values in a vector, a, are between x and y, inclusive:
check_values <- function(a,x,y){
	if(max(a)>y){
		stop(paste("The maximum of 'a' needs to be <", y, sep=""))
	}
	if(min(a)<x){
		stop(paste("The minimum of 'a' needs to be >", x,sep=""))
	}
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

mutate_snps <- function(xx, mu=mu){
	# MLG <- multilocus genotype
	ll <- length(xx)
	to.mutate_snps <- sample(x=c(0,1), ll, prob=c(1-mu,mu), replace=TRUE)
	if(sum(to.mutate_snps>0)){
		xx[which(to.mutate_snps==1)] <- sapply(xx[which(to.mutate_snps==1)], switch_values)
	}
	return(xx)
}

# Given the frequency of crossover events, f.cross, sample whether one is to occur 
#     (1=yes, crossover; 0=no, non-crossover).
cross_vs_noncross <- function(f.cross){
	sample(c(1,0), 1, prob=c(f.cross, 1-f.cross))
}


