#' @title Simulate haploid parental genotypes along a chromosome
#' 
#' @description This function creates divergent genotypes between two parents. The first parent's genotypes
#' are denoted by a vector of zeros of length \code{L} while the second parent's genotypes are a vector of ones.
#' 
#' @param \code{snps} a vector gving the locations of snps along a chromosome
#' 
#' @return an object of class \code{parent.genomes} giving the parental ancestry at each snp. 
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
#' make_parents(snps=c(1:10)*10)

make_parents <- function(snps){
	if(length(snps)<2){
		stop("The number of loci to simulate needs to be >= 2")
	}

	# Code the parent genomes (0 or 1)
	p1 <- rep(0, length(snps))
	p2 <- rep(1, length(snps))

	out <- list(snps=snps, p1=p1, p2=p2)
	class(out) <- c("list", "parent.genomes")
	return(out)
}

#' @title Simulate recombination points along a chromosome.
#' 
#' @description Simulate recombination events given the vector of snp locations and recombination rate(s). 
#' Haldane's function is used to simulate the probabilty of getting an odd number of crossover events between
#' any two snps given their physical distance and the specified recombination rate.
#' 
#' @param \code{scale} either vector of length 1 specifying the genome-wide recombination rate (Morgans/bp)
#' or a vector of length \code{snps-1} specifying the recombination rate between all snps. In either
#' case rates must be positive.
#'
#' @param \code{snps} a vector gving the locations of snps along a chromosome
#' 
#' @return a vector of length \code{snps-1} specifying the location of recombiantion points.
#' 
#' @seealso \code{\link{make_parents}}
#' 
#' @author Tyler D. Hether
#' 
#' @export recombine_index
#' 
#' @examples
#' 
#' set.seed(1234567)
#' # Simulating recombination across 1000 loci randomly spaced 
#' # along a 200kbp chromosome
#' chromSize <- 2e5 		# Chromosome size
#' l <- 1e3 				# number of loci to simulate
#' c <- 1e-04 				# Morgans/bp
#' snps <- sort(sample(size=l, 1:chromSize, replace=FALSE))
#' # Find recombination points between all the snps
#' recomb_points <- recombine_index(scale=c, snps=snps)
#' recomb_points

recombine_index <- function(scale, snps){

	if(length(scale)==1){
		c <- rep(scale, length(snps)-1)
	} else {
		c <- scale
	}
	if(!all(scale >=0)){
		stop("all elements in scale need to be a positive rate (Morgans / bp)")
	}		
	if(length(snps)<2){
		stop("the vector of snps needs to be greater than 1 in length")
	}
	if(length(c)!=(length(snps)-1)){
		stop(paste("scale needs to either a single rate or vector of length ", length(snps)-1 ," giving the rates between each snp", sep=""))
	}
	if(sum(snps-sort(snps))!=0){
		warning("Input snps were not sorted")
		snps <- sort(snps)
	}
	haldane <- function(d){ (1/2)*(1- exp(-2*d))}
	displacement <- sapply(2:length(snps), function(x){
		return(snps[x]-snps[x-1])
		})
	if(!all(displacement>0)){
		stop("The physical distance between all snps must be >0")
	}

	# 1==recombination occurs, 0==no recombination
	out <- sapply(1:length(displacement), function(i) rbinom(n=1, size=1, prob=haldane(displacement[i]*c[i])))
	return(out)
}


#' @title Simulate recombination (crossovers or non-crossovers) along a chromosome
#' 
#' @description This function simulates recombination between two parent chromosomes generated from class
#' \code{parent.genomes}.
#' 
#' @param parents an object of class \code{parent.genomes} specifying the parental ancestry 
#' at each simulated snp.
#' 
#' @param r.index a vector of length \code{snps-1} specifying whether a recombination is to be 
#' simulated (1) or not (0) in between two adjacent snps. 
#' 
#' @param mu.rate a numeric between 0 and 1 (inclusive) specifying the per snp mutation rate. 
#' 
#' @param f.cross a numeric between 0 and 1 (inclusive) giving the frequency of recombination 
#' events that result in crossing over. This is same as 1 minus the frequency of non-crossovers.
#' 
#' @param f.convert a numeric between 0 and 1 (inclusive) that gives the frequency of gene conversion 
#' during recombination.
#' 
#' @param length.conversion an integer specifying the mean (and variance) of a given gene conversion 
#' tract (in bps).
#' 
#' @details This function simulates crossover events (with or without gene conversions)
#' and non-crossover events given the parental snps and location of recombination events.
#' First, parental chromosomes are duplicated into 2 pairs of sister chromatids. During this step
#' mutations occur at rate \code{mu.rate} that switch parental alleles. Second, recombination is
#' simulated by systematically going down the length of the chromosome. At recombination points, idenitfied in \code{r.index},
#' two non-sister chromatids are picked at random. Third, a CO or NCO is simulated; this is done
#' at a frequency of f.cross or 1-f.cross, respectively. If a CO occurs, the genotypes of the 3'
#' end of the non-sister chromatids are simply switched. If a NCO occurs, no such switch takes place. 
#' Note that is is impossible to detect NCO events without an associated gene conversion tract. 
#' Following every CO or NCO event gene conversion is simulated at a rate of f.convert. During a 
#' gene conversion, one of the two unpicked chromatids' genotypes are switched. The length of the 
#' gene conversion tract is sampled from a poission distribution with mean (variance) of length.conversion.
#' Note that the gene conversion tract will stop at the end of the chromosme if necessary. 
#' 
#' @return an object of class \code{recombine} that contains the following data:
#' \itemize{
#' 	\item{parents}{ input data}
#' 	\item{r.index}{ input data}
#' 	\item{m.rate}{ input data}
#' 	\item{f.cross}{ input data}
#' 	\item{f.convert}{ input data}
#' 	\item{length.conversion}{ input data}
#' 	\item{chromatids.mutated_snps}{ a list giving the genotypes of 4 chromatids (two pairs of sister chromatids) following mutation}
#' 	\item{chromatids.recombined}{ a list giving the genotypes of 4 chromatids following recombination}
#' }
#' 
#' 
#' @seealso \code{\link{recombine_index}}
#' 
#' @author Tyler D. Hether
#' 
#' @export recombine
#' 
#' @examples
#' set.seed(1234567)        # For reproducibility
#' l <- 50                  # number of snps to simulate
#' c <- 3e-05               # recombination rate between snps (Morgan/bp)
#' snps <- c(1:l)*1e4       # snps are evenly spaced 10kbp apart
#' p <- make_parents(snps)  # make the parents
#' #
#' # The recombination indeces are:
#' r_points <- recombine_index(scale=c, snps=snps) 
#' #
#' # Now recombine the parents:
#' recomb_sim <- recombine(parents=p, r.index=r_points, mu.rate=1e-05, 
#' 	f.cross=0.6, f.convert=0.2, length.conversion=15000) 
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
	# type = 1
	for(i in 1:(l-1)){
		# Cases where there is a recombination event, pick
		# two of the four chromatids to recombine:

		if(i %in% which(r.index==1)){
			# If recobomination with sister chromatids is allowed:
			# picked.chromatids <- sample(c(1:4), 2, replace=FALSE)

			# If recobomination only with non-sister chromatids is allowed:
			vals <- as.numeric(as.data.frame(chromatids.recombined)[i,])

			# In some rare cases of frequent gene conversion, vals will not have two 0s and two 1s. 
			# If that's the case, sample 1 or 2 and 3 or 4. 
			if(sum(vals)==2){
				picked.chromatids <- c(sample(which(vals==0), 1), sample(which(vals==1), 1))
			} else {
				picked.chromatids <- c(sample(c(1,2), 1), sample(c(3,4), 1))
			}

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

					# find which snps are to the 3' of i, x
					x = which(p$snp>p$snp[i])
					# and the snps that are 5' of (length_of_conversion+i), y
					y = which(p$snps<(p$snps[i]+length_of_conversion))
					# If there are snps in both cats then they will be converted:
					if(length(intersect(x,y))!=0){
						chromatids.recombined[[to_convert]][c(intersect(x,y))] <- sapply(chromatids.recombined[[to_convert]][c(intersect(x,y))], function(j){switch_values(j)})
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
#' @description This function populates simulated data of class \code{recombine} with
#' read counts given a specified coverage and assignment probability.
#' 
#' @param simdata an object of class \code{recombine}
#' 
#' @param p.assign a numeric between 0 and 1 (inclusive) that gives the probability of 
#' correct sequencing assignment. A value of 1 means that no sequencing error or ancestral 
#' polymorphism is to be simulated. 
#' 
#' @param coverage the mean (and variance) of sequencing read coverage to be simlated.  \code{coverage} is 
#' sampled from a poisson distribution.
#'
#' @return an object of class \code{simulated.coverage} providing the a list of simulated sequencing reads for 
#' each parental type for each of the four haploid recombinants. There are four lists in the output -- one for
#' each recombinant. Within each list the read counts for parent 0 and parent 1 are given as well as the snp locations
#' and given (simulated) states. 
#' 
#' @seealso \code{\link{recombine}}
#' 
#' @author Tyler D. Hether 
#' 
#' @export simulate_coverage
#' 
#' @examples
#' set.seed(1234567)        # For reproducibility
#' l <- 50                  # number of snps to simulate
#' c <- 3e-05               # recombination rate between snps (Morgan/bp)
#' snps <- c(1:l)*1e4       # snps are evenly spaced 10kbp apart
#' p_a <- 0.97              # assignment probability
#' coverage <- 1            # mean coverage
#' p <- make_parents(snps)  # make the parents
#' #
#' # The recombination indeces are:
#' r_points <- recombine_index(scale=c, snps=snps) 
#' #
#' # Now recombine the parents:
#' recomb_sim <- recombine(parents=p, r.index=r_points, mu.rate=1e-05, 
#' 	f.cross=0.6, f.convert=0.2, length.conversion=15000) 
#' # And finally simulated read coverage for each haploid recombinant
#' sim_cov <- simulate_coverage(simdata=recomb_sim, p.assign=p_a, coverage=coverage) # simulate sequencing coverage
#' sim_cov

simulate_coverage <- function(simdata, p.assign, coverage){
	if(!inherits(simdata, "recombine")){
		stop(paste("Object", simdata, " needs to be of class 'recombine'.", sep=""))
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

	class(out) <- c("list", "simulated.coverage")
	return(out)
}

#' @title Simulate recombination across a given number of meiosis events
#' 
#' @description This is a wrapper function of many other functions that simulates a given number of
#' tetrads each with four haploid spores that are recombinants between two parents.  
#' 
#' @param n.tetrads An integer specifying the number of tetrads to simulate.
#' 
#' @param \code{scale} either vector of length 1 specifying the genome-wide recombination rate (Morgans/bp)
#' or a vector of length \code{snps-1} specifying the recombination rate between all snps. In either
#' case rates must be positive.
#'
#' @param \code{snps} a vector gving the locations of snps along a chromosome
#'
#' @param p.assign a numeric between 0 and 1 (inclusive) that gives the probability of 
#' correct sequencing assignment. A value of 1 means that no sequencing error or ancestral 
#' polymorphism is to be simulated. 
#' 
#' @param mu.rate a numeric between 0 and 1 (inclusive) specifying the per snp mutation rate. 
#' 
#' @param r.index a vector of length \code{snps-1} specifying whether a recombination is to be 
#' simulated (1) or not (0) in between two adjacent snps. 
#' 
#' @param f.cross a numeric between 0 and 1 (inclusive) giving the frequency of recombination 
#' events that result in crossing over. This is same as 1 minus the frequency of non-crossovers.
#' 
#' @param f.convert a numeric between 0 and 1 (inclusive) that gives the frequency of gene conversion 
#' during recombination.
#' 
#' @param length.conversion an integer specifying the mean (and variance) of a given gene conversion 
#' tract (in bps).
#'
#' @param chr.name a numeric or character value specifying the chromosome name (default to "I") 
#' 
#' @return A data.frame of class \code{tetrad}. Each row contains the following information:
#' \itemize{
#' 	\item{Tetrad}{ The tetrad ID}
#' 	\item{Spore}{ The spore ID (1:4)}
#' 	\item{Chr}{ The chromosome name}
#'  \item{Snp}{ The snp location (in bps)}
#'  \item{p0}{ The number of reads that mapped to parent 0}
#'  \item{p1}{ The number of reads that mapped to parent 1}
#'  \item{states_given}{ The simulated states (0 or 1)}
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
#' set.seed(1234567)        # For reproducibility
#' n_tetrads <- 100         # number of tetrads (meiosis events)
#' l <- 50                  # number of snps to simulate
#' c <- 3e-05               # recombination rate between snps (Morgan/bp)
#' snps <- c(1:l)*2e4       # snps are evenly spaced 20kbp apart
#' p_a <- 0.97              # assignment probability
#' coverage <- 1            # mean coverage
#' # Now simulate
#' sim100 <- sim_tetrad(n.tetrads=n_tetrads, scale=c, snps=snps, 
#' 	p.assign=p_a, mu.rate=0, f.cross=0.8, f.convert=0.3, 
#' 	length.conversion=2e3, coverage=2)
#' sim100

sim_tetrad <- function(n.tetrads, scale, snps, p.assign, mu.rate, f.cross, f.convert, 
        length.conversion, coverage, chr.name="I"){
	p <- make_parents(snps)
    
    res <- lapply(1:n.tetrads, function(Z, ...){

        r <- recombine_index(scale, snps)
        # print(sum(r)/snps[l])
        
        recomb_sim <- recombine(parents=p, r.index=r, mu.rate=mu.rate, f.cross=f.cross, 
                f.convert=f.convert, length.conversion=length.conversion)
        sim_reads <- simulate_coverage(simdata=recomb_sim, p.assign=p.assign, coverage=coverage)
        class(sim_reads) <- list("individual.tetrad")
        return(sim_reads)
    })
    
    # Convert to proper dataframe:
    dat0 <- lapply(1:length(res), function(x){

      one <- res[[x]][[1]]
      two <- res[[x]][[2]]
      three <- res[[x]][[3]]
      four <- res[[x]][[4]]

      return(rbind(
        data.frame(Tetrad=x, Spore=1, Chr=chr.name, Snp=one$snps, p0=one$p0.assign, p1=one$p1.assign, states_given=one$states_given),
        data.frame(Tetrad=x, Spore=2, Chr=chr.name, Snp=two$snps, p0=two$p0.assign, p1=two$p1.assign, states_given=two$states_given),
        data.frame(Tetrad=x, Spore=3, Chr=chr.name, Snp=three$snps, p0=three$p0.assign, p1=three$p1.assign, states_given=three$states_given),
        data.frame(Tetrad=x, Spore=4, Chr=chr.name, Snp=four$snps, p0=four$p0.assign, p1=four$p1.assign, states_given=four$states_given)
      ))

      })
    out <- do.call(rbind, dat0)

    class(out) <- c("data.frame", "tetrad")
    return(out)
}







# tmpres <- lapply(1:length(out), function(x){
# 	tmp <- out[[x]][[1]]$states_given

# 	return(sum(abs(tmp[2:l] - tmp[1:(l-1)]))/(range(snps)[2]-range(snps)[1]))
# 	})








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
#' set.seed(1234567) # For reproducibility
#' # simulate a recombination hotspot between the 100th and 101st snp
#' l=100; scale = 0.01; snps=1:l
#' n.spores <- 500 # number of spores to simulate
#' spores <- sim_en_masse(n.spores=n.spores, scale=scale, snps=snps, 
#' p.assign=.999, mu.rate=0.001, f.cross=0.5, 
#'     f.convert=0.5, length.conversion=10, coverage=1)

sim_en_masse <- function(n.spores, scale, snps, p.assign, mu.rate, f.cross, f.convert, length.conversion, coverage){

    out <- lapply(1:n.spores, function(Z, ...){

        r <- recombine_index(scale, snps)
        p <- make_parents(snps)
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


