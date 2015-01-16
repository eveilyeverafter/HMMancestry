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

	out <- list(snps=1:L, p1=p1, p2=p2)
	class(out) <- c("list", "parent.genomes")
	return(out)

}

# Function to make sure values in a vector, a, are between x and y, inclusive:
check_values <- function(a,x,y){
	if(max(a)>y){
		stop(paste("The maximum of 'a' needs to be <", y, sep=""))
	}
	if(min(a)<x){
		stop(paste("The minimum of 'a' needs to be >", x,sep=""))
	}
}

# Simulate recombination events given the vector transistion probs and the number of loci
recombine_index <- function(p_trans){
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

mutate_snps <- function(xx, mu=mu){
	# MLG <- multilocus genotype
	ll <- length(xx)
	to.mutate_snps <- sample(x=c(0,1), ll, prob=c(1-mu,mu), replace=TRUE)
	if(sum(to.mutate_snps>0)){
		xx[which(to.mutate_snps==1)] <- sapply(xx[which(to.mutate_snps==1)], switch_values)
	}
	return(xx)
}


# f.cross <- 0.5 # frequency of crossovers during recombination
# f.convert <- 1 # frequency of gene conversion during recombination
# length.conversion <- 20 # average length of gene conversion track
# Note: NCO can occur two ways: but crossing over with sister chromotids (which is only noticed when there are
	# mutations present to differentiate sister chromotids) and the alternate ways in which Holliday Junctions resolve.
# Note: mutations during simulated meiosis create situations that look like 1 basepair non-crossover events.

# recombine_proposed(parents=p, r.index=r, mu.rate=0, f.cross=0.5, f.convert=1, length.conversion=20)
recombine <- function(parents, r.index, mu.rate=0, f.cross=0.5, f.convert=0, length.conversion=20){
	if(!inherits(parents, "parent.genomes")){
		stop(paste("Object ", parents, " needs to be of class parent.genomes", sep=""))
	}
	l <- length(parents$p1)

	# S-phase and mutation
	chromotids.orig <- list(p1_1=parents$p1, p1_2=parents$p1, p2_1=parents$p2, p2_2=parents$p2)
	chromotids.mutated_snps <- lapply(chromotids.orig, function(i, ...){
		mutate_snps(i, mu=mu.rate)
		})
	
	# Recombine chromotids:
	chromotids.recombined <- chromotids.mutated_snps
	for(i in 1:(l-1)){
		# Cases where there is a recombination event, pick
		# two of the four chromotids to recombine:
		if(i %in% which(r.index==1)){
			picked.chromotids <- sample(c(1:4), 2, replace=FALSE)
			# Recombine them:
			chromotids <- data.frame(one=chromotids.recombined[[picked.chromotids[1]]], 
					   two=chromotids.recombined[[picked.chromotids[2]]])
			one <- chromotids[(i+1):l,1]
			two <- chromotids[(i+1):l,2]
				# Does recombination result in a crossover (1) or non-crossover (0)?
				to_cross <- cross_vs_noncross(f.cross=f.cross)
				if(to_cross==1){
					chromotids[(i+1):l,] <- cbind(two, one)
				}			
				# Is there gene conversion at this recombination point?
				if(sample(c(1,0),1,prob=c(f.convert,1-f.convert))==1){
					# simulate gene conversion on one of the non-picked chromotids:
					to_convert <- sample(c(1:4)[-picked.chromotids],1) # non-picked chromotid
					length_of_conversion <- rpois(n=1,lambda=length.conversion)
					# if the length of the track extends pass the chromosome then stop at the end of the chromosome.
					# Otherwise, create a gene conversion track of length length_of_conversion starting at the 
					# recombination point. 
					if( (i+1+length_of_conversion) > l){
						chromotids.recombined[[to_convert]][c((i+1):l)] <- sapply(chromotids.recombined[[to_convert]][c((i+1):l)], function(j){switch_values(j)})
					} else {
						chromotids.recombined[[to_convert]][c((i+1):(i+1+length_of_conversion))] <- sapply(chromotids.recombined[[to_convert]][c((i+1):(i+1+length_of_conversion))], function(j){switch_values(j)})
					}

				}
			# Save results:
			chromotids.recombined[[picked.chromotids[1]]] <- chromotids[,1]
			chromotids.recombined[[picked.chromotids[2]]] <- chromotids[,2]
		}
	}

	out <- list(parents=parents, r.index=r.index, mu.rate=mu.rate, 
			chromotids.mutated_snps=chromotids.mutated_snps, picked.chromotids=picked.chromotids, 
			chromotids.recombined=chromotids.recombined)
	class(out) <- c("list", "recombine")
	return(out)	
}


# Given the frequency of crossover events, f.cross, sample whether one is to occur (1=yes, crossover; 0=no, non-crossover).
cross_vs_noncross <- function(f.cross){
	sample(c(1,0), 1, prob=c(f.cross, 1-f.cross))
}


recombine_dep <- function(parents, r.index, mu.rate){
	if(!inherits(parents, "parent.genomes")){
		stop(paste("Object ", parents, " needs to be of class parent.genomes", sep=""))
	}
	l <- length(parents$p1)

	# S-phase and mutation
	chromotids.orig <- list(p1_1=parents$p1, p1_2=parents$p1, p2_1=parents$p2, p2_2=parents$p2)
	chromotids.mutated_snps <- lapply(chromotids.orig, function(i, ...){
		mutate_snps(i, mu=mu.rate)
		})
	
	# Recombine chromotids:
	chromotids.recombined <- chromotids.mutated_snps
	for(i in 1:(l-1)){
		# Cases where there is a recombination event, pick
		# two of the four chromotids to recombine:
		if(i %in% which(r.index==1)){
			picked.chromotids <- sample(c(1:4), 2, replace=FALSE)
			# Recombine them:
			chromotids <- data.frame(one=chromotids.recombined[[picked.chromotids[1]]], 
					   two=chromotids.recombined[[picked.chromotids[2]]])
			one <- chromotids[(i+1):l,1]
			two <- chromotids[(i+1):l,2]
			chromotids[(i+1):l,] <- cbind(two, one)	
			# Save results:
			chromotids.recombined[[picked.chromotids[1]]] <- chromotids[,1]
			chromotids.recombined[[picked.chromotids[2]]] <- chromotids[,2]
		}
	}

	out <- list(parents=parents, r.index=r.index, mu.rate=mu.rate, 
			chromotids.mutated_snps=chromotids.mutated_snps, picked.chromotids=picked.chromotids, 
			chromotids.recombined=chromotids.recombined)
	class(out) <- c("list", "recombine")
	return(out)	
}

# Simulate sequencing of x coverage across simulated data: 
simulate_coverage <- function(a, p_assign, coverage){
	if(!inherits(a, "recombine")){
		stop(paste("Object 'a' needs to be of class 'recombine'.", sep=""))
	}

	simulated.reads.total <- lapply(a$chromotids.recombined, function(i) return(rpois(i, coverage)))
	
	# This block returns the number of simulated reads for each parental type, 0 and 1, for each snp position.
	out <- lapply(1:4, function(i){
		correct_reads <- rbinom(size=simulated.reads.total[[i]], prob=p_assign, n=simulated.reads.total[[i]])
		incorrect_reads <- simulated.reads.total[[i]] - correct_reads 

		p0.assign <- sapply(1:length(a$chromotids.recombined[[i]]), function(x){
				if(a$chromotids.recombined[[i]][x]==0){
					return(correct_reads[x])
				}
				if(a$chromotids.recombined[[i]][x]==1){
					return(incorrect_reads[x])
				}
			})
		p1.assign <- sapply(1:length(a$chromotids.recombined[[i]]), function(x){
				if(a$chromotids.recombined[[i]][x]==1){
					return(correct_reads[x])
				}
				if(a$chromotids.recombined[[i]][x]==0){
					return(incorrect_reads[x])
				}
			})
		return(list(p0.assign=p0.assign, p1.assign=p1.assign))

	})
	out[['snps']] <- a$parents$snps

	class(out) <- c("list", "snp.recom")
	return(out)
}

















