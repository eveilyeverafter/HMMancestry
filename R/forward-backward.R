#' @title Inferring genotypic states using the Forward-Backward algorithm when all spores of a 
#' tetrad are genotyped
#' 
#' @description Uses the forward-backward algorithm to estimate ancestral genotypes 
#' along a given chromosome for a given genotyped tetrad or simulated data.
#' 
#' @param snp.dat an object of class \code{snp.recom} that gives the bi-allelic read counts
#' (e.g., from bowtie2 or from \code{\link{simulate_coverage}}) of each parental type 
#' for each spore along a given chromosome
#' 
#' @param tetrad.id (optional) a character or numeric giving the tetrad identity. 
#' If snp.dat is of class \code{forward.backward} then \code{tetrad.id} will be inherited 
#' from that class.
#' 
#' @param chr.name (optional) a character (typical) or numeric giving the chromosome identity. 
#' If snp.dat is of class \code{forward.backward} then \code{chr.name} will be inherited 
#' from that class.
#' 
#' @param p.assign a numeric between 0 and 1 (exclusive) that specifies the probability of 
#' correct sequencing assignment (i.e., 1 - sequencing error rate) plus the probability of 
#' new a new mutation at a given snp (see details)
#' 
#' @param p.trans the per snp (check that) recombination (i.e., transition) probability. 
#' 
#' @details \code{tetrad_estimate_anc_fwd_back} attempts to estimate parental genotypic 
#' 'states' along a chromosome given empirical or simulated F2 cross data. Next-generation 
#' data inherits both sequencing error and missing data -- especially when sequencing 
#' coverage is low. In these two cases (sequence error and missing data) the parental state 
#' is ambiguous or unknown, respectively. \code{tetrad_estimate_anc_fwd_back} takes these
#' uncertainities into account. 
#' 
#' To do - Add details of the actual algorithm
#' 
#' @return to do
#' 
#' @references Hohenlohe, P.A., S. Bassham, M. Currey, and W.A. Cresko. 2012. Extensive linkage 
#' disequilibrium and parallel adaptive divergence across threespine stickleback genomes. 
#' Phil. Trans. R. Soc. B, 395-408, doi: 10.1098/rstb.2011.0245
#' 
#' @seealso \code{\link{simulate_coverage}}
#' 
#' @author Tyler D. Hether
#' 
#' @export tetrad_est_fwd_back
#' 
#' @examples
#' set.seed(1234567) # For reproducability
#' l <- 10 # number of loci to simulate
#' rec <- 0.01 # recombination rate between each snp
#' r <- recombine_index(rep(rec, l-1)) # recombination rate between each snp (vector form)
#' # probability of correct sequencing assignment (1-sequence error rate) 
#' # + probabilty of mutation that makes snps identical by state
#' p_a <- .999 
#' p <- make_parents(l) # make the parent
#' recomb_sim <- recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=1, length.conversion=10) # recombine parents
#' sim_reads <- simulate_coverage(simdata=recomb_sim, p.assign=p_a, coverage=1) # simulate sequencing coverage
#' # Use the forward-backward algorithm to get the posterior probability of parent '0' ancestry and infer states
#' fbres <- tetrad_est_fwd_back(snp.dat=sim_reads, tetrad.id=2, chr.name="II", p.assign=p_a, p.trans=rec)
#' # Need S3 methods for printing and plots results

tetrad_est_fwd_back <- function(snp.dat, tetrad.id=1, chr.name="I", p.assign, p.trans){

    # Run this entire block of code for each of the four spores:
    all_spores <- lapply(1:4, function(x, ...){
        spore_number <- x
  
        snp.locations <- snp.dat$snps

        check_est_input(snp.dat, spore_number, chr.name, snp.locations, p.assign, p.trans)
        # 
        displace <- snp.locations[-c(1)]-snp.locations[1:(length(snp.locations)-1)]
          
        n_snps <- length(snp.locations)

        # 1) establish matrices/vectors:
        emissions <- forward <- backward <- posterior <- matrix(0, ncol=2, nrow=n_snps)
        scale <- numeric(n_snps)
        #
        # 2) calculate emission probabilities for all states at all positions
        # data at position i are counts n_ij of allele matching parent j, with n_i = sum(n_ij)
        # emissions[i][j] = (n_i choose n_ij) * (1 - eps)^(n_ij) * (eps)^(n_i - n_ij)
        # if no data at SNP i, emissions[i][j] = 1 for all j
        #
        k0 <- snp.dat[[spore_number]]$p0.assign
        k1 <- snp.dat[[spore_number]]$p1.assign
        n_i <- k0+k1
        p0 <- choose(n_i,k0)*(p.assign^k0)*((1-p.assign)^(n_i-k0))
        p1 <- choose(n_i,k1)*(p.assign^k1)*((1-p.assign)^(n_i-k1))    
        emissions <- t(rbind(p0, p1))
        #
        # 3) calculate pi[k] (vector of state probabilities at first position -- here use 1/k for each)
        #
        pi_initial <- c(0.5, 0.5)
        # 4) calculate forward probabilities and scaling factor:
        forward[1,] <- pi_initial*emissions[1,]
        scale[1] <- sum(forward[1,])
        forward[1,] <- forward[1,]/scale[1]         #re-scale forward probabilities by their sum to avoid underflow
        for(i in 2:n_snps){
            transition_matrix_i <- matrix(c(1-(displace[i-1]*p.trans), displace[i-1]*p.trans, displace[i-1]*p.trans, 1-(displace[i-1]*p.trans)), byrow=TRUE, ncol=2)
            forward[i,1] <- sum(transition_matrix_i[1,1]*forward[i-1,1], transition_matrix_i[2,1]*forward[i-1,2]) * emissions[i,1]
            forward[i,2] <- sum(transition_matrix_i[1,2]*forward[i-1,1], transition_matrix_i[2,2]*forward[i-1,2]) * emissions[i,2]
            # where transition_matrix_i is scaled by the distance between position i-1 and i
            scale[i] = sum(forward[i,])
            forward[i,] = forward[i,]/scale[i]
        }

        # 5) calculate backward probabilities:
        backward[n_snps,] = 1/scale[n_snps]
        for(i in (n_snps-1):1){
            transition_matrix_i <- matrix(c(1-(displace[i]*p.trans), displace[i]*p.trans, displace[i]*p.trans, 1-(displace[i]*p.trans)), byrow=TRUE, ncol=2)
            backward[i,1] <- sum(transition_matrix_i[1,1]*backward[i+1,1], transition_matrix_i[2,1]*backward[i+1,2]) * emissions[i+1,1]
            backward[i,2] <- sum(transition_matrix_i[1,2]*backward[i+1,1], transition_matrix_i[2,2]*backward[i+1,2]) * emissions[i+1,2]
            # transition matrix scaled by distance from i to i+1
            # note the i+1 for the emissions, not i
            backward[i,] = backward[i,]/scale[i]
        }
                
        # 6) caculate posteriors: 
        for(i in 1:n_snps){
            posterior[i,] <- (forward[i,]*backward[i,])/(forward[i,1]*backward[i,1] + forward[i,2]*backward[i,2])
        }
        # Bug: in rare cases the denominator of the above formula returns 0. This
        # happens when there is the foward and backward probabilities each return 0 for alternate states.
        # The result is "NaN" produced as a result. 
        # The function check_tie (infer-recom.R) will treat these cases as missing data. 


        # 7) total likelihood:
        lnL = sum(log(scale))

        # 8) Call states 0 or 1 based on the posterior probs. If it's a tie, pick one at random:
        states_inferred <- apply(posterior, 1, function(i){
            if(check_tie(i)==0){
                return(which(i==max(i))-1)
            }
            if(check_tie(i)==1 | NaN %in% i){
                # The first is typically true for low coverage data and address bug #6.
                # The second typically occurs with high coverage data and addresses bug #4.
                return(sample(c(0,1), 1, prob=c(0.5,0.5)))
            }
        })

        out <- list(snp.dat=snp.dat, spore_number=spore_number, tetrad.id=tetrad.id, chr.name=chr.name, 
                 snp.locations=snp.locations, p.assign=p.assign, p.trans=p.trans,
                 emissions=emissions, forward=forward, backward=backward, scale=scale, 
                 posterior=posterior, states_inferred=states_inferred, lnL=lnL)

        # class(out) <- c("list", "forward.backward")
        return(out)
    })
    class(all_spores) <- c("list", "forward.backward")
    return(all_spores)
}


#' @title Infer genotypic states using the Forward-Backward algorithm when spores are randomly
#' sampled en masse
#'
#' @description to do 
#' 
#' @param to do 
#' 
#' @return to do
#' 
#' @references Hohenlohe, P.A., S. Bassham, M. Currey, and W.A. Cresko. 2012. Extensive linkage 
#' disequilibrium and parallel adaptive divergence across threespine stickleback genomes. 
#' Phil. Trans. R. Soc. B, 395-408, doi: 10.1098/rstb.2011.0245
#' 
#' @seealso \code{\link{tetrad_est_fwd_back}}
#' 
#' @author Tyler D. Hether
#' 
#' @export est_fwd_back
#' 
#' @examples
#' set.seed(1234567) # For reproducability
#' # simulate a recombination hotspot between the 100th and 101st snp
#' rec <- c(rep(0.001, 99), 0.1, rep(0.001, 99))
#' n.spores <- 500 # number of spores to simulate
#' spores <- sim_en_masse(n.spores=n.spores, l=200, rec=rec, 
#'  p.assign=.999, mu.rate=0.001, f.cross=0.5, 
#'     f.convert=0.5, length.conversion=10, coverage=1)
#' 
#' # Run the fb algorithm to estimate the parental states:
#' Allspores <- lapply(1:n.spores, function(Z){
#'         fbres2 <- est_fwd_back(single.snp.dat=spores[[Z]], 
#'  spore_number=Z, chr.name="I", p.assign=0.999, p.trans=mean(rec))
#'         return(fbres2)
#'     })
#' 
#' df <- do.call(rbind,lapply(Allspores, function(i){
#'  return(as.numeric(i$states_inferred))}))
#' 
#' plot(apply(t(apply(df, 1, id_hotspots)), 2, sum), 
#'  type="l", xlab="snp", ylab="Number of recombination events")

est_fwd_back <- function(single.snp.dat, spore_number=1, chr.name="I", p.assign, p.trans){
 
    snp.locations <- single.snp.dat$snps

    check_est_input_single(single.snp.dat, chr.name, snp.locations, p.assign, p.trans)
    # 
    displace <- snp.locations[-c(1)]-snp.locations[1:(length(snp.locations)-1)]
          
    n_snps <- length(snp.locations)

    # 1) establish matrices/vectors:
    emissions <- forward <- backward <- posterior <- matrix(0, ncol=2, nrow=n_snps)
    scale <- numeric(n_snps)
    #
    # 2) calculate emission probabilities for all states at all positions
    # data at position i are counts n_ij of allele matching parent j, with n_i = sum(n_ij)
    # emissions[i][j] = (n_i choose n_ij) * (1 - eps)^(n_ij) * (eps)^(n_i - n_ij)
    # if no data at SNP i, emissions[i][j] = 1 for all j
    #
    k0 <- single.snp.dat$p0.assign
    k1 <- single.snp.dat$p1.assign
    n_i <- k0+k1
    p0 <- choose(n_i,k0)*(p.assign^k0)*((1-p.assign)^(n_i-k0))
    p1 <- choose(n_i,k1)*(p.assign^k1)*((1-p.assign)^(n_i-k1))    
    emissions <- t(rbind(p0, p1))
    #
    # 3) calculate pi[k] (vector of state probabilities at first position -- here use 1/k for each)
    #
    pi_initial <- c(0.5, 0.5)
    # 4) calculate forward probabilities and scaling factor:
    forward[1,] <- pi_initial*emissions[1,]
    scale[1] <- sum(forward[1,])
    forward[1,] <- forward[1,]/scale[1]         #re-scale forward probabilities by their sum to avoid underflow
    for(i in 2:n_snps){
        transition_matrix_i <- matrix(c(1-(displace[i-1]*p.trans), displace[i-1]*p.trans, displace[i-1]*p.trans, 1-(displace[i-1]*p.trans)), byrow=TRUE, ncol=2)
        forward[i,1] <- sum(transition_matrix_i[1,1]*forward[i-1,1], transition_matrix_i[2,1]*forward[i-1,2]) * emissions[i,1]
        forward[i,2] <- sum(transition_matrix_i[1,2]*forward[i-1,1], transition_matrix_i[2,2]*forward[i-1,2]) * emissions[i,2]
        # where transition_matrix_i is scaled by the distance between position i-1 and i
        scale[i] = sum(forward[i,])
        forward[i,] = forward[i,]/scale[i]
    }

    # 5) calculate backward probabilities:
    backward[n_snps,] = 1/scale[n_snps]
    for(i in (n_snps-1):1){
        transition_matrix_i <- matrix(c(1-(displace[i]*p.trans), displace[i]*p.trans, displace[i]*p.trans, 1-(displace[i]*p.trans)), byrow=TRUE, ncol=2)
        backward[i,1] <- sum(transition_matrix_i[1,1]*backward[i+1,1], transition_matrix_i[2,1]*backward[i+1,2]) * emissions[i+1,1]
        backward[i,2] <- sum(transition_matrix_i[1,2]*backward[i+1,1], transition_matrix_i[2,2]*backward[i+1,2]) * emissions[i+1,2]
        # transition matrix scaled by distance from i to i+1
        # note the i+1 for the emissions, not i
        backward[i,] = backward[i,]/scale[i]
    }
                
    # 6) caculate posteriors: 
    for(i in 1:n_snps){
        posterior[i,] <- (forward[i,]*backward[i,])/(forward[i,1]*backward[i,1] + forward[i,2]*backward[i,2])
    }
    # Bug: in rare cases the denominator of the above formula returns 0. This
    # happens when there is the foward and backward probabilities each return 0 for alternate states.
    # The result is "NaN" produced as a result. 
    # The function check_tie (infer-recom.R) will treat these cases as missing data. 


    # 7) total likelihood:
    lnL = sum(log(scale))

    # 8) Call states 0 or 1 based on the posterior probs. If it's a tie, pick one at random:
    states_inferred <- apply(posterior, 1, function(i){
        if(check_tie(i)==0){
            return(which(i==max(i))-1)
        }
        if(check_tie(i)==1 | NaN %in% i){
            # The first is typically true for low coverage data and address bug #6.
            # The second typically occurs with high coverage data and addresses bug #4.
            return(sample(c(0,1), 1, prob=c(0.5,0.5)))
        }
    })

    out <- list(single.snp.dat=single.snp.dat, spore_number=spore_number, chr.name=chr.name, 
             snp.locations=snp.locations, p.assign=p.assign, p.trans=p.trans,
             emissions=emissions, forward=forward, backward=backward, scale=scale, 
             posterior=posterior, states_inferred=states_inferred, lnL=lnL)

    class(out) <- c("list", "single.forward.backward")
    return(out)
}

# Minor functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
check_est_input <- function(snp.dat, spore_number, chr.name, snp.locations, p.assign, p.trans, ...){
    if (!inherits(snp.dat, "snp.recom"))
        stop("Object must be of class 'snp.recom'. See 'make_snp.data()'")
    check_values(p.assign, 0, 1)
    check_values(p.trans, 0, 1)
    if(length(snp.dat[[1]]$p0.assign)!=length(unique(snp.locations))){
        stop("Each snp needs 1 unique numeric id.")
    }
    if(!spore_number %in% c(1:4)){
        stop("spore_number must be an the integer 1, 2, 3, or 4")
    }
}

check_est_input_single <- function(single.snp.dat, chr.name, snp.locations, p.assign, p.trans){
    if (!inherits(single.snp.dat, "single.spore"))
        stop("Object must be of class 'single.spore'")
    check_values(p.assign, 0, 1)
    check_values(p.trans, 0, 1)
    if(length(single.snp.dat$p0.assign)!=length(unique(snp.locations))){
        stop("Each snp needs 1 unique numeric id.")
    }
}




#' @method print estimate_anc_fwd_back
#' @export 
print.estimate_anc_fwd_back <- function(x, ...){
    ## make pretty output format
    ## print to screen
}

