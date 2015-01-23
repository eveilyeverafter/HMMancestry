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
#' @details \code{tetrad_estimate_anc_fwd_back} attempts to estimate 
#' parental genotypic 'states' along a chromosome given empirical or 
#' simulated F2 cross data were all four spores of a tetrad are genotyped. 
#' Next-generation data inherits both sequencing error and missing data 
#' -- especially when sequencing coverage is low. In these two cases 
#' (sequence error and missing data) the parental state is ambiguous 
#' or unknown, respectively. \code{tetrad_estimate_anc_fwd_back} takes 
#' these uncertainties into account.
#' 
#' The three steps taken in \code{tetrad_estimate_anc_fwd_back} are: 
#' \enumerate{
#'  \item Calculate forward probabilities for each state (5' to 3')
#'  \item Calculate backward probabilities for each state (3' to 5')
#'  \item From these two probabilities, calculate the posterior probability
#'  that a snp location is of a given parental state.
#' } 
#' 
#' The two element vector of forward probabilities for each parental state 
#' \equ{f_{i}} at each position, i, are calculated as:
#' \dequ{f_{i} = e_{i}T_{i}f_{i-1}}
#' 
#' where \equ{e_{i}} is the emission probabilities for each state (see below), \equ{T_{i}}
#' is a 2-by-2 matrix describing the transition (recombination) probabilities
#' between two states, and \equ{f_{i-1}} is the forward probability at the 
#' previous position along the chromosome (5' of position i). It is assumed 
#' that each state is equally likely to occur at the first position. 
#' 
#' 
#' The emmission probabilities are calculated independently for 
#' each snp and depends on the sequence reads assigned to parent "0" 
#' and parent "1" and \code{p.assign}. The emission probabilities are
#' calculated using the binomial equation. For example, the emission
#' probability for parental state "0" at snp position i is: 
#' 
#' \dequ{e_{i,0} = {n \choose{k0}} p^{k0} (1-p)^{n-k0}}
#' 
#' where n is the sum of the number of reads from both parents and p 
#' is \code{p.assign}. 
#' 
#' Like the emission probabilities, the transition probabilities are also
#' calculated for each snp position. This accounts for the displacement 
#' between snps.  For example if the per base recombination rate is 
#' 0.001 and two snps are 10 base pairs apart, the transition probability of
#' is 0.01 (i.e., the probability of not recombining would be 0.99).
#' 
#' To avoid underflow, we rescale the forward (and backward) probabilities
#' each iteration to that they sum to unity. 
#' 
#' The backward probabilities are calculated similarly but in the 3' to 5' 
#' direction. It is assumed that the backward probability is equally likely 
#' in each state. Following Durbin et al (1998) the backward probability at
#' snp position i is: 
#' 
#' \dequ{b_{i} = T_{i}e_{i+1}b_{i+1}}
#' 
#' Again, these probabilities are rescale to avoid underflow. 
#' 
#' The posterior probability that the state at position i is k, \equ{\pi_{i}=k} 
#' given the observed sequence read counts for each parent at position i, \equ{x_{i}}
#' is calculated by:  
#' 
#' \dequ{P(\pi_{i} = k | x_{i}) = \frac{f_{i} b{i}}{\sum\limits_{k=0}^{1} f_{i}^{(k)} b{i}^{(k)}}}
#' 
#' where the super scripts in the denominator are vector indices. 
#' 
#' Finally, we can determine the log likelihood of the whole sequence 
#' of observations by summing up the log of all the scale factors in the 
#' forward probability calculation. 
#' 
#' @return A list of class \code{forward.backward} containing 4 elements corresponding to each 
#' of the four tetrads. Each element is itself a list containing:
#' \describe{
#'     \item The input snp.data
#'     \item spore_number an integer specifying the spore number
#'     \item tetrad.id an integer specifying the tetrad number
#'     \item chr.name a numeric or character specifying the name of the chromosome
#'     \item snp.locations a vector specifying the snp locations along the chromosome
#'     \item p.assign a numeric between 0 and 1 specifying the probability of correct 
#'     sequence assignment
#'     \item p.trans a numeric between 0 and 1 specifying the per base pair recombination 
#'     rate (transition probability)
#'     \item emissions a matrix of 2 columns by length(snp.locations) rows specifying
#'     the emission probabilities (see details)
#'     \item forward a matrix of 2 columns by length(snp.locations) rows giving the scaled
#'     forward probabilities
#'     \item backward a matrix of 2 columns by length(snp.locations) rows giving the scaled
#'     backward probabilities
#'     \item scale a vector of length snp.locations giving the forward scaling factors
#'     \item posterior a matrix of 2 columns by length(snp.locations) rows giving the posterior
#'     probabilities for each parental state
#'     \item states_inferred a vector of length snp.locations giving the state with the
#'     highest posterior. In the event of a tie, a state is picked randomly. 
#'     \item lnL a numeric giving the total likelihood of the data (see details).
#' }
#'  
#' @references Hohenlohe, P.A., S. Bassham, M. Currey, and W.A. Cresko. 2012. Extensive linkage 
#' disequilibrium and parallel adaptive divergence across threespine stickleback genomes. 
#' Phil. Trans. R. Soc. B, 395-408, doi: 10.1098/rstb.2011.0245
#'
#' Drubin, R. S. Eddy, A. Krogh, and G. Mitchison. 1998. Biological Sequence Analysis: Probabilistic Models
#' of proteins and nucleic acids. Cambridge University Press, Cambridge CB2 8RU, UK.
#' 
#' @seealso \code{\link{simulate_coverage}}, \code{\link{est_fwd_back}}
#' 
#' @author Tyler D. Hether
#' 
#' @export tetrad_est_fwd_back
#' 
#' @examples
#' # Example 1: a single tetrad:
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
#' fbres
#' #
#' #
#' # Example 2: 500 simulated tetrads
#' set.seed(1234567) # For reproducability
#' l <- 10
#' rec <- rep(0.01, l-1)
#' n.tetrads <- 500 # number of spores to simulate
#' res <- sim_tetrad(n.tetrads=n.tetrads, l=l, rec=rec, p.assign=0.999, 
#'    mu.rate=0, f.cross=0.8, f.convert=0.9, length.conversion=2, coverage=2.5)
#' # loop through the results list and estimate parental states for each of the n.tetrads
#' fbres2 <- lapply(1:length(res), function(Z){
#'         tetrad_est_fwd_back(snp.dat=res[[Z]], tetrad.id=Z, chr.name="I", p.assign=0.999, p.trans=mean(rec))   
#'     })
#' fbres2

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
        scale <- scaleb <- numeric(n_snps)
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
            T <- matrix(c(1-(displace[i-1]*p.trans), displace[i-1]*p.trans, displace[i-1]*p.trans, 1-(displace[i-1]*p.trans)), byrow=TRUE, ncol=2)
            e <- diag(emissions[i,])
            forward[i,] <- e %*% T %*% forward[i-1,]
            scale[i] = sum(forward[i,])
            forward[i,] = forward[i,]/scale[i]
        }

        # 5) calculate backward probabilities:
        backward[n_snps,] = c(1,1) #1/scale[n_snps]
        for(i in (n_snps-1):1){
            T <- matrix(c(1-(displace[i]*p.trans), displace[i]*p.trans, displace[i]*p.trans, 1-(displace[i]*p.trans)), byrow=TRUE, ncol=2)
            e <- diag(emissions[i+1,])
            backward[i,] <- T %*% e %*% backward[i+1,]
            # backward[i+1,] %*% T %*% e # old (incorrect) way
            # Hardcoded check using algorithm on page 60 of Durbin et al. 1998
            # sum(T[1,1]*emissions[i+1,1]*backward[i+1,1], T[2,1]*emissions[i+1,2]*backward[i+1,2])
            # sum(T[2,1]*emissions[i+1,1]*backward[i+1,1], T[2,2]*emissions[i+1,2]*backward[i+1,2])
            scaleb[i] = sum(backward[i,])
            backward[i,] = backward[i,]/scaleb[i]
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
#'             spore_number=Z, chr.name="I", p.assign=0.999, 
#'             p.trans=mean(rec))
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
        T <- matrix(c(1-(displace[i-1]*p.trans), displace[i-1]*p.trans, displace[i-1]*p.trans, 1-(displace[i-1]*p.trans)), byrow=TRUE, ncol=2)
        e <- diag(emissions[i,])
        forward[i,] <- e %*% T %*% forward[i-1,]
        scale[i] = sum(forward[i,])
        forward[i,] = forward[i,]/scale[i]
    }

    # 5) calculate backward probabilities:
    backward[n_snps,] = c(1,1) #1/scale[n_snps]
    for(i in (n_snps-1):1){
        T <- matrix(c(1-(displace[i]*p.trans), displace[i]*p.trans, displace[i]*p.trans, 1-(displace[i]*p.trans)), byrow=TRUE, ncol=2)
        e <- diag(emissions[i+1,])
        backward[i,] <- T %*% e %*% backward[i+1,]
        # backward[i+1,] %*% T %*% e # old (incorrect) way
        # Hardcoded check using algorithm on page 60 of Durbin et al. 1998
        # sum(T[1,1]*emissions[i+1,1]*backward[i+1,1], T[2,1]*emissions[i+1,2]*backward[i+1,2])
        # sum(T[2,1]*emissions[i+1,1]*backward[i+1,1], T[2,2]*emissions[i+1,2]*backward[i+1,2])
        scaleb[i] = sum(backward[i,])
        backward[i,] = backward[i,]/scaleb[i]
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
    if (!inherits(snp.dat, "snp.recom") & !inherits(snp.dat, "individual.tetrad"))
        stop("Object must be of class 'snp.recom' or 'individual.tetrad'")
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

