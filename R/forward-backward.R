#' Fits forward and backward algorithm
#'
#'

check_est_input <- function(snp_dat, spore_number, chr_name, snp_locations, p_assign, p_trans, ...){

    if (!inherits(snp_dat, "snp.recom"))
        stop("Object must be of class 'snp.recom'. See 'make_snp_data()'")
    check_values(p_assign, 0, 1)
    check_values(p_trans, 0, 1)
    if(length(snp_dat[[1]]$p0.assign)!=length(unique(snp_locations))){
        stop("Each snp needs 1 unique numeric id.")
    }
    if(!spore_number %in% c(1:4)){
        stop("spore_number must be an the integer 1, 2, 3, or 4")
    }

}

estimate_anc_fwd_back <- function(snp_dat, spore_number, chr_name, p_assign, p_trans){

    snp_locations <- snp_dat$snps

    check_est_input(snp_dat, spore_number, chr_name, snp_locations, p_assign, p_trans)
    # 
    displace <- snp_locations[-c(1)]-snp_locations[1:(length(snp_locations)-1)]
      
    n_snps <- length(snp_locations)

    # 1) establish matrices/vectors:
    emissions <- forward <- backward <- posterior <- matrix(0, ncol=2, nrow=n_snps)
    scale <- numeric(n_snps)
    #
    # 2) calculate emission probabilities for all states at all positions
    # data at position i are counts n_ij of allele matching parent j, with n_i = sum(n_ij)
    # emissions[i][j] = (n_i choose n_ij) * (1 - eps)^(n_ij) * (eps)^(n_i - n_ij)
    # if no data at SNP i, emissions[i][j] = 1 for all j
    #
    k0 <- snp_dat[[spore_number]]$p0.assign
    k1 <- snp_dat[[spore_number]]$p1.assign
    n_i <- k0+k1
    p0 <- choose(n_i,k0)*(p_assign^k0)*((1-p_assign)^(n_i-k0))
    p1 <- choose(n_i,k1)*(p_assign^k1)*((1-p_assign)^(n_i-k1))    
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
        transition_matrix_i <- matrix(c(1-(displace[i-1]*p_trans), displace[i-1]*p_trans, displace[i-1]*p_trans, 1-(displace[i-1]*p_trans)), byrow=TRUE, ncol=2)
        forward[i,1] <- sum(transition_matrix_i[1,1]*forward[i-1,1], transition_matrix_i[2,1]*forward[i-1,2]) * emissions[i,1]
        forward[i,2] <- sum(transition_matrix_i[1,2]*forward[i-1,1], transition_matrix_i[2,2]*forward[i-1,2]) * emissions[i,2]
        # where transition_matrix_i is scaled by the distance between position i-1 and i
        scale[i] = sum(forward[i,])
        forward[i,] = forward[i,]/scale[i]
    }

    # 5) calculate backward probabilities:
    backward[n_snps,] = 1/scale[n_snps]
    for(i in (n_snps-1):1){
        transition_matrix_i <- matrix(c(1-(displace[i]*p_trans), displace[i]*p_trans, displace[i]*p_trans, 1-(displace[i]*p_trans)), byrow=TRUE, ncol=2)
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

    # 8) Call states 0 or 1 based on the posterior probs.
    states_inferred <- apply(posterior, 1, function(i){
        if(NaN %in% i){
            warning("Something might not be right. NaNs are present in data.")
        }
        if(check_tie(i)==0){
            return(which(i==max(i))-1)
        }
        if(check_tie(i)==1){
            return(NA)
        }
      })



    out <- list(snp_dat=snp_dat, spore_number=spore_number, chr_name=chr_name, 
             snp_locations=snp_locations, p_assign=p_assign, p_trans=p_trans,
             emissions=emissions, forward=forward, backward=backward, scale=scale, 
             posterior=posterior, states_inferred=states_inferred, lnL=lnL)

    class(out) <- c("list", "forward.backward")
    return(out)

}



#' @method print fwd.back
#' @export 
print.fwd.back <- function(x, ...){
    ## make pretty output format
    ## print to screen
}

