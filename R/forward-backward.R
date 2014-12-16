#' Fits forward and backward algorithm
#'
#'
estimate_anc_fwd_back <- function(snp_dat, p_assign, p_trans,...){

    check_est_input(snp_dat, p_assing, p_trans)

    p_1 <- lapply(snp_dat$parent_1, function(x)
                  anc_fwd_back(x, p_assign, p_trans))
    p_2 <- lapply(sp_dat$parent)

    ## collect all the stuff
    ## build an output
    
}


check_est_input <- function(snp_dat, p_assign, p_trans, ...){

    if (!inherits(snp_dat, "snp.recom"))
        stop("Object must be of class 'snp.recom'. See 'make_snp_data()'")

    if (p_assign >= 1 | p_assign < 0)
        stop("p_assign must be between 0 and 1")

    if (p_trans >= 1 | p_trans < 0)
        stop("p_trans must be between 0 and 1")
}

#' @method print fwd.back
#' @export 
print.fwd.back <- function(x, ...){
    ## make pretty output format
}

