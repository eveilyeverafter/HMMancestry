

# # Make a temp test dataset of 250 individual diploids
# set.seed(1234567) # For reproducability
# # simulate a recombination hotspot between the 100th and 101st snp
# l=100; scale = 0.01; snps=1:l
# n.spores <- 500 # number of spores to simulate
# spores <- sim_en_masse(n.spores=n.spores, l=l, scale=scale, snps=snps, 
# p.assign=.999, mu.rate=0.001, f.cross=0.5, 
#     f.convert=1, length.conversion=10, coverage=1)

# # make diploids
# dat <- lapply(1:250, function(x){
#     return(data.frame(Ind=x,Chr="I", Snp=spores[[x]]$snps, 
#         p0=spores[[x]]$p0+spores[[x+250]]$p0,
#         p1=spores[[x]]$p1+spores[[x+250]]$p1))
#     })


# diploid.dat <- dat[[8]] # pick one

#' @export est_fwd_back_diploid
est_fwd_back_diploid <- function(diploid.dat, p.assign, scale){
    p.trans = scale
    snp.dat <- diploid.dat
    # colnames(diploid.dat) == c("Ind", "Chr", "Snp", "p0", "p1")
    # if(!inherits(snp.dat, "data.frame")) {
    #     if(inherits(snp.dat, "tetrad")) {
    #         snp.dat <- tetrad_to_df(snp.dat)
    #         }
    #     if(inherits(snp.dat, "en.masse")){
    #         snp.dat <- en_masse_to_df(snp.dat)
    #         }
    # }
    # if(!inherits(snp.dat, "data.frame")){
    #     stop("Object snp.dat needs to be of type data.frame 
    #     with columns c('Tetrad', 'Spore', 'Chr', 'Snp', 'p0', 'p1')")
    # } else {
        # tetrad.id <- snp.dat[,1]
        spore_number <- snp.dat[,1]
        chr.name <- snp.dat[,2]
        snp.locations <- snp.dat[,3]        
    # }  

    # check_est_input(p.assign, p.trans)
    # 
    displace <- snp.locations[-c(1)]-snp.locations[1:(length(snp.locations)-1)]
          
    n_snps <- length(snp.locations)

    # 1) establish matrices/vectors:
    emissions <- forward <- backward <- posterior <- matrix(0, ncol=3, nrow=n_snps)
    scale <- scaleb <-  numeric(n_snps)
    #
    # 2) calculate emission probabilities for all states at all positions
    # data at position i are counts n_ij of allele matching parent j, with n_i = sum(n_ij)
    # emissions[i][j] = (n_i choose n_ij) * (1 - eps)^(n_ij) * (eps)^(n_i - n_ij)
    # if no data at SNP i, emissions[i][j] = 1 for all j
    #
    k0 <- snp.dat[,4]
    k1 <- snp.dat[,5]
    n_i <- k0+k1
    p0 <- choose(n_i,k0)*(p.assign^k0)*((1-p.assign)^(n_i-k0))
    ph <- choose(n_i, k0)*(0.5^k0)*(0.5^k1)
    p1 <- choose(n_i,k1)*(p.assign^k1)*((1-p.assign)^(n_i-k1))    
    emissions <- t(rbind(p0, ph, p1))
    #
    # 3) calculate pi[k] (vector of state probabilities at first position -- here use 1/k for each)
    #
    pi_initial <- c(0.25, 0.5, 0.25)
    # 4) calculate forward probabilities and scaling factor:
    forward[1,] <- pi_initial*emissions[1,]
    scale[1] <- sum(forward[1,])
    forward[1,] <- forward[1,]/scale[1]         #re-scale forward probabilities by their sum to avoid underflow
    for(i in 2:n_snps){
        p <- haldane(displace[i-1]*p.trans)
        T <- matrix(c(1-p-p^2, p, p^2,
                      p, 1-2*p, p,
                      p^2, p, 1-p-p^2),3,3,byrow=TRUE)

        # T <- matrix(c(1-(haldane(displace[i-1]*p.trans)), haldane(displace[i-1])*p.trans, haldane(displace[i-1])*p.trans, 1-(haldane(displace[i-1])*p.trans)), byrow=TRUE, ncol=2)
        e <- diag(emissions[i,])
        forward[i,] <- e %*% T %*% forward[i-1,]
        scale[i] = sum(forward[i,])
        forward[i,] = forward[i,]/scale[i]
    }

    # 5) calculate backward probabilities:
    backward[n_snps,] = c(1,1,1) #1/scale[n_snps]
    for(i in (n_snps-1):1){
        p <- haldane(displace[i]*p.trans)  
        T <- matrix(c(1-p-p^2, p, p^2,
              p, 1-2*p, p,
              p^2, p, 1-p-p^2),3,3,byrow=TRUE)      
        # T <- matrix(c(1-(haldane(displace[i]*p.trans)), haldane(displace[i]*p.trans), haldane(displace[i]*p.trans), 1-(haldane(displace[i]*p.trans))), byrow=TRUE, ncol=2)
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
        posterior[i,] <- (forward[i,]*backward[i,])/(forward[i,1]*backward[i,1] + forward[i,2]*backward[i,2] + forward[i,3]*backward[i,3])
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

    out <- data.frame(Ind=spore_number, Chr=chr.name, 
             Snp=snp.locations, emiss0=emissions[,1], emiss1=emissions[,2],
             forward0=forward[,1], forward1=forward[,2], backward0=backward[,1],
             backward1=backward[,2], Fscale=scale, Bscale=scaleb,
             posterior0=posterior[,1], posterior1=posterior[,2],
             states_inferred=states_inferred, lnL=rep(lnL, length(n_snps)))

    class(out) <- c("data.frame")
    return(out)
}

