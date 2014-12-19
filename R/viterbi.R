#' Fits Viterbi algorithm
#'
estimate_anc_viterbi <- function(snp_dat, p_assign, p_trans, ...){

   check_est_input(snp_dat, p_assing, p_trans)
   
}


viterbi <- function(snp_dat, spore_number, chr_name, snp_locations, p_assign, p_trans){
 
    check_est_input(snp_dat, spore_number, chr_name, snp_locations, p_assign, p_trans)
    # Distance between each snp
    displace <- snp_locations[-c(1)]-snp_locations[1:(length(snp_locations)-1)]
    
    # n <- snp_locations
    n_snps <- length(snp_locations)

    k0 <- snp_dat[[spore_number]]$p0.assign
    k1 <- snp_dat[[spore_number]]$p1.assign
    n_i <- k0+k1
    p0 <- choose(n_i,k0)*(p_assign^k0)*((1-p_assign)^(n_i-k0))
    p1 <- choose(n_i,k1)*(p_assign^k1)*((1-p_assign)^(n_i-k1)) 
    B = log2(rbind(p0, p1)) # where Bij stores the prob of observing o_j from state s_i

    initial <- log2(c(0.5, 0.5))
    T1 <- T2 <- matrix(NA, nrow=2, ncol=n_snps)

    # t = 1
    for(i in 0:1){
        T1[(i+1),1] <- initial[(i+1)]+B[(i+1),1]
        T2[(i+1),1] <- i
    }
    # 1 > t >= T
    k = 1:2 # states "1" and "2"
    for(i in 2:n_snps){
        # Per snp recombination rate estimate
        Pi <- matrix(c(1-(displace[i-1]*p_trans), displace[i-1]*p_trans, displace[i-1]*p_trans, 1-(displace[i-1]*p_trans)), byrow=TRUE, ncol=2)
        A = log2(Pi)
        for(j in 1:2){
            T1[j,i] <- max(T1[k,i-1]+A[k,j]+B[k,i])
            tmp <- names(which(T1[k,i-1]+A[k,j]+B[k,i]==T1[j,i]))
            # if(length(tmp)==2){
            #     warning(paste("Transition from either state '0' or '1' and emitting state ", j-1 , " is equally likely at snp ", i, ".  Sequence coverage may be too low for Viterbi.", sep=""))
            # }
            if(length(tmp)==2){
                T2[j,i] <- 0.5
            }
            if(length(tmp)==1){
                if(tmp=="p0"){
                    T2[j,i] <- 0  
                }
                if(tmp=="p1"){
                    T2[j,i] <- 1
                }  
            }   
        }
    }
    z_T <- apply(T1, 2, max)
    states <- numeric(n_snps)

    for(i in length(z_T):1){
        states[i] <- T2[which(T1[,i]==z_T[i]),i]
    }

    out <- list(snp_dat=snp_dat, spore_number=spore_number, chr_name=chr_name,
             snp_locations=snp_locations, p_assign=p_assign, p_trans=p_trans, T1=T1, T2=T2, states=states)
    class(out) <- c("list", "viterbi")
    return(out)
    }






