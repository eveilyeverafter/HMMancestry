#' Fits Viterbi algorithm
#'
estimate_anc_viterbi <- function(snp_dat, p_assign, p_trans, ...){

   check_est_input(snp_dat, p_assing, p_trans)
   
}


       BinomViterbi <- function(sporeNumber,...){
                n=pos.dat.frameTotal[,1+sporeNumber]
                # p=0.99 # prob of sucess
                kY=pos.dat.frameY[,1+sporeNumber]
                kD=pos.dat.frameD[,1+sporeNumber]
                pY <- choose(n,kY)*(p^kY)*((1-p)^(n-kY))
                pD <- choose(n,kD)*(p^kD)*((1-p)^(n-kD))
                B = log2(rbind(pY, pD)) # where Bij stores the prob of observing o_j from state s_i
                # t_p <- 90.5/12800000 # per bp recombination rate
                # Distance between each snp
                displace <- pos.dat.frameTotal$pos[-c(1)]-pos.dat.frameTotal$pos[1:(length(pos.dat.frameTotal$pos)-1)]
 
                # A = log2(Pi)
                initial <- log2(c(0.5, 0.5))
                T1 <- T2 <- matrix(NA, nrow=2, ncol=length(kY))

                # t = 1
                for(i in 1:2){
                    T1[i,1] <- initial[i]+B[i,1]
                    T2[i,1] <- i
                }
                # 1 > t >= T
                k = 1:2
                for(i in 2:length(kY)){
                    # Per snp recombination rate estimate
                    Pi <- matrix(c(1-(displace[i-1]*t_p), displace[i-1]*t_p, displace[i-1]*t_p, 1-(displace[i-1]*t_p)), byrow=TRUE, ncol=2)
                    A = log2(Pi)
                    for(j in 1:2){

                        T1[j,i] <- max(T1[k,i-1]+A[k,j]+B[k,i])
                        tmp <- names(which(T1[k,i-1]+A[k,j]+B[k,i]==T1[j,i]))
                        if(tmp=="pY"){
                            T2[j,i] <- 1  
                        }
                        if(tmp=="pD"){
                            T2[j,i] <- 2
                        }     
                    }
                }
                z_T <- apply(T1, 2, max)
                states <- numeric(length(z_T))

                for(i in length(z_T):2){
                    states[i] <- T2[which(T1[,i]==z_T[i]),i]
                }
                states[1] <- states[2]
                return((states*-1)+2)
        }


        if(useViterbi==TRUE){
            states <- mclapply(1:4, BinomViterbi, mc.cores=4)
        }
