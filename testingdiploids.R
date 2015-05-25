


library(multicore)
library(plyr)
library(dplyr)
library(hmmspore)
library(ggplot2)
library(gridExtra)
library(lattice)
library(microbenchmark)

# # Make a temp test dataset of 250 individual diploids
# set.seed(1234567) # For reproducability
# l=100; scale = 0.1; snps=c(1:l)
# n.spores <- 250 # number of spores to simulate
# # spores <- sim_en_masse(n.spores=n.spores, l=l, scale=scale, snps=snps, 
# # p.assign=.965, mu.rate=0, f.cross=1, 
# #     f.convert=0, length.conversion=10, coverage=5)


set.seed(1234567) # For reproducability
l=1000; scale = 0.01; snps=(1:l);
n.tetrads <- 100 # number of spores to simulate
res <- sim_tetrad(n.tetrads=n.tetrads, l=l, scale=scale, snps=snps, p.assign=0.95, 
   mu.rate=0.001, f.cross=0.8, f.convert=0.2, length.conversion=10, coverage=0.25)


# make diploids
dat <- lapply(1:(n.tetrads/2), function(x){
  true_states1 <- res[[x]][[1]]$states_given
  true_states2 <- res[[x+(n.tetrads/2)]][[1]]$states_given

    return(data.frame(Ind=x,Chr="I", Snp=res[[x]]$snps, 
        p0=res[[x]][[1]]$p0.assign+res[[x+(n.tetrads/2)]][[1]]$p0.assign,
        p1=res[[x]][[1]]$p1.assign+res[[x+(n.tetrads/2)]][[1]]$p1.assign,
        states_known=true_states1+true_states2))
    })


# dat <- dat[[1]] # pick one


# RR <- est_fwd_back_diploid(dat, 0.999, 0.01)
# CPP <- c_est_fwd_back_diploid(dat[,3], dat[,4], dat[,5], 0.999, 0.01)
# # ^^^ same within round error tol.

# microbenchmark(
# est_fwd_back_diploid(dat, 0.999, 0.01),
# c_est_fwd_back_diploid(dat[,3], dat[,4], dat[,5], 0.999, 0.01)
# )


# Try to estimate the scale and p.assign parameters from the data:

dat <- do.call(rbind, dat)
colnames(dat) <- c("Sample_ID", "Chr", "Snp", "p0", "p1","states_known")

pars <-   expand.grid(p.assign=seq(from=0.85, to=0.97, length.out=21), 
                      scale=seq(from=0.0001, to=0.2, length.out=21))

RES <- mclapply(1:dim(pars)[1], function(y){
    res_c <- ddply(dat, .(Sample_ID, Chr), function(x){
       # est_fwd_back(snp.dat=x[,1:6], p.assign=pars[y,1], p.trans= pars[y,2])
       c_est_fwd_back_diploid(x[,3], x[,4], x[,5], pars[y,1], pars[y,2])
       })
    # print(dim(res_c))
    # Return the lnL for each unique tetrad, spore, & chr
    res_lnl <- ddply(res_c, .(Sample_ID, Chr), function(z){
      return(z$lnL[1])
      })

    return(data.frame(p_assign=pars[y,1], scale=pars[y,2], lnL=res_lnl))
  }, mc.cores=8)

LnL <- unlist(lapply(RES, function(x){
    return(sum(x$lnL.V1))
  }))

RES[[which(LnL==max(LnL))]]



MAX <- head(RES[[which(LnL==max(LnL))]])[1,1:2]
# haldane <- function(d){ (1/2)*(1- exp(-2*d))}
# plot(haldane(c(1:10000)*0.0001),type="l")


LnL <- data.frame(p_assign=pars[,1],scale=pars[,2],LnL)


# pdf("./likelihood_suface_coarse.pdf")


contourplot(LnL~p_assign*scale, data=LnL, region=TRUE,cuts=10)


# dev.off()


resMax <- ddply(dat, .(Sample_ID, Chr), function(x){
   # est_fwd_back(snp.dat=x[,1:6], p.assign=0.9973333, p.trans= 0.012)
   cbind(c_est_fwd_back_diploid(x[,3], x[,4], x[,5], as.numeric(MAX[1]), as.numeric(MAX[2])),states_given=x[,6])
   })

resMax[,c(1,3,4,5,20,22)]




MAX

error_frequency <- length(which((resMax[,20] - resMax[,22])!=0))/dim(resMax)[1]
error_frequency










