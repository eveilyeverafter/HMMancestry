


library(multicore)
library(plyr)
library(dplyr)
library(hmmspore)
library(ggplot2)
library(gridExtra)
library(lattice)
library(microbenchmark)

# Make a temp test dataset of 250 individual diploids
set.seed(1234567) # For reproducability
# simulate a recombination hotspot between the 100th and 101st snp
l=1000; scale = 0.01; snps=1:l
n.spores <- 500 # number of spores to simulate
spores <- sim_en_masse(n.spores=n.spores, l=l, scale=scale, snps=snps, 
p.assign=.999, mu.rate=0.001, f.cross=1, 
    f.convert=0, length.conversion=10, coverage=10)

# make diploids
dat <- lapply(1:250, function(x){
    return(data.frame(Ind=x,Chr="I", Snp=spores[[x]]$snps, 
        p0=spores[[x]]$p0+spores[[x+250]]$p0,
        p1=spores[[x]]$p1+spores[[x+250]]$p1))
    })


dat <- dat[[1]] # pick one


RR <- est_fwd_back_diploid(dat, 0.999, 0.01)
CPP <- c_est_fwd_back_diploid(dat[,3], dat[,4], dat[,5], 0.999, 0.01)
# ^^^ same within round error tol.

microbenchmark(
est_fwd_back_diploid(dat, 0.999, 0.01),
c_est_fwd_back_diploid(dat[,3], dat[,4], dat[,5], 0.999, 0.01)
)





















# How well does the max likelihood function work in finding the true
# values of scale and p.assign in simulation? 

set.seed(1234567) # For reproducability
l=1000; scale = 0.013; snps=1:l;
n.tetrads <- 48 # number of spores to simulate
res <- sim_tetrad(n.tetrads=n.tetrads, l=l, scale=scale, snps=snps, p.assign=0.9964737, 
   mu.rate=0, f.cross=1, f.convert=0, length.conversion=10, coverage=1)


# Convert to proper dataframe:
dat <- lapply(1:length(res), function(x){

	one <- res[[x]][[1]]
	two <- res[[x]][[2]]
	three <- res[[x]][[3]]
	four <- res[[x]][[4]]

	return(rbind(
		data.frame(Tetrad=x, Spore=1, Chr="I", Snp=one$snps, p0=one$p0.assign, p1=one$p1.assign),
		data.frame(Tetrad=x, Spore=2, Chr="I", Snp=two$snps, p0=two$p0.assign, p1=two$p1.assign),
		data.frame(Tetrad=x, Spore=3, Chr="I", Snp=three$snps, p0=three$p0.assign, p1=three$p1.assign),
		data.frame(Tetrad=x, Spore=4, Chr="I", Snp=four$snps, p0=four$p0.assign, p1=four$p1.assign)
	))

	})
dat <- do.call(rbind, dat)

# Parameters range to investigate:
pars <-   expand.grid(p.assign=seq(from=0.995, to=0.997, length.out=20), 
                      scale=seq(from=0.001, to=0.02, length.out=20))

RES <- mclapply(1:dim(pars)[1], function(y){
    res_c <- ddply(dat, .(Tetrad, Spore, Chr), function(x){
       # est_fwd_back(snp.dat=x[,1:6], p.assign=pars[y,1], p.trans= pars[y,2])
       # est_fwd_back(snp.dat=x, p.assign=pars[y,1], scale=pars[y,2])
       c_est_fwd_back(x[,4], x[,5], x[,6], pars[y,1], pars[y,2])
       })
    # Return the lnL for each unique tetrad, spore, & chr
    res_lnl <- ddply(res_c, .(Tetrad, Spore, Chr), function(z){
      return(z$lnL[1])
      })

    return(data.frame(p_assign=pars[y,1], scale=pars[y,2], lnL=res_lnl))
  }, mc.cores=8)

LnL <- unlist(lapply(RES, function(x){
    return(sum(x$lnL.V1))
  }))

RES[[which(LnL==max(LnL))]]
# haldane <- function(d){ (1/2)*(1- exp(-2*d))}
# plot(haldane(c(1:10000)*0.0001),type="l")

LnL <- data.frame(p_assign=pars[,1],scale=pars[,2],LnL)
contourplot(LnL~p_assign*scale, data=LnL, region=TRUE)



