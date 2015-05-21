# testing the simulation and fb algorithms
set.seed(1234567)
l <- 1000
rec <- 0.01
p_a <- .7
p <- make_parents(l)
r <- recombine_index(rep(rec, l-1))
a <- recombine(parents=p, r.index=r, mu.rate=0.001)
sim_reads <- simulate_coverage(a=a, p_assign=.999, coverage=.2)

fbres <- lapply(c(1:4), function(x,...){
    out <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=x, 
    chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
    return(out)
    })

# Note, if true p_assign >> estimated p_assingn that sensitivity decreases.
# but if true p_assign << estimated p_assign that you get high false positives. 


# black = true ancestral assignment (identical by state) from simulation
# non-red color =posterior prob from fb algorithm
# red = location of mutated snps
mut_locs <- c(which(a$chromotids.mutated_snps$p1_1==1),
			  which(a$chromotids.mutated_snps$p1_2==1),
			  which(a$chromotids.mutated_snps$p2_1==0),
			  which(a$chromotids.mutated_snps$p2_2==0))

par(mfrow=c(4,1))
calls <- list(a$chromotids.recombined$p1_1, a$chromotids.recombined$p1_2,
	a$chromotids.recombined$p2_1, a$chromotids.recombined$p2_2
	)

for(i in 1:4){
	plot(fbres[[i]]$snp_locations, fbres[[i]]$posterior[,1], xlab="snp", ylab="posterior probability of parent '0' ancestry", 
			main=paste("FB spore ", i, sep=""), ylim=range(0,1), type="l", col=c("blue", "orange", "green", "purple"))
	
	points(fbres[[i]]$snp_locations, sapply(calls[[i]], switch_values), col="black", type="p", pch=".")
	abline(v=mut_locs, col="red", lty=3)
}


# # Testing the viterbi algorithm
# V1_1 <- viterbi(snp_dat=sim_reads, spore_number=1, chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
# V1_2 <- viterbi(snp_dat=sim_reads, spore_number=2, chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
# V2_1 <- viterbi(snp_dat=sim_reads, spore_number=3, chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
# V2_2 <- viterbi(snp_dat=sim_reads, spore_number=4, chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)

# dev.new()
# par(mfrow=c(4,1))

# plot(V1_1$snp_locations, V1_1$states, xlab="snp", ylab="ancestry assignment", main="Vitberi spore 1", ylim=range(0,1), type="l", col="blue")
# points(V1_1$snp_locations, a$chromotids.recombined$p1_1, col="black", type="p", pch=".")
# abline(v=mut_locs, col="red", lty=3)

# plot(V1_2$snp_locations, V1_2$states, xlab="snp", ylab="ancestry assignment", main="spore 2", ylim=range(0,1), type="l", col="orange")
# points(V1_2$snp_locations, a$chromotids.recombined$p1_2, col="black", type="p", pch=".")
# abline(v=mut_locs, col="red", lty=3)

# plot(V2_1$snp_locations, V2_1$states, xlab="snp", ylab="ancestry assignment", main="spore 3", ylim=range(0,1), type="l", col="green")
# points(V2_1$snp_locations, a$chromotids.recombined$p2_1, col="black", type="p", pch=".")
# abline(v=mut_locs, col="red", lty=3)

# plot(V2_2$snp_locations, V2_2$states, xlab="snp", ylab="ancestry assignment", main="spore 4", ylim=range(0,1), type="l", col="purple")
# points(V2_2$snp_locations, a$chromotids.recombined$p2_2, col="black", type="p", pch=".")
# abline(v=mut_locs, col="red", lty=3)


# sum(abs(a$chromotids.recombined$p2_2-V2_2$states))
















