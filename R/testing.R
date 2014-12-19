# testing the simulation and fb algorithms
set.seed(1234567)
l <- 1000
rec <- 0.01
p_a <- .95
p <- make_parents(l)
r <- recombine_index(rep(rec, l-1))
a <- recombine(parents=p, r.index=r, mu.rate=0)
sim_reads <- simulate_coverage(a=a, p_assign=p_a, coverage=0.5)

fbres1 <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=1, 
	chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
fbres2 <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=2, 
	chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
fbres3 <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=3, 
	chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)
fbres4 <- estimate_anc_fwd_back(snp_dat=sim_reads, spore_number=4, 
	chr_name="I", snp_locations=c(1:l), p_assign=p_a, p_trans=rec)

# black = true ancestral assignment (identical by state) from simulation
# non-red color =posterior prob from fb algorithm
# red = location of mutated snps
mut_locs <- c(which(a$chromotids.mutated_snps$p1_1==1),
			  which(a$chromotids.mutated_snps$p1_2==1),
			  which(a$chromotids.mutated_snps$p2_1==0),
			  which(a$chromotids.mutated_snps$p2_2==0))

par(mfrow=c(4,1))
plot(fbres1$snp_locations, fbres1$posterior[,1], xlab="snp", ylab="posterior probability of parent '0' ancestry", main="spore 1", ylim=range(0,1), type="l", col="blue")
points(fbres1$snp_locations, sapply(a$chromotids.recombined$p1_1, switch_values), col="black", type="p", pch=".")
abline(v=mut_locs, col="red")

plot(fbres2$snp_locations, fbres2$posterior[,1], xlab="snp", ylab="posterior probability of parent '0' ancestry", main="spore 2", ylim=range(0,1), type="l", col="orange")
points(fbres2$snp_locations, sapply(a$chromotids.recombined$p1_2, switch_values), col="black", type="p", pch=".")
abline(v=mut_locs, col="red")


plot(fbres3$snp_locations, fbres3$posterior[,1], xlab="snp", ylab="posterior probability of parent '0' ancestry", main="spore 3", ylim=range(0,1), type="l", col="green")
points(fbres3$snp_locations, sapply(a$chromotids.recombined$p2_1, switch_values), col="black", type="p", pch=".")
abline(v=mut_locs, col="red")


plot(fbres4$snp_locations, fbres4$posterior[,1], xlab="snp", ylab="posterior probability of parent '0' ancestry", main="spore 4", ylim=range(0,1), type="l", col="purple")
points(fbres4$snp_locations, sapply(a$chromotids.recombined$p2_2, switch_values), col="black", type="p", pch=".")
abline(v=mut_locs, col="red")

