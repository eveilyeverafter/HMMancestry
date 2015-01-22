source("helper-hmmspore.R")

context("Testing tetrad dependent code")

## Testing make_parents function

p1 <- 10 # good
p1v <- c(p1, p1) # should fail
p2 <- 1 # should fail
p3 <- "10" # should fail

test_that("Initializing parents returns correct format.", {
	expect_that(make_parents(p1), is_a("parent.genomes"))
	expect_equal(length(make_parents(p1)$snps), 10)
	expect_that(make_parents(c(p1v)), throws_error("The number of loci to simulate needs to be a vector of length 1"))
	expect_that(make_parents(p2), throws_error("The number of loci to simulate needs to be >= 2"))
	expect_that(make_parents(p3), throws_error("The number of loci to simulate needs to be numeric"))
	})

## Testing recombine_index function

p4 <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1) # good
p5 <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.4, 0.3, 0.2, 0.1) # should fail
p6 <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, "0.1") # should fail

set.seed(1)
test_that("Recombination point sampling works as expected.", {
	expect_that(recombine_index(p4), equals(c(0,0,0,1,0,1,1,0,0)))
	expect_that(recombine_index(p5), throws_error("The maximum of 'a' needs to be <0.5"))
	expect_that(recombine_index(p6), throws_error("The vector 'p.trans' needs to be numeric"))
	})

## Testing recombine function

set.seed(1)
l <- 50 
rec <- 0.2
r <- recombine_index(rep(rec, l-1))
p <- make_parents(l)

test_that("Main recombination function's params work.", {
		expect_that(length(recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=1, length.conversion=10)), equals(8)) # works
		expect_that(recombine(parents=p, r.index=r, mu.rate=-1, f.cross=.5, f.convert=1, length.conversion=10), throws_error("non-positive probability")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate="a", f.cross=.5, f.convert=1, length.conversion=10), throws_error("in 1 - mu : non-numeric argument to binary operator")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=1.1, f.cross=.5, f.convert=1, length.conversion=10), throws_error("non-positive probability")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=0, f.cross=-1, f.convert=1, length.conversion=10),throws_error("non-positive probability")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=0, f.cross="a", f.convert=1, length.conversion=10),throws_error("1 - f.cross : non-numeric argument to binary operator")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=0, f.cross=1.1, f.convert=1, length.conversion=10), throws_error("non-positive probability")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=-1, length.conversion=10),throws_error("non-positive probability")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert="a", length.conversion=10), throws_error("in 1 - f.convert : non-numeric argument to binary operator")) # fails
		expect_that(recombine(parents=p, r.index=r, mu.rate=0, f.cross=.5, f.convert=1.1, length.conversion=10),throws_error("non-positive probability")) # fails
})








