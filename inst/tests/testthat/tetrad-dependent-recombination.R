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
	expect_that(make_parents(p3), throws_error("The number2 of loci to simulate needs to be numeric"))
	})

