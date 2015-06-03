

#' @export est_maxLnL
est_maxLnL <- function(dat, ploidy="diploid", initial_p_assign=0.99, initial_scale=1e-04, tolerance=1e-04, n_iterations=30)
{
	# initial iteration
	n = 0 

	# This is the tolerance. If the estimates of x1 and y1 are 
	# less than tol distance apart than stop the algorithm
	tol <- tolerance 			

	# diff is the difference in lnl between the old and new values. 
	# It is initially larger than tol
	diff <- 1		

	# dx is the distance between the focal x1 point and points x0 or x2
	# where x is the p_assign parameter that is to be estimated
	dx <- 1e-04

	# dy is the distancw between the focal y1 point and points y0 or y2
	# where y is the scale parameter that is to be estimated
	dy <- 1e-06

	x <- initial_p_assign; y <- initial_scale # (x,y) is the starting (p_assign, scale) point for the algorithm. 

	# dat is a df with 5 columns: "Ind", "Chr", "Snp", "p0", "p1"
	colnames(dat) <- c("Ind", "Chr", "Snp", "p0", "p1")

	repeat{
		n <- n + 1
		# print(n)
		# For each iteration there are 6 (p_assign, scale) points that need estimated lnls
		# before a jump can be proposed. These points are below. Note that 
		# x is bounded between 0.5 and 1 (exclusive) and 
		# y must always be positive.

		points <- list(
		"x0y0" = c(max((x-dx), 0.50000001),max((y-dy), 1e-09)),	
		"x1y0" = c(x, max((y-dy), 1e-09)), 
		"x0y1" = c(max((x-dx), 0.50000001),y), 
		"x1y1" = c(x,y), 
		"x2y1" = c(min((x+dx), 0.99999999),y), 
		"x1y2" = c(x, (y+dy))
		)

		lnls <- lapply(points, function(fun){
			# Call the hmm and get the lnl for each unique individual by chromosome combination
			res_c <- ddply(dat, .(Ind, Chr), function(x){
		    	   if(ploidy=="haploid") c_est_fwd_back(x[,3], x[,4], x[,5], fun[1], fun[2])
		           if(ploidy=="diploid") c_est_fwd_back_diploid(x[,3], x[,4], x[,5], fun[1], fun[2])
		           })
		    # Return the lnL for each unique individual and chromosome
		     res_lnl <- ddply(res_c, .(Ind, Chr), function(z){
		          return(z$lnL[1])
		          })
		     # Sum over the above lnl to get the total lnl for all data
		     lnl <- sum(res_lnl$V1)
		     return(lnl)	
		})

		# Here's the Hessian matrix
		H <- matrix(c(
			(lnls$x0y1 + lnls$x2y1 -2*lnls$x1y1)/(dx*dx),
			(lnls$x0y0 + lnls$x1y1 - lnls$x0y1 - lnls$x1y0)/(dx*dy),
			(lnls$x0y0 + lnls$x1y1 - lnls$x0y1 - lnls$x1y0)/(dx*dy),
			(lnls$x1y0 + lnls$x1y2 -2*lnls$x1y1)/(dy*dy)
			),2,2,byrow=TRUE)
		# Here's the gradient
		Grad <- c(
			(lnls$x2y1 - lnls$x1y1)/dx, 
			(lnls$x1y2 - lnls$x1y1)/dy
			)
		# The new x,y point
		xy <- c(x,y) - solve(H)%*%Grad

		# x and y are bounded: x between 0.5 and 1 (exclusive) and y between 0 and 1.
		# Adjust the new xy point if needed to their dx (dy) displacements:
		if(xy[1]<=0.5)
		{
			xy[1] = x-dx
		}
		if(xy[1]>=1)
		{
			xy[1] = x+dx
		}
		if(xy[2]<=0)
		{
			xy[2] = y-dy
		}
		if(xy[2]>=1)
		{
			xy[2] = y+dy
		}

		# The Euclidean distance between the old and new point is:
		diff <- dist(matrix(c(x,y,xy[1],xy[2]),2,2, byrow=TRUE))

		# Stop the algorithm if the number of iterations is exceeded or 
		# the max lnl values do not change much (within tolerance)
		if(n==n_iterations || diff<tol)
		{
			if(n==n_iterations)
			{
				print(paste("maximum lnl estimate not reached after ", n, "iterations", sep=""))
				break
			}
			print(paste("maximum lnl estimate reached after ", n, " iterations", sep=""))
			break
		}
		# Update parameters
		x <- xy[1]
		y <- xy[2]
	}

	return(data.frame(p_assign_hat=xy[1], scale_hat=xy[2]))
}

