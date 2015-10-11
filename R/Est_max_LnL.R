

#' @export est_maxLnL
est_maxLnL <- function(dat, ploidy="diploid", initial_p_assign="NULL", initial_scale=5e-05, tolerance=1e-04, n_coarse_steps=5, n_iterations=30, plot=FALSE)
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
	dx <- dx_master <-  1e-04

	# dy is the distancw between the focal y1 point and points y0 or y2
	# where y is the scale parameter that is to be estimated
	dy <- dy_master <-  1e-05

	x <- initial_p_assign; y <- initial_scale # (x,y) is the starting (p_assign, scale) point for the algorithm. 

	# dat is a df with 5 columns: "Ind", "Chr", "Snp", "p0", "p1"
	colnames(dat) <- c("Ind", "Chr", "Snp", "p0", "p1")
	dat <- dat[,1:5]

	# If x or y == NULL, do a coarse scale grid to estimate appropriate starting values.
	# This is either a 5 point 1D grid (if either x or y is NULL) or a 5x5 2D grid if 
	# both x and y are NULL.
	if(x=="NULL" || y=="NULL")
	{
		if(x=="NULL")
		{
			x <- seq(from=0.95, to=0.999, length.out=n_coarse_steps)
		} 
		if(y=="NULL")
		{
			y <- seq(from=1e-06, to=1e-04, length.out=n_coarse_steps)
		}
		pars <-   expand.grid(x=x, y=y) 			# This is the grid
    
	    # Estimate the most likely set of starting values 
	    coarse_res <- lapply(1:dim(pars)[1], function(yy){
	        res_c <- ddply(dat, .(Ind, Chr), function(xx){
	           if(ploidy=="haploid")
	           {
	           		return(c_est_fwd_back(xx[,3], xx[,4], xx[,5], pars[yy,1], pars[yy,2]))
	           	}
	           if(ploidy=="diploid")
	           {
	           		return(c_est_fwd_back_diploid(xx[,3], xx[,4], xx[,5], pars[yy,1], pars[yy,2]))
	           	}
	           
	           })
	        # Return the lnL for each unique tetrad, spore, & chr
	        res_lnl <- ddply(res_c, .(Ind, Chr), function(z){
	          return(z$lnL[1])
	          })
	      	return(data.frame(p_assign=pars[yy,1], scale=pars[yy,2], lnL=res_lnl))
	    })

	    LnL <- unlist(lapply(coarse_res, function(x){
	        return(sum(x$lnL.V1))
	      }))
	    
	    # The maximum likelihood estimate for p_assign and scale are:     
	    MAX <- head(coarse_res[[which(LnL==max(LnL))]])[1,1:2]
		print(paste("Coarse estimates are ", MAX[1], " and ", MAX[2], sep=""))
		x <- as.numeric(MAX[1])
		y <- as.numeric(MAX[2])

	    # For visualizing:
	    if(plot==TRUE)
	    {
		    LnL <- data.frame(p_assign=pars[,1],scale=pars[,2],LnL)

		    if(initial_p_assign=="NULL" && initial_scale!="NULL")
		    {
	    		print(plot(LnL$p_assign, LnL$LnL, type="b", xlab="p_assign", ylab=paste("LnL given scale = ", y, sep="")))
		    }
		    if(initial_p_assign!="NULL" && initial_scale=="NULL")
		    {
	    		print(plot(LnL$scale, LnL$LnL, type="b", xlab="scale", ylab=paste("LnL given p_assign = ", x, sep="")))
		    }
		    if(initial_p_assign=="NULL" && initial_scale=="NULL")
		    {
	    		print(contourplot(LnL~p_assign*scale, data=LnL, region=TRUE, cuts=6))
		    }
	    }

	}

	# Perform the hill climbing procedure to fine-tune the parameters
	repeat{
		if(n_iterations==0) break
		n <- n + 1
		print(paste("Processing iteration ", n, " out of a max of ", n_iterations, ".",sep=""))
		dx <- dx_master/n
		dy <- dy_master/n

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
			res_c <- ddply(dat[,c(1:5)], .(Ind, Chr), function(xx){
		    	   if(ploidy=="haploid")
		    	   {
		    	   		return(c_est_fwd_back(xx[,3], xx[,4], xx[,5], fun[1], fun[2]))
		    	   }
		           if(ploidy=="diploid")
		           {
		           		return(c_est_fwd_back_diploid(xx[,3], xx[,4], xx[,5], fun[1], fun[2]))
		           }
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
			xy[1] = max(x-1e-06, 0.50001)
		}
		if(xy[1]>=1)
		{
			xy[1] = min(x+1e-06, 0.9999)
		}
		if(xy[2]<=0)
		{
			xy[2] = max(y-1e-06, 0.00001)
		}
		if(xy[2]>=1)
		{
			xy[2] = min(y+1e-06, 0.9999)
		}

		# print(xy)

		# See if the new xy point has a higher lnl than all 6 points. 
		# If it is greater, accept the move. 
		new_xy_lnl <- ddply(dat[,c(1:5)], .(Ind, Chr), function(xx){
		   if(ploidy=="haploid")
		   {
		   		return(c_est_fwd_back(xx[,3], xx[,4], xx[,5], xy[1], xy[2]))
		   }
		   if(ploidy=="diploid")
		   {
		   		return(c_est_fwd_back_diploid(xx[,3], xx[,4], xx[,5], xy[1], xy[2]))
		   }
		   })
		# Return the lnL for each unique individual and chromosome
		res_new_lnl <- ddply(new_xy_lnl, .(Ind, Chr), function(z){
		     return(z$lnL[1])
		     })
		# Sum over the above lnl to get the total lnl for all data
		new_lnl <- sum(res_new_lnl$V1)

# Troubleshooting:
# xs <- c(lapply(points, "[", 1L), xy[1]); ys <- c(lapply(points, "[", 2L), xy[2]); DAT <- data.frame(xs=unlist(xs), ys=unlist(ys), l=as.numeric(c(unlist(lnls), new_lnl))); ggplot(data=DAT, aes(x=xs, y=ys, colour=factor(l))) + geom_point()


		# Is the new xy's lnl > all 6 points? If so, accept it. 
		if(new_lnl > rev(sort(unlist(lnls)))[1])
		{
			# The Euclidean distance between the old and new point is:
			diff <- dist(matrix(c(x,y,xy[1],xy[2]),2,2, byrow=TRUE))
		
		} else {
			# If one or more of the 6 points have a greater lnl than the 
			# proposed xy point, then pick the point with the highest lnl 
			xy <- points[[as.numeric(which(unlist(lnls)==max(unlist(lnls))))]]
			# print("One of the 6 points has a greater lnl than the proposed point")
		}
		
		# Update parameters
		x <- xy[1]
		y <- xy[2]

		# Stop the algorithm if the number of iterations is exceeded or 
		# the max lnl values do not change much (within tolerance)
		if(n==n_iterations || diff<tol)
		{
			if(n==n_iterations)
			{
				print(paste("maximum lnl estimate not reached after ", n, " iterations", sep=""))
				break
			}
			print(paste("maximum lnl estimate reached after ", n, " iterations", sep=""))
			break
		}
			
	}

	return(data.frame(p_assign_hat=xy[1], scale_hat=xy[2],n_iterations=n))
}

