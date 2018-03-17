########## AMS 206: R code for problem 2(B)#########

y <- c( 409., 400., 406., 399., 402., 406., 401., 403., 401., 403., 398., 
        403., 407., 402., 401., 399., 400., 401., 405., 402., 408., 399., 399., 
        402., 399., 397., 407., 401., 399., 401., 403., 400., 410., 401., 407., 
        423., 406., 406., 402., 405., 405., 409., 399., 402., 407., 406., 413., 
        409., 404., 402., 404., 406., 407., 405., 411., 410., 410., 410., 401., 
        402., 404., 405., 392., 407., 406., 404., 403., 408., 404., 407., 412., 
        406., 409., 400., 408., 404., 401., 404., 408., 406., 408., 406., 401., 
        412., 393., 437., 418., 415., 404., 401., 401., 407., 412., 375., 409., 
        406., 398., 406., 403., 404. )

print( n <- length( y ) )

# [1] 100

log.likelihood <- function( mu, sigma, nu, y ) {
  
  n <- length( y )
  
  ll <- ( n * lgamma( ( nu + 1 ) / 2 ) - n * log( sigma ) 
          - n * lgamma( nu / 2 ) - n * log( nu ) / 2  - ( ( nu + 1 ) / 2 ) *
            sum( log( 1 + ( ( y - mu ) / sigma )^2 / nu ) ) )
  
  return( ll )
  
}

log.likelihood( 404.3, 3.7, 4.0, y )

# [1] -251.496

log.likelihood.vector.input <- function( theta, y ) {
  
  n <- length( y )
  
  mu <- theta[ 1 ]
  
  sigma <- theta[ 2 ]
  
  nu <- theta[ 3 ]
  
  ll <- ( n * lgamma( ( nu + 1 ) / 2 ) - n * log( sigma ) 
          - n * lgamma( nu / 2 ) - n * log( nu ) / 2  - ( ( nu + 1 ) / 2 ) *
            sum( log( 1 + ( ( y - mu ) / sigma )^2 / nu ) ) )
  
  return( ll )
  
}

print( mu.0 <- mean( y ) )

# [1] 404.59

nu.0 <- 6.

print( sigma.0 <- sd( y ) * sqrt( ( nu.0 - 2 ) / nu.0 ) )

# [1] 5.280158

theta.0 <- c( mu.0, sigma.0, nu.0 )

print( log.likelihood.maximization <- optim( theta.0, 
                                             log.likelihood.vector.input, y = y, method = 'L-BFGS-B',
                                             control = list( fnscale = -1 ), hessian = T,
                                             lower = c( 0, 0, 0 ) ) )

# $par
# [1] 404.278406   3.700285   3.008206

# $value
# [1] -250.8514

# $counts
# function gradient 
#       15       15 

# $convergence
# [1] 0

# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

# $hessian
#            [,1]       [,2]       [,3]
# [1,] -4.7208971  0.1850944  0.1138392
# [2,]  0.1850944 -7.7689323  1.9642652
# [3,]  0.1138392  1.9642652 -1.8484565

print( theta.hat.mle <- log.likelihood.maximization$par )

# [1] 404.278406   3.700285   3.008206

mu.hat.mle <- theta.hat.mle[ 1 ]

sigma.hat.mle <- theta.hat.mle[ 2 ]

nu.hat.mle <- theta.hat.mle[ 3 ]

print( covariance.matrix.theta.hat <- 
         solve( - log.likelihood.maximization$hessian ) )

#            [,1]       [,2]       [,3]
# [1,] 0.21288379 0.01146797 0.02529714
# [2,] 0.01146797 0.17662445 0.18839652
# [3,] 0.02529714 0.18839652 0.74274971

print( correlation.matrix.theta.hat <- 
         cov2cor( covariance.matrix.theta.hat ) )

#            [,1]       [,2]       [,3]
# [1,] 1.00000000 0.05914118 0.06361785
# [2,] 0.05914118 1.00000000 0.52014714
# [3,] 0.06361785 0.52014714 1.00000000

print( se.hat.theta.hat <- 
         sqrt( diag( covariance.matrix.theta.hat ) ) )

# [1] 0.4613933 0.4202671 0.8618293

####################### end of numerical analysis; ######################
########## graphical analysis from here to the end of the file ##########

n.grid <- 25

c( theta.hat.mle[ 1 ] - 3 * se.hat.theta.hat[ 1 ], 
   theta.hat.mle[ 1 ] + 3 * se.hat.theta.hat[ 1 ] )

# [1] 402.8942 405.6626

mu.grid <- seq( theta.hat.mle[ 1 ] - 3 * se.hat.theta.hat[ 1 ], 
                theta.hat.mle[ 1 ] + 3 * se.hat.theta.hat[ 1 ], length = n.grid )

c( theta.hat.mle[ 2 ] - 3 * se.hat.theta.hat[ 2 ], 
   theta.hat.mle[ 2 ] + 3 * se.hat.theta.hat[ 2 ] )

# [1] 2.439484 4.961086

sigma.grid <- seq( theta.hat.mle[ 2 ] - 3 * se.hat.theta.hat[ 2 ], 
                   theta.hat.mle[ 2 ] + 3 * se.hat.theta.hat[ 2 ], length = n.grid )

c( theta.hat.mle[ 3 ] - 3 * se.hat.theta.hat[ 3 ], 
   theta.hat.mle[ 3 ] + 3 * se.hat.theta.hat[ 3 ] )

# [1] 0.4227177 5.5936944

nu.grid <- seq( theta.hat.mle[ 3 ] - 3 * se.hat.theta.hat[ 3 ], 
                theta.hat.mle[ 3 ] + 3 * se.hat.theta.hat[ 3 ], length = n.grid )

#########################################################################

########### beginning of code block you need to modify twice, ###########
########### to make the contour and perspective plots for     ###########
########### ( mu, nu ) and for ( sigma, nu )                  ###########

log.likelihood.mu.sigma.grid <- matrix( NA, n.grid, n.grid )

for ( i in 1:n.grid ) {
  
  for ( j in 1:n.grid ) {
    
    log.likelihood.mu.sigma.grid[ i, j ] <- 
      log.likelihood.vector.input( c( mu.grid[ i ],
                                      sigma.grid[ j ], nu.hat.mle ), y )
    
  }
  
}

pdf( 'plot1_THT2.pdf' )  # you need to pick a filename,


par( mfrow = c( 2, 2 ) )

contour( mu.grid, sigma.grid, log.likelihood.mu.sigma.grid,
         xlab = 'mu', ylab = 'sigma', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( mu.grid, sigma.grid, log.likelihood.mu.sigma.grid,
       xlab = 'mu', ylab = 'sigma', theta = 30, phi = 30,
       zlab = 'log likelihood', col = 'red' )

contour( mu.grid, sigma.grid, exp( log.likelihood.mu.sigma.grid ),
         xlab = 'mu', ylab = 'sigma', nlevels = 10, col = 'blue',
         lwd = 2 )

persp( mu.grid, sigma.grid, exp( log.likelihood.mu.sigma.grid ),
       xlab = 'mu', ylab = 'sigma', theta = 30, phi = 30,
       zlab = 'likelihood', col = 'red' )

par( mfrow = c( 1, 1 ) )

dev.off( )        # you need to uncomment and run these two lines
#                     to send any subsequent plots to your screen

########### end of code block you need to modify twice        ###########