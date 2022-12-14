
stirling_approx <- function(n) 1/2 * log( 2 * pi ) + ( 1/2 + n ) * log(n) - n

detect_outliers <- function(x, method = "iqr", ... )
{
  
  if( method == "iqr" ) iqr_outliers( x = x, ...)
  if( method == "z_score") z_outliers( x = x, ...)
  if( method == "ueda" ) ueda_outliers( x = x, ... )
  
}

iqr_outliers <- function( x, alpha = 1.5 )
{
  q1  <- as.numeric( quantile(x, .25) )
  q3  <- as.numeric( quantile(x, .75) )
  iqr <- q3 - q1
  
  lower_bound <- q1 - alpha * iqr
  upper_bound <- q3 + alpha * iqr
  
  x_old <- sort( x )
  x_new <- x_old[ x_old > lower_bound & x_old < upper_bound ]
  
  s_detected <- sum( !(x_old > lower_bound & x_old < upper_bound) )
  
  list( alpha = alpha, s_detected = s_detected, x = x_old, x_new = x_new, flag = as.numeric( !x_old %in% x_new) )
}

z_outliers <- function( x, alpha = 2 )
{
  x_sd   <- sd( x, na.rm = T )
  x_mean <- mean( x, na.rm = T ) 
  
  lower_bound <- x_mean - alpha * x_sd
  upper_bound <- x_mean + alpha * x_sd
  
  x_old <- sort( x )
  x_new <- x_old[ x_old > lower_bound & x_old < upper_bound ]
  
  s_detected <- sum( !(x_old > lower_bound & x_old < upper_bound) )
  
  list( alpha = alpha, s_detected = s_detected, x = x_old, x_new = x_new, flag = as.numeric( !x_old %in% x_new) )
}

ueda_outliers <- function( x, s_max = NULL )
{
  n_obs <- length( x )
  x_ord <- scale( sort(x) )[,1]
  
  if( is.null(s_max) ) s_max <- ifelse( n_obs %% 2 == 0, (n_obs / 2) - 1, 0.5 * (n_obs - 1) ) else s_max <- s_max
  
  ut_matrix  <- matrix(0, ncol = s_max + 1, nrow = s_max + 1 )
  n_elements <- length( ut_matrix ) 
  
  for( e in 1:n_elements )
  {
    # update row and column indices
    i <- e %% (s_max + 1) + as.numeric( e %% (s_max + 1) == 0 ) * (s_max + 1) 
    j <- floor((e - 1) / (s_max + 1) ) + 1 #((ii - 1) %/% (s_max + 1)) + 1
    
    # total number of proposed outliers
    s <- (i-1) + (j-1)
    
    # remove proposed outliers from either end, based on current i and j indices
    # compute new length and standard deviation for truncated set
    zs <- x_ord[(1 + (i - 1) ):(n_obs - (j - 1) )]
    ns <- length( zs )
    sd <- sqrt( sum( ( zs - mean(zs) )^2) / ns )
    
    # Compute ut and assign to current i,j location in matrix
    ut_matrix[i,j] <- ns * log( sd ) + ( sqrt(2) * stirling_approx(ns) ) / ns * s
  } 
  
  # Results and outputs
  mat_index    <- which( ut_matrix == min(ut_matrix), arr.ind = T)
  n_left_right <- mat_index - 1
  s_detected   <- sum( n_left_right )
  
  # Remove detected outliers  
  x_old <- sort( x )
  x_new <- x_old[(1 + n_left_right[1] ):(n_obs - n_left_right[2]) ]
  
  list( ut = round(ut_matrix, 3), ut_min = round(ut_matrix[mat_index], 3), s_detected = s_detected, 
        x = x_old, x_new = x_new, flag = as.numeric( !x_old %in% x_new) )
  
}
