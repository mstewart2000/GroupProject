########################################################################################################################
#### MASTER CODE defineing the initial feedback matrix based on the given feedback regime ##############################
#### Feedbacks scenarios are as follows (i.e. based on S_{ij} values from Fukami & Nakajima, 2013, Ecology Letters) ####
########################################################################################################################
#### This function takes in an integer between 1 - 5 representing different feedback regimes as follows: 
#### (1) NO_FEEDBACK = 1;         All S_{ij} = 1
#### (2) POSITIVE_FEEDBACK = 2;   Conspecific S_{ii} are between 1 - 1.5 && Heterospecific S_{ij} = 1  
#### (3) NEGATIVE_FEEDBACK = 3;   Conspecific S_{ii} are between 0.5 - 1 && Heterospecific S_{ij} = 1  
#### (4) COMPLEX_FEEDBACK = 4;    Conspecific S_{ii} are between 0.5 - 1 && Heterospecific S_{ij} are between 0.5 - 1.5
#### (5) MIXED_FEEDBACK = 5;      All S_{ij} are between 0.5 - 1.5
########################################################################################################################


#### This function takes in: 
####      (1) FEED_REGIME: an integer between 1 - 5 representing different feedback regimes
####      (2) size: an integer representing the size of the species pool (i.e. number of species)
defineFeed <- function(FEED_REGIME, size = pool.size){

  # Reiterate input variables (i.e. place-holders)
  feed_type <- FEED_REGIME
  
  # If "1", NO_FEEDBACK (i.e. generate a 2-dimension array with all elements = 1)
  if (feed_type == 1){
    feed <- array(1, dim=c(size, size))
  }

  # If "2", POSITIVE_FEEDBACK (i.e. generate a 2-dimension array with off-diagonal elements = 1 & diagonal elements between 1 - 1.5 from a uniform distribution)
  if (feed_type == 2){
    feed <- array(1, dim=c(size, size))
    for (i in 1:size){
      feed[i, i] <- runif(1, 1, 1.5)
    }
  }
  
  # If "3", NEGATIVE_FEEDBACK (i.e. generate a 2-dimension array with off-diagonal elements = 1 & diagonal elements between 0.5 - 1 from a uniform distribution)
  if (feed_type == 3){
    feed <- array(1, dim=c(size, size))
    for (i in 1:size){
      feed[i, i] <- runif(1, 0.5, 1)
    }
  }
  
  # If "4", COMPLEX_FEEDBACK (i.e. generate a 2-dimension array with off-diagonal elements between 0.5 - 1.5 & diagonal elements between 0.5 - 1 from a uniform distribution)
  if (feed_type == 4){
    feed <- array(runif(size^2, 0.5, 1.5), dim=c(size, size))
    for (i in 1:size){
      feed[i, i] <- runif(1, 0.5, 1)
    }
  }
  
  # If "5", MIX_FEEDBACK (i.e. generate a 2-dimension array with all elements between  0.5 - 1.5 from a uniform distribution)
  if (feed_type == 5){
    feed <- array(runif(size^2, 0.5, 1.5), dim=c(size, size))
  }
  
  # Add column and row of NaN
  # This allows indexing "No-species arrival scenario" easier in later parts of the main simulation code (i.e. func_DemoDynAssembly.R)
  feed <- cbind(feed, rep(NA, size))
  feed <- rbind(feed, rep(NA, size+1))
  
  # Output of the function is a list with one element: "feed"
  output <- list(feed = feed)
  return(output)
  
} # End of function body