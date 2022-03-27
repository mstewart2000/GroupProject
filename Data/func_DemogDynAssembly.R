########################################################################################################################
#### MASTER CODE wrapping-up of the entire simulation for the SIMULTANEOUS IMMIGRATION-RERPODUCTION model ##############
#### This code is the one used to run future simulations in a systematic way ###########################################
#### This version of the code randomly generates species pool and environmental quality, instead of allowing input #####
#### At 2017/03/16, the code was updated to allow inputing species pool and environmental quality ######################
#### At 2017/03/19, the code was updated to allow inputing specific demongraphically dynamic feedback matrix ###########
#### At 2018/11/08, the code was updated to allow inputing specific slope and bidirectional matrices ###################
#### For commenting, most 'objects' and 'functions' are indicated with ' ', and "variable" are indicated with " " ###### 
########################################################################################################################


#### This function takes in: 
####      (1) time: the time length of the simulation (default is 100) 
####      (2) pool.size: the number of species in the regional species pool (default is 50) 
####      (3) patch.num: the number of independent patches within the meta-community landscape (default = 1)
####      (4) patch.size: the number of cell/sites in one independent 'patch' (default is 1024) 
####      (5) imm.rate: the probability of a species immigrating from the regional species pool, which is the same for all species (default is 0.0) 
####      (6) num.founders: the number of propagules a species produces when it is immigrating from the regional species pool (default is 10) 
####      (7) wavelet.abiotic: the environment wavelet used for the underlying environmental quality (between 1 (very variable) - 4 (very correlated) spatial structure )
####      (8) feed.regime: the static feedback regime used for S_{ij} (1 = none, 2 = positive, 3 = negative, 4 = complex, 5 = mixed) (see also "defineFeed.R")
####      (9) disp.var: the variance used when simulating dispersal distance (default is 5) 
####      (10) r.prob: the probability of a species having internal seed reproduction (default is 1)
####      (11) fecundity: the number of propagules a species produces during internal seed reproduction (default is 5)
####      (12) death.rate: mortality rate for each species (default is 0.4)
####      (13) prop.full: the proportion of patch.size that is being occupied at the beginning of the simulation (default is 1)
####      (14) dyn.increment: the linear increment for each time/age step for demographically dynamic feedbacks (default is 0.025)
####      (15) dyn.direction: the scenario for demographically dynamic feedbacks(default is 'mag')
####      (16) brown.var: the variance used to generate random walk movements when  dyn.direction is 'brownian' (default is 0.01)
####      (17) pool.input: a TRUE/FALSE value indicating whether specific 'species trait values' are supplied or not (default is FALSE, meaning the function will choose randomly by itself)
####      (18) env.input: a TRUE/FALSE value indicating whether specific 'environmental quality values' are supplied or not (default is FALSE, meaning the function will choose randomly by itself)
####      (19) feedback.input: a TRUE/FALSE value indicating whether specific 'feedback matrix' are supplied or not (default is FALSE, meaning the function will choose randomly by itself)
####      (20) slope.input: a TRUE/FALSE value indicating whether specific 'slope matrix' for the demongraphically dynamic feedback is supplied or not (default is FALSE, meaning the function will choose randomly by itself)
####      (21) bidirection.input: a TRUE/FALSE value indicating whether specific 'bidirectionality matrix' (a matrix of -1 & 1) for the demongraphically dynamic feedback is supplied or not (default is FALSE, meaning the function will choose randomly by itself)
####      (22) pool.list: a vector of specific 'species trait values' (default is NA, meaning the function will choose randomly by itself)
####      (23) env.list: a vector of specific 'environmental quality values' (default is NA, meaning the function will choose randomly by itself)
####      (24) feedback.matrix: a specific 'feedback matrix' (default is NA, meaning the function will choose randomly by itself)
####      (25) slope.matrix: a specific 'slope matrix' used to modified the feedback matrix (default is NA, meaning the function will choose randomly by itself)
####      (26) bidirection.matrix: a specific 'bidirectionality matrix' (matrix of -1 & 1) used to modified the feedback matrix (default is NA, meaning the function will choose randomly by itself)
func_DemogDynAssembly_IRcombine <- function(time = 100, 
                                            pool.size= 50, 
                                            patch.num = 1,
                                            patch.size = 1024, 
                                            imm.rate = 0.01, 
                                            num.founders = 10, 
                                            wavelet.abiotic, 
                                            feed.regime, 
                                            disp.var = 100, 
                                            r.prob = 1,
                                            fecundity = 10,
                                            death.rate = 0.4,
                                            prop.full = 0,
                                            dyn.increment = 0.025,
                                            dyn.direction = 'mag',
                                            brown.var = 0.01,
                                            pool.input = FALSE,
                                            env.input = FALSE,
                                            feedback.input = FALSE, 
                                            slope.input = FALSE, 
                                            bidirection.input = FALSE, 
                                            pool.list = NA, 
                                            env.list = NA,
                                            feedback.matrix = array(NA),
                                            slope.matrix = array(NA),
                                            bidirection.matrix = array(NA)){

  
  ##############################################################
  ##############################################################
  #### PARAMETER/SAVING SPACE SETUP !! #########################
  ##############################################################
  ##############################################################
  
  ##############################################################
  #### Reiterate input variables (i.e. place-holders, all-caps, connected via underscore)
  TIME <- time
  POOL_SIZE <- pool.size
  PATCH_NUM <- patch.num
  PATCH_SIZE <- patch.size
  IMM_RATE <- imm.rate
  NUM_FOUNDERS <- num.founders
  WAVELET_ABIOTIC <- wavelet.abiotic
  FEED_REGIME <- feed.regime
  DISPERSAL_VAR <- disp.var
  REPRODUCTION_PROB <- r.prob
  FECUNDITY <- fecundity
  DEATH_RATE <- death.rate
  PROP_FULL <- prop.full
  DYN_INCREMENT <- dyn.increment
  DYN_DIRECTION <- dyn.direction
  BROWN_VAR <- brown.var
  POOL_INPUT <- pool.input
  ENV_INPUT <- env.input
  FEEDBACK_INPUT <- feedback.input
  SLOPE_INPUT <- slope.input
  BIDIRECTION_INPUT <- bidirection.input
  POOL_LIST <- pool.list 
  ENV_LIST <- env.list
  FEEDBACK_MATRIX <- feedback.matrix
  SLOPE_MATRIX <- slope.matrix
  BIDIRECTION_MATRIX <- bidirection.matrix
  
  
  ##############################################################
  #### Dummy-proof code to ensure dimensions were provided correctly (added for versions after 2017/03/16)
  #### If POOL_INPUT == TRUE, the input vector (i.e. POOL_LIST) needs to have the same length as POOL_SIZE
  #### Check if species trait list data was input correctly, exit function if not
  if((POOL_INPUT == TRUE) & (length(POOL_LIST) != POOL_SIZE)) {
    print("Error in species trait list input! Length of input species trait vector not as claimed"); stop()
  }
  
  #### If ENV_INPUT == TRUE, the input vector (i.e. ENV_LIST) needs to have the same length as PATCH_SIZE
  #### Check if environment quality list data was input correctly, exit function if not
  if((ENV_INPUT == TRUE) & (length(ENV_LIST) != PATCH_SIZE)) {
    print("Error in environment quality list input! Length of input environment quality vector not as claimed"); stop()
  }
  
  #### If FEEDBACK_INPUT == TRUE, the input matrix (i.e. FEEDBACK_MATRIX) needs to have the same dimension as (POOL_SIZE+1)
  #### Check if demongraphically dynamic feedback matrix was input correctly, exit function if not
  #### 'FEEDBACK_MATRIX' must be provided with caution! The dimension of this matrix should be 'POOL_SIZE'+1 since it should include a 'fake species' with NA values for indexing 'no species scenario'
  if((FEEDBACK_INPUT == TRUE) & (dim(FEEDBACK_MATRIX)[1] != (POOL_SIZE + 1))) {
    print("Error in demongraphically dynamic feedback matrix input! Dimension of input feedback matrix not as claimed (e.g. Did you added NA to index 'no species scenario'?)"); stop()
  }  
  
  #### If SLOPE_INPUT == TRUE, the input matrix (i.e. SLOPE_MATRIX) needs to have the same dimension as (POOL_SIZE+1)
  #### Check if slope matrix was input correctly, exit function if not
  #### 'SLOPE_MATRIX' must be provided with caution! The dimension of this matrix should be 'POOL_SIZE'+1 since it should include a 'fake species' 
  #### For simplicity, here you MUST just input a matrix with dimension 'POOL_SIZE'+1 that is full of values, DO NOT add in NA for the 'fake speces' 
  if((SLOPE_INPUT == TRUE) & (dim(SLOPE_MATRIX)[1] != (POOL_SIZE + 1))) {
    print("Error in slope matrix input! Dimension of input feedback matrix not as claimed (e.g. Did you added extra values for 'no species scenario'?)"); stop()
  }
  
  #### If BIDIRECTION_INPUT == TRUE, the input matrix (i.e. BIDIRECTION_MATRIX) needs to have the same dimension as (POOL_SIZE+1)
  #### Check if bidirectionality matrix was input correctly, exit function if not
  #### 'BIDIRECTION_MATRIX' must be provided with caution! The dimension of this matrix should be 'POOL_SIZE'+1 since it should include a 'fake species' 
  #### For simplicity, here you MUST just input a matrix with dimension 'POOL_SIZE'+1 that is full of values (-1 or 1), DO NOT add in NA for the 'fake speces' 
  if((BIDIRECTION_INPUT == TRUE) & (dim(BIDIRECTION_MATRIX)[1] != (POOL_SIZE + 1))) {
    print("Error in bidirectionality matrix input! Dimension of input feedback matrix not as claimed (e.g. Did you added extra values for 'no species scenario'?)"); stop()
  }
  
  
  ##############################################################
  #### Setup the MAIN MATRIX FOR SAVING "SPECIES_ID" AT EACH PATCH THROUGH TIME (THIS IS THE MAIN OUTPUT OF THIS FUNCTION)
  ## Generate 'patch_time': a 2-dimension array with "site_ID" as rows (1st-dim.) and "time_step" as columns (2nd-dim.) (i.e. PATCH_SIZE and TIME, respectively)
  ## Generate 'full_patch_time': a 3-dimension array with "site_ID" as rows (1st-dim.) and "time_step" as columns (2nd-dim.) for each "independent meta-community patches" (3rd-dim.) (i.e. PATCH_SIZE, PATCH_NUM, and TIME, respectively)
  patch_time <- array(dim=c(PATCH_SIZE, TIME))
  full_patch_time <- array(dim=c(PATCH_SIZE, TIME, PATCH_NUM))
  
  
  ##############################################################
  #### Setup the MAIN MATRIX FOR SAVING "INDIVIDUAL_AGE" AT EACH PATCH THROUGH TIME (THIS IS THE NEW STUFF OF THIS MODEL)
  #### These arrays are initiated with zeros, such that if nothing colonzies the site, it will have zero value since it was originated
  ## Generate 'age': a zero vector of plant age/time occupied by an individual for each site (this is temporal snapshot slice)
  ## Generate 'age_array': a zero 3-dimension arrary of plant age for each site (1st-dim.) through time (2nd-dim.) for each independent patch (3rd-dim.)
  age <- rep(0, PATCH_SIZE)
  age_array <- array(0, dim=c(PATCH_SIZE, TIME, PATCH_NUM))
  
  
  ##############################################################
  #### Setup the MAIN MATRIX FOR SAVING FEEDBACK PLAYERS' IDENTITY
  #### This matrix tells you which species generated the feedback (i.e. what was there at the last time step)
  ## Generate 'soil_feed': a 3-dimension array saving the "species_ID" at each site (1st-dim.) over time (2nd-dim.) for each independent patch (3rd-dim.)
  soil_feed <- array(0, dim=c(PATCH_SIZE, TIME, PATCH_NUM))
  
  
  ##############################################################
  #### Saving spaces related to IMMIGRATION
  ## Generate 'num.imm': a vector recording the "number of species" that were allowed to immigrate from the regional species pool
  ## Generate 'founders_time': a list recording the "species_ID" that were allowed to immigrate from the regional species pool
  num.imm <- c()
  founders_time <- list()
  
  
  ##############################################################
  #### Setup the MAIN MATRIX FOR COMPETITION (i.e. c_{ijk} value in Fukami & Nakajima, 2011, Ecology Letters)
  #### For competition at the regional immigration period
  ## Generate 'competeI': a 2-dimension array recording how many multiple landing calculations for immigrants (i.e. competition events) need to be done for this time step (1st-dim.) and patch (2nd-dim.)
  ## Generate 'FEEDBACK_time': a list recording the feedback strength at each time step during the immigration process
  competeI <- array(dim=c(TIME, PATCH_NUM))
  FEEDBACK_time <- list()
  
  
  ##############################################################
  #### Setup the MAIN MATRIX FOR COMPETITION (i.e. c_{ijk} value in Fukami & Nakajima, 2011, Ecology Letters)
  #### For competition at the internal reproduction period
  ## Generate 'competeII': a 2-dimension array recording how many multiple landing calculations for seeds (i.e. competition events) need to be done for this time step (1st-dim.) and patch (2nd-dim.)
  ## Generate 'seed_FEEDBACK_time': a list recording the feedback strength at each time step during the reproduction process
  competeII = array(dim=c(TIME, PATCH_NUM))
  seed_FEEDBACK_time <- list()
  
  
  ##############################################################
  #### Setup FINAL WINNER COMPETITION VALUE saving space
  #### This is recreated in the if-statements for "w/immigrants" and "w/reproducers"
  #### We generate them here because this object would otherwise not exist in scenarios where there were no immigration/reproduction, causing an error when determining winners
  #### Set them as 'numeric(0)' because you don't want it to really have value when there was no immigration/reproduction (i.e. you want it to exist but have empty entry)
  ## Generate 'FINAL_het_compete': a vector for immigrant's winning C_{ijkt} value at each local site (i.e. length is PATCH_SIZE)
  ## Generate: 'FINAL_het_seed_compete': a vector for reproducing seed's winning C_{ijkt} value at each local site (i.e. length is PATCH_SIZE)
  FINAL_het_compete <- numeric(0)
  FINAL_het_seed_compete <- numeric(0)
  
  
  ##############################################################
  #### Setup saving spaces for species identity distribution at different SIMULATION STAGES
  ## Generate 'patch_BeforeImm': a vector of length 'PATCH_SIZE' to save the renewed 'patch' before immigration
  ## Generate 'patch_AfterRecruit': a vector of length 'PATCH_SIZE' to save the renewed 'patch' after immigration/reproduction
  ## Generate 'patch_AfterDeath':  a vector of length 'PATCH_SIZE' to save the renewed 'patch' after mortality
  ## Generate 'empty_sites': a 2-dimension array recording how many sites remained empty after the simulation for each independent patch (1st-dim.) at each time step (2nd-dim.)
  ## Generate 'life_cycle_status': a huge 4-dimension array recording "species_ID change" for each site (1st-dim.), at different stages (2nd-dim.), for each independent patch (3rd-dim.) through time (4th-dim.)
  patch_BeforeImm <- c()
  patch_AfterRecruit <- c()
  patch_AfterDeath <- c()
  empty_sites <- array(dim=c(PATCH_NUM, TIME))
  life_cycle_status <- array(dim=c(PATCH_SIZE, 3, PATCH_NUM, TIME))
  
  
  ##############################################################
  #### Setup of SPECIES POOL (i.e. POOL_SIZE, Z_{i} in Fukami & Nakajima, 2011, Ecology Letters)
  #### After the 2017/03/16 version, depending on 'POOL_INPUT', 'pool' could be randomly generated or specifically supplied (i.e. POOL_LIST)
  ## Generate 'pool': a vector of trait value for each species (from a uniform random between 0 - 1) based on regional pool size 
  if(POOL_INPUT == TRUE) {pool <- POOL_LIST}
  if(POOL_INPUT == FALSE) {pool <- runif(POOL_SIZE)}
  
  
  ##############################################################
  #### Setup ENVIRONMENTAL QUALITY / HETEROGENEITY
  #### After the 2017/03/16 version, depending on 'ENV_INPUT', 'patch_cond' could be randomly generated or specifically supplied (i.e. ENV_LIST)
  ## Generate 'patch_cond': a 2-dimension array with "site_ID" (1st-dim.) and "number of independent patches" (2nd-dim.) (i.e. PATCH_SIZE and PATCH_NUM, respectively)
  patch_cond <- array(dim=c(PATCH_SIZE, PATCH_NUM))
  
  
  #### Write in the environmental quality (i.e. H_{K} in Fukami & Nakajima, 2011, Ecology Letters)
  #### Based on Zee & Fukami, 2015, Ecography, we set the environment quality based on wavelets
  #### The 'wavelets.RData' has four matrix to choose from, each with different degree of spatial auto-correlation (i.e. WAVELET_ABIOTIC, between 1 - 4)
  #### Within a fixed spatial auto-correlation, the matrix is 100 (100 randomly-generated environment) * 1024 (i.e. PATCH_SIZE)
  if(ENV_INPUT == TRUE) {patch_cond <- ENV_LIST}
  if(ENV_INPUT == FALSE) {patch_cond <- WaveletTransformed_H[ceiling(runif(1)*100), 1:PATCH_SIZE, WAVELET_ABIOTIC]}
  # patch_cond <- rep(0.5, PATCH_SIZE)      ### Fix values
  # patch_cond <- runif(PATCH_SIZE)         ### Random environmental values
  # patch_cond <- sort(runif(PATCH_SIZE))   ### Cline environmental values
  
  
  ##############################################################
  #### INITIATE SPECIES DISTRIBUTION at the first time step 
  #### This depends on what proportion of one patch is initially occupied (i.e. PROP_FULL)
  #### Create a vector of species_ID (length = number of occupied sites, the 1st 'sample' function), combine it with zeros, and reshuffle (the 2nd 'sample' function)
  ## Generate 'InitSpp': a vector of initial species distribution with "Species_ID" on each "site"
  InitSpp <- sample(c(sample(POOL_SIZE, round(PATCH_SIZE * PROP_FULL), replace=TRUE), 
                      rep(0, (PATCH_SIZE - round(PATCH_SIZE * PROP_FULL)))))
  
  
  #### Save the initial species distribution to a place-holder (i.e. 'patch' and 'patch_Init')
  ## Generate 'patch': a vector of with the length of "#_of_sites" recording the "species_ID" at each location 
  patch <- InitSpp
  patch_Init <- InitSpp
  
  
  #### Save the initial species distribution to a main saving matrix (i.e. 'patch_time' and 'full_patch_time')
  patch_time[, 1] <- patch
  full_patch_time[, 1, ] <- patch
  soil_feed[, 1, ] <- patch
  
  
  #### FIXED BUG ON Oct. 17, 2018
  #### Save the initial species distribution to 'life_cycle_status'
  #### Honestly, it doesn't really make sense to count Immigration/Recruitment/Death at t=1
  #### We are, however, still doing this to avoid NA in the matrix
  #### For simplicity, let's first set this all equal to the initial species distribution
  life_cycle_status[, , , 1] <- cbind(InitSpp,
                                      InitSpp, 
                                      InitSpp)
  
  
  ##############################################################
  #### UNDERLYING FEEDBACK STRENGTH MATRIX
  #### After the 2017/03/19 version, depending on 'FEEDBACK_INPUT', 'feed' could be randomly generated or specifically supplied (i.e. FEEDBACK_MATRIX)
  #### 'FEEDBACK_MATRIX' must be provided with caution! The dimension of this matrix should be 'POOL_SIZE'+1 since it should include a 'fake species' with NA values for indexing 'no species scenario'
  ## Generate 'feed': a static underlying pairwise feedback matrix, depending on the given feedback scenario (i.e. FEED_REGIME)
  if(FEEDBACK_INPUT == TRUE) {feed <- FEEDBACK_MATRIX}
  if(FEEDBACK_INPUT == FALSE) {feed <- defineFeed(FEED_REGIME, size=POOL_SIZE)$feed}
  
  
  #### After the 2018/11/08 version, depending on 'SLOPE_INPUT' and 'BIDIRECTION_INPUT', matrices associated with the PDF scenario could be randomly generated or specifically supplied
  #### 'SLOPE_MATRIX' and 'BIDIRECTION_MATRIX' must be provided with caution! The dimension of this matrix should be 'POOL_SIZE'+1 since it should include a 'fake species' 
  #### For simplicity, here you MUST just input a matrix with dimension 'POOL_SIZE'+1 that is full of values (-1 or 1), DO NOT add in NA for the 'fake speces' 
  #### Setup demongraphically dynamic feedback depending on the static feedback matrix (i.e. 'feed'), and other parameters specifying the shape (i.e. DYN_INCREMENT, DYN_DIRECTION, BROWN_VAR)
  ## Generate 'demogDyn': a 3-dimension array for pairwise demographically dynamic feedback through time 
  demogDyn <- generateDemogFeedbacks(feed,
                                     depth = (100 * (1 / DEATH_RATE)), 
                                     increment = DYN_INCREMENT, 
                                     direction = DYN_DIRECTION,
                                     slope.array = SLOPE_MATRIX, 
                                     bidirect.array = BIDIRECTION_MATRIX,
                                     brownian = BROWN_VAR)$dyn.feedback
  
  ##############################################################
  ##############################################################
  #### END OF PARAMETER/SAVING SPACE SETUP !! ##################
  ##############################################################
  ##############################################################
  
  
  
  
  
  
  ##############################################################
  ##############################################################
  #### START SIMULATION !! #####################################
  ##############################################################
  ##############################################################
  
  for (t in 2:TIME){
    
    for (j in 1:PATCH_NUM){
      
      
      ###########################################################
      #### IMMIGRATION PROCESS!! ################################
      ###########################################################
      
      #### Preparation before immigration at this time step
      #### At this step, 'patch' is the "species_ID distribution" of the previous time step (i.e. 'patch' = 'full_patch_time[, t-1, j]' = 'patch_BeforeImm')
      ## Generate 'patch_BeforeImm': the "species_ID at each location" before immigration at this time step (you create this even if there is no immigration)
      patch <- full_patch_time[, t-1, j]
      patch_BeforeImm <- full_patch_time[, t-1, j]
    
      
      #### Choose immigrants via comparing random number (length == #_of_species, i.e. POOL_SIZE) to immigration probability (i.e. IMM_RATE)
      ## Generate 'imm_spp': a vector of "species_ID" that were allowed to immigrate from the regional species pool
      imm_spp <- which(runif(POOL_SIZE) < IMM_RATE)
      
      
      #### Record immigrating species' information
      #### That is, "species ID" in 'founders_time[[t]]'; "number of species" in 'num.imm[t]'; "trait value" in 'imm_alpha'
      founders_time[[t]] <- imm_spp
      num.imm[t] <- length(imm_spp)
      imm_alpha <- pool[imm_spp]     
      
      
      #### Renew 'FINAL_het_compete' as numeric(0) (hence with length == 0) before immigration establishment
      #### This is because, in cases where there were no immigrants, the following if-statement would not run, and the code would use 'FINAL_het_compete' from last time step when determining final winner 
      #### In case 'FINAL_het_compete' was not empty last time step, we would have used the wrong 'FINAL_het_compete' to renew 'patch'
      #### This wasn't an issue before since in previous versions, 'patch' was modified before leaving the immigration process (but now we do not update patch until the later outside the immigration for-loop)
      FINAL_het_compete <- numeric(0)
      
      
      #############################
      #### IMMIGRATION ESTABLISHMENT !
      #### If there are any species that were allowed to immigrate from the species pool...
      if (length(imm_spp) > 0){
        
        #### Record the "species_ID" that landed on the site (if any, set as zero otherwise)
        ## Generate 'tempEstab': a 2-dimension array with "patch_ID" (1st-dim.) and "number of immigrating species" (2nd-dim.)
        tempEstab <- array(0, dim=c(PATCH_SIZE, length(imm_spp)))
        
        
        ##############
        #### For each 'immigrating species'...
        for (i in 1:length(imm_spp)){
          
          #### Select sites for each immigrating propagule (i.e. in total NUM_FOUNDERS of sites) to land (i.e. projected landing site)
          ## Generate 'landing': a vector of 'site_ID'
          landing  <- ceiling(runif(NUM_FOUNDERS) * PATCH_SIZE)
          
          
          #### Write "species_ID" (i.e. 'imm_spp[i]') separately into each column of 'tempEstab' (irregardless whether it is currently empty or not)
          tempEstab[landing, i] <- imm_spp[i] 
          
          
          # #### The code below was the source of errors!
          # #### This was never as issue in previous verions because we never asked "which(tempTotEstab != 0)" like we did here when determining the final winner
          # #### When the patch is occasionally completely full, "which(patch[landing] == 0)" and 'tmp.land' is 'numeric(0)'
          # #### The resulting 'tempEstab' would be completely zero and 'mult_landing_array' would be emply, both causing issue when later indexing the final winner
          # #### Since we use 'patch[patch == 0]' in the last code when renewing 'patch', one solution for this error is just to scarifies run speed and populate 'tempEstab' as it would be 
          # ##############################################################
          # #### Remove already occupied sites from the projected landing site (i.e. propagules cannot colonize already occupied sites)
          # ## Generate 'tmp.land': a vector of 'site_ID' that would really allow propagule to land
          # tmp.land <- landing[which(patch[landing] == 0)]
          # 
          # 
          # #### If this species' projected landing site is empty/available, write its "species_ID" (i.e. 'imm_spp[i]') in 'tempEstab'
          # if (length(tmp.land) > 0){
          #   # tempEstab[landing, i] <- imm_spp[i] 
          #   tempEstab[tmp.land, i] <- imm_spp[i] 
          # }
          # ##############################################################
          
        } 
        # End of for-loop for populating 'tempEstab' matrix (i.e. finished determining where immigrating propagules land for each species)
        ##############
        
        
        #### Get the "site_ID" that have MULTIPLE IMMIGRATING PROPAGULES with simultaneous arrivals
        mult_landing <-  which(rowSums(tempEstab > 0) > 1)   
        
        
        #### Record how many multiple landing calculations need to be done for this time step
        competeI[t, j] <- length(mult_landing)
        

        ##############
        #### If there is indeed multiple propagules landing at a single site...
        if (length(mult_landing) > 0){
          
          #### Set up a temporary trait value table with the fake species 
          #### Include a "NA" value for indexing "no species" case (it's "species_ID" being "regional richness + 1") (i.e. a fake species)
          ## Generate 'tempPool': a vetor of "species_trait value" with the fake species
          tempPool <- c(pool, NA)
          
          
          #### Subset the 'tempEstab' array to only consider multiple landing 
          #### This 'drop=FALSE' is to avoid R from automatically downscaling a single row matrix to a vector 
          ## Generate 'mult_landing_array': a 2-dimension array with "patch_ID where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.)
          mult_landing_array <- tempEstab[mult_landing, , drop=FALSE]
          
          
          #### Note that for "sites where multiple landing happened", not all immigrating species would arrive (e.g. only n out of POOL_SIZE might arrive)
          #### If this is the case, modify that element in 'mult_landing_array[site, species]' (which is the "species_ID") to the fake species (i.e. it's "species_ID" being "regional richness + 1")
          mult_landing_array[mult_landing_array == 0] <- (length(pool) + 1)
         
          
          #### Get the environmental condition at those sites where competition happens due to multi-landing
          compete_cond <- patch_cond[mult_landing]
          
          
          #### Expand a matrix to facilitate calculation of competitive effects
          ## Generate 'condition_array': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.), and environment quality as enteries
          condition_array <- matrix(rep(compete_cond, dim(mult_landing_array)[2]), 
                                    length(compete_cond), dim(mult_landing_array)[2])
          
          
          #### DETERMINE FEEDBACKS !
          #### Temporally save the feedback strength values, that happen at those sites for those species, from the underlying feedback matrix
          ## Generate 'FEEDBACK': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.)
          FEEDBACK <- array(dim=dim(mult_landing_array))
          
          
          ##############
          for(i in 1:length(mult_landing)){
            
            #### If that site (i.e. mult_landing[i]) at last time step (i.e. t-1) was not occupied (i.e. soil_feed == 0), THERE IS NO FEEDBACK 
            #### Assign corresponding feedback for that "site_ID" as 1, but leave fake species as NA (i.e. the original setup value)
            if (soil_feed[mult_landing[i], t-1, j] == 0){
              FEEDBACK[i, mult_landing_array[i, ] != (POOL_SIZE + 1)] <- 1
            }
            
            #### If that site (i.e. mult_landing[i]) at last time step (i.e. t-1) was occupied (i.e. soil_feed != 0), THERE IS FEEDBACK!
            #### Assign corresponding feedback for that "site_ID" as: S_{ijt}:
            #### 'i' = which species lands and receives feedback (i.e. mult_landing_array[i, ])
            #### 'j' = which species was originally there and thus imposes feedback (i.e. soil_feed[mult_landing[i], t-1, j], recall that mult_landing is "site_ID")
            #### 't' = how long has the previous species been there (i.e. age_array[mult_landing[i], t-1, j]+1
            #### Push up by age by 1 since plants that create feedback are at least 1 year old by default (i.e. 'demogDyn' can only be indexed starting from [, , i>1])
            #### That is, the + 1 is to allow indexing of 'demogDyn' while remaining meaningful age_array setup (i.e. at age 0) at the same time
            if (soil_feed[mult_landing[i], t-1, j] != 0){
              FEEDBACK[i, ] <- demogDyn[mult_landing_array[i, ], soil_feed[mult_landing[i], t-1, j], (age_array[mult_landing[i], t-1, j] + 1)]    
            }
            
          } # End of for-loop extracting feedback strength from underlying feedback matrix
          ##############

            
          #### Save the 'FEEDBACK' matrix into the 'FEEDBACK_time' list (one list element for each time step)
          FEEDBACK_time[[t]] <- FEEDBACK
      
          
          #### COMPETITION !
          #### Create a matrix of "species_trait" value (it has same values for each column, except in the case where the species didn't arrive at the site, then it is NA, the fake species value)
          #### Note that this matrix has the same dimension as 'condition_array', but 'condition_array' has the same value for each row, not for each column
          ## Generate 'alpha_compete': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.)
          alpha_compete <- matrix(tempPool[mult_landing_array], nrow=length(mult_landing), ncol=length(imm_spp))
        
          
          #### Create a matrix to calculate competition advantage based on: C_{ijk} = (1 - |H_{k} - Z_{i}|) * S_{ij} in Fukami & Nakajima, 2011, Ecology Letters
          #### Here, competition advantage is calculated as: C_{ijkt} = (1 - |H_{k} - Z_{i}|) * S_{ijt}
          ## Generate 'het_compete': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.)
          het_compete <- (1 - abs(condition_array - alpha_compete)) * FEEDBACK
          
          
          #### Create a vector to record the "winner's C_{ijkt} value" for each site (i.e. the maximum value for each row)
          ## Generate 'het_win': a vector with length equal to "#_of_sites where multiple landing happened"
          het_win <- apply(het_compete, 1, max, na.rm=TRUE)
          
          
          #### Know if the species was the winner at that site or not by comparing winner's C_{ijkt} values at each site (i.e. 'het_win') to the realized value for all species/site combination (i.e. 'het_compete')
          ## Generate 'het_win_ID': a 2-dimension array with (dim equal that of 'het_compete')
          het_win_ID <- (het_compete == het_win)
          
          
          #### Determine the "winner's species_trait value" from 'alpha_compete' for each site 
          #### Determine the "winner's species_ID" from 'pool' and 'winning_trait' for each site 
          ## Generate 'winning_trait': a single-column 2-dimension array for each site where multiple landing happened, with "species_trait values" as entries
          ## Generate 'winning_spp': a single-column 2-dimension array for each site where multiple landing happened, with "species_ID" as entries
          winning_trait <- array(dim=c(length(mult_landing), 1))
          winning_spp <- array(dim=c(length(mult_landing), 1))
          
          for (i in 1:length(mult_landing)){
            winning_trait[i, 1] <- alpha_compete[i, which(het_win_ID[i, ] == TRUE)]
            winning_spp[i, 1] <- which(pool == winning_trait[i, ])
          }
          
        } # End of if-statement for multiple propagules landing
        ##############
                
        
        #### Record the "species_ID" of the winning immigrant, but only let it establish after comparing to reproductive winners
        #### If there is only single landing and no multiple landing, then the "species_ID" is just simply the 'rowSums' 
        #### This is because the enteries of 'tempEstab' are 'species_ID' from the very beginning
        #### If there is multiple landing, update that element in 'tempTotEstab' from 'rowSums' to the 'winning_spp'
        ## Generate 'tempTotEstab': a vector of "species_ID" that established after successful immigration
        tempTotEstab <- rowSums(tempEstab)
        if (length(mult_landing) > 0){
          tempTotEstab[mult_landing] <- winning_spp
        }
        
        
        # ##############################################################
        # #### Deactivated in this 'simultaneous immigration-reproduction version' (i.e. do not allow winners at the reproduction stage immediately colonize 'patch', do so only after comparing with immigrant winners)
        # #### Save 'species_ID distribution' (i.e. patch) after immigration establishment
        # #### Only those that were originally empty can be established (or remained empty, depending in there was any landing at all)
        # patch[patch == 0] <- tempTotEstab[patch == 0]
        # ##############################################################
        
        
        #### Record the winning immigrant's trait value
        ## Generate 'tempTotEstab_trait': a vector of "species_ID" that established after successful immigration
        tempTotEstab_trait <- array(dim=length(tempTotEstab))
        tempTotEstab_trait[tempTotEstab != 0] <- pool[tempTotEstab[tempTotEstab != 0]]
        
        
        #### Recalculate the C_{ijkt} values for the immigrant winners at all local sites
        #### We have to do this because is the previous code, we only calculated C_{ijkt} for sites with multiple landing
        #### When there are no species, skip the calculation and set its C_{ijkt} as 0
        ## Generate 'FINAL_FEEDBACK': a vector saving feedback strength (i.e. S_{ijt}) for "sites where there is a winning immigrant
        ## Generate 'FINAL_het_compete': a vector saving C_{ijkt} for all local sites
        ## Generate 'focal_cell': a 'site_ID' for non-zero immigrants that need its C_{ijkt} calculated
        FINAL_FEEDBACK <- array(dim=length(which(tempTotEstab != 0)))
        FINAL_het_compete <- rep(0, length(tempTotEstab))
        
        for(i in 1:length(which(tempTotEstab != 0))){
          focal_cell <- which(tempTotEstab != 0)[i]
          if (soil_feed[focal_cell, t-1, j] == 0){
            FINAL_FEEDBACK[i] <- 1
          }
          
          if(soil_feed[focal_cell, t-1, j] != 0){
            FINAL_FEEDBACK[i] <- demogDyn[tempTotEstab[focal_cell], soil_feed[focal_cell, t-1, j], (age_array[focal_cell, t-1, j] + 1)]    
          }
        }
        
        FINAL_het_compete[which(tempTotEstab != 0)] <- (1 - abs(patch_cond[which(tempTotEstab != 0)] - tempTotEstab_trait[which(tempTotEstab != 0)])) * FINAL_FEEDBACK    
        
      } # End of if-statement for any immigration 
      #############################

      
      # ##############################################################
      # #### Deactivated in this 'simultaneous immigration-reproduction version' (i.e. there is no such thing as 'patch_AfterImm', there is only 'patch_AfterRecruit')
      # #### Save the renewed 'patch' (i.e. after immigration) to 'patch_AfterImm' 
      # #### At this step, 'patch' is the "species_ID distribution" after immigration (i.e. 'patch' = 'patch_AfterImm') (you save this even if there is no immigration)
      # patch_AfterImm <- patch
      # ##############################################################
      
      ###########################################################
      #### END OF IMMIGRATION PROCESS !! ########################
      #### Skip the immigration process if no immigrants ########
      ###########################################################
      
      
      
      
      ###########################################################
      #### REPRODUCTION AND DISPERSAL PROCESS !! ################
      ###########################################################
      
      #### Note: (1) The unit for reproduction/dispersal is individuals (could be different individuals of the same species at different sites)
      #### Note: (2) The unit for immigration is species within the species pool
      
      
      #### Select individuals to reproduction/dispersal, this only happens at site which remained occupied after the immigration step
      ## Generate 'occ_patch': a vector recording which "sites_ID" were occupied, and thus potentially could reproduce
      occ_patch <- which(patch > 0) 
      
      
      #### Choose sites where the occupied individuals get to reproduce via comparing this random number to reproduction probability (i.e. REPRODUCTION_PROB)
      ## Generate 'repDRAW': a vector of random numbers (length == #_of_occupied_site, i.e. occ_patch) from unifor distribution 
      ## Generate 'rep_location' a vector recording the "site_ID" that were chosen to reproduce/disperse
      ## Generate 'reproducers': a vector recording the "species_ID" that were chosen to reproduce/disperse
      repDRAW <- runif(length(occ_patch))      
      rep_locations <- which(patch!=0)[repDRAW < REPRODUCTION_PROB]
      reproducers <- patch[patch != 0][repDRAW < REPRODUCTION_PROB]
      
    
      #### Renew 'FINAL_het_seed_compete' as numeric(0) (hence with length == 0) before reproduction establishment
      #### This is because, in cases where there were no reproducers, the following if-statement would not run, and the code would use 'FINAL_het_seed_compete' from last time step when determining final winner 
      #### In case 'FINAL_het_seed_compete' was not empty last time step, we would have used the wrong 'FINAL_het_seed_compete' to renew 'patch'
      #### This wasn't an issue before since in previous versions, 'patch' was modified before leaving the reproduction process (but now we do not update patch until the later outside the reproduction for-loop)
      FINAL_het_seed_compete <- numeric(0)
      
      
      #### If there are any individual that were allowed to reproduce from existing population..
      if(length(reproducers) != 0){
        
        #### Simulate "disperal distance" for each seed ("#_of_seeds" is FECUNDITY) of each reproducing individual (i.e. 'reproducers') from a normal distribution with given variance (i.e. DISPERSAL_VAR)
        ## Generate 'disp_distances': a matrix of integer "dispersal distance" for each reproducing individual (1st-dim.) and each of its propagules/seeds (2nd-dim.)
        disp_distances <- matrix(round(rnorm(length(reproducers) * FECUNDITY, mean=0, DISPERSAL_VAR)), 
                                 nrow=length(reproducers), ncol=FECUNDITY)
        
        
        #### Determine landing "site_ID" of dispersing reproduction propagules (i.e. current "site_ID" plus dispersing distance)
        ## Generate 'dispersed_locations': a matrix for landing location for each reproducing individual (1st-dim.) and each of its propagules/seeds (2nd-dim.)
        dispersed_locations <- array(dim=c(length(reproducers), FECUNDITY))
        for (i in 1:length(reproducers)){
          dispersed_locations[i, ] <- disp_distances[i, ] + rep_locations[i]
        }
        
        
        #### Adjust for boundary conditions since it is a circular array
        #### If 'dispersed_locations' < 1, it continues its counting from c(PATCH_SIZE, PATCH_SIZE-1, PATCH_SIZE-2, etc.)
        #### If 'dispersed_locations' > PATCH_SIZE, it continues its counting from c(1, 2, 3, etc.)
        dispersed_locations[dispersed_locations < 1] <- PATCH_SIZE - abs(dispersed_locations[dispersed_locations < 1])
        dispersed_locations[dispersed_locations > PATCH_SIZE] <- dispersed_locations[dispersed_locations > PATCH_SIZE] - PATCH_SIZE
    
               
        #### Record the "species_ID" that landed on the site as reproduction propagules (if any, set as zero otherwise)
        #### Note the dimension of 'seed' is "site_ID" * "each reproducing individual"
        #### Note the dimension of 'dispersed_distance' is "each reproducing individual" * "fecundity"
        ## Generate 'tempEstab': a 2-dimension array with "patch_ID" (1st-dim.) and "number of reproducing individuals" (2nd-dim.) 
        ## This resembles 'tempEstab' for immigration 
        seeds <- array(0, dim=c(length(patch), length(reproducers)))
        
        #### For each reproducing individual, place its propagules (#_of_seeds == FECUNDITY)
        #### Note that one individual may place its propagules at the same "site_ID", if so, colSums(seeds>0) would be less than FECUNDITY
        for (i in 1:length(reproducers)){
          seeds[dispersed_locations[i, ], i] <- reproducers[i]
        }
    
             
        #### Get the "site_ID" that have MULTIPLE REPRODUCTION PROPAGULES with simultaneous arrivals
        ## Generate 'mult_seed_location': a vector of "site_ID" where multiple landing occuried
        ## This resembles 'mult_landing' for immigration 
        mult_seed_location <- which(rowSums(seeds > 0) > 1)
        
        #### Record how many multiple landing calculations need to be done for this time step
        ## This resembles 'competeI' for immigration 
        competeII[t, j] <- length(mult_seed_location)
        
        
        ##############
        #### If there is indeed multiple propagules landing at a single site...
        if (length(mult_seed_location) > 0){
          
          #### Subset the 'seeds' array to only consider multiple landing
          #### This 'drop=FALSE' is to avoid R from automatically downscaling a single row matrix to a vector
          ## Generate 'mult_seed_array': a 2-dimension array with "patch_ID where multiple landing happened" (1st-dim.) and "number of reproducing individuals" (2nd-dim.)
          ## This resembles 'mult_landing_array' for immigration
          mult_seed_array <- seeds[mult_seed_location, , drop=FALSE]
          
          
          #### Set up a temporary trait value table with the fake species
          #### Include a "NA" value for indexing "no species" case (it's "species_ID" being "regional richness + 1") (i.e. a fake species)
          ## Generate 'tempPool': a vetor of "species_trait value" with the fake species
          ## This resembles 'tempPool' for immigration
          tempPool <- c(pool, NA)
          
          
          #### Note that for "sites where multiple landing happened", not all reproducing individuals would arrive (e.g. only n out of reproducers might arrive)
          #### If this is the case, modify that element in 'mult_seed_array[site, individual]' (which is the "species_ID") to the fake species (i.e. it's "species_ID" being "regional richness + 1")
          mult_seed_array[mult_seed_array == 0] <- (length(pool) + 1)
          
          
          #### Get the environmental condition at those sites where competition happens due to multi-landing
          compete_cond <- patch_cond[mult_seed_location]
          
          
          #### Expand a matrix to facilitate calculation of competitive effects
          ## Generate 'condition_array': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.), and environment quality as enteries
          condition_array <- matrix(rep(compete_cond, dim(mult_seed_array)[2]), 
                                    length(compete_cond), dim(mult_seed_array)[2])  
          
          
          #### DETERMINE FEEDBACKS !
          #### Temporally save the feedback strength values, that happen at those sites for those species, from the underlying feedback matrix
          #### Note: Individual who arrived only after immigration immediately become mature and can reproduce, however, they do not immediately cause feedback in this code
          ## Generate 'seed_FEEDBACK': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of immigrating species" (2nd-dim.)
          ## This resembles 'FEEDBACK' for immigration
          seed_FEEDBACK <- array(dim=dim(mult_seed_array))
          
          
          ##############
          for(i in 1:length(mult_seed_location)){
            
            #### If that site (i.e. mult_seed_location[i]) at last time step (i.e. t-1) was not occupied (i.e. soil_feed == 0), THERE IS NO FEEDBACK
            #### Assign corresponding feedback for that "site_ID" as 1, but leave fake species as NA (i.e. the original setup value)
            if (soil_feed[mult_seed_location[i], t-1, j] == 0){
              seed_FEEDBACK[i, mult_seed_array[i, ] != (POOL_SIZE + 1)] <- 1
            }
            
            #### If that site (i.e. mult_seed_location[i]) at last time step (i.e. t-1) was occupied (i.e. soil_feed != 0), THERE IS FEEDBACK!
            #### Assign corresponding feedback for that "site_ID" as: S_{ijt}:
            #### 'i' = which species lands and receives feedback (i.e. mult_seed_array[i, ])
            #### 'j' = which species was originally there and thus imposes feedback (i.e. soil_feed[mult_seed_location[i], t-1, j], recall that mult_seed_location is "site_ID")
            #### 't' = how long has the previous species been there (i.e. age_array[mult_seed_location[i], t-1, j]+1
            #### Push up by age by 1 since plants that create feedback are at least 1 year old by default (i.e. 'demogDyn' can only be indexed starting from [, , i>1])
            #### That is, the + 1 is to allow indexing of 'demogDyn' while remaining meaningful age_array setup (i.e. at age 0) at the same time
            if (soil_feed[mult_seed_location[i], t-1, j] != 0){
              seed_FEEDBACK[i, ] <- demogDyn[mult_seed_array[i, ], soil_feed[mult_seed_location[i], t-1, j], (age_array[mult_seed_location[i], t-1, j] + 1)]
            }
            
          } # End of for-loop extracting feedback strength from underlying feedback matrix
          ##############

                    
          #### Save the 'seed_FEEDBACK' matrix into the 'seed_FEEDBACK_time' list (one list element for each time step)
          ## This resembles 'FEEDBACK_time' for immigration 
          seed_FEEDBACK_time[[t]] <- seed_FEEDBACK    
          
          
          #### COMPETITION !
          #### Create a matrix of "species_trait" value (it has same values for each column, except in the case where the species didn't arrive at the site, then it is NA, the fake species value)
          #### Note that this matrix has the same dimension as 'condition_array', but 'condition_array' has the same value for each row, not for each column
          ## Generate 'alpha_compete': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of reproducing individuals" (2nd-dim.)
          ## This resembles 'alpha_compete' for immigration
          alpha_compete <- matrix(tempPool[mult_seed_array], nrow=length(mult_seed_location), ncol=length(reproducers))
          
          
          #### Create a matrix to calculate competition advantage based on: C_{ijk} = (1 - |H_{k} - Z_{i}|) * S_{ij} in Fukami & Nakajima, 2011, Ecology Letters
          #### Here, competition advantage is calculated as: C_{ijkt} = (1 - |H_{k} - Z_{i}|) * S_{ijt}
          ## Generate 'het_seed_compete': a matrix with "sites where multiple landing happened" (1st-dim.) and "number of reproducing individuals" (2nd-dim.)
          ## This resembles 'het_compete' for immigration
          het_seed_compete <- (1 - abs(condition_array - alpha_compete)) * seed_FEEDBACK
          
          
          #### Create a vector to record the "winner's C_{ijkt} value" for each site (i.e. the maximum value for each row)
          ## Generate 'het_seed_win': a vector with length equal to "#_of_sites where multiple landing happened"
          ## This resembles 'het_win' for immigration
          het_seed_win <- apply(het_seed_compete, 1, max, na.rm=TRUE)
          
          
          #### Know if the individual was the winner at that site or not by comparing winner's C_{ijkt} values at each site (i.e. 'het_seed_win') to the realized value (i.e. 'het_seed_compete')
          ## Generate 'het_seed_win_ID': a 2-dimension array with (dim equal that of 'het_seed_compete')
          ## This resembles 'het_win_ID' for immigration
          het_seed_win_ID <- (het_seed_compete == het_seed_win)
          
          
          #### Determine the "winner's species_trait value" for that individual from 'alpha_compete' for each site 
          #### Determine the "winner's species_ID" for that individual from 
          ## Generate 'winning_seed_trait': a single-column 2-dimension array for each site where multiple landing happened, with "species_trait values" as entries
          ## Generate 'winning_seed': a single-column 2-dimension array for each site where multiple landing happened, with "species_ID" as entries
          ## These resemble 'winning_trait' and 'winning_spp' for immigration
          winning_seed_trait <- array(dim=c(length(mult_seed_location), 1))
          winning_seed <- array(dim=c(length(mult_seed_location), 1))
          
          for (i in 1:length(mult_seed_location)){
            
            #### For reproduction, multiple seeds from different individuals of the same species can arrive, and be determined as winner simultaneously
            ## Generate 'winner_ind': a vector to know the "individual_ID" of the winner(s)
            ## Generate 'winner_ind_want': one number representing THE winning individual (if multiple seeds of the same species arrived, randomly pick one, although they are all the same)
            winner_ind <- which(het_seed_win_ID[i, ] == TRUE)
            winner_ind_want <- ifelse(length(winner_ind) == 1, winner_ind, sample(winner_ind, 1))
            winning_seed_trait[i, 1] = alpha_compete[i, winner_ind_want]
            winning_seed[i, 1] <- mult_seed_array[i, winner_ind_want]
            
          }
  
        } # End of if-statement for multiple propagules landing
        ##############
        
        
        #### Record the "species_ID" of the successful reproducing species, but only let it establish after comparing to immigrant winners
        #### If there is only single landing and no multiple landing, then the "species_ID" is just simply the 'rowSums' 
        #### This is because the enteries of 'seeds' are 'species_ID' from the very beginning
        #### If there is multiple landing, update that element in 'tempTotEstab' from 'rowSums' to the 'winning_spp'
        ## Generate 'tempTotSeed': a vector of "species_ID" that established after successful immigration
        ## This resembles 'tempTotEstab' for immigration
        tempTotSeed <- rowSums(seeds)
        if (length(mult_seed_location) > 0){
          tempTotSeed[mult_seed_location] <- winning_seed
        }
        
        
        # ##############################################################
        # #### Deactivated in this 'simultaneous immigration-reproduction version' (i.e. do not allow winners at the reproduction stage immediately colonize 'patch', do so only after comparing with immigrant winners)
        # #### Save 'species_ID distribution' (i.e. patch) after reproduction/dispersal seed establishment
        # #### Only those that were originally empty can be established (or remained empty, depending in there was any landing at all)
        # patch[patch == 0] <- tempTotSeed[patch == 0]
        # ##############################################################
        
        
        #### Record the winning reproductive seed's trait value for each local site (leave as NA if no seed landed)
        ## Generate 'tempTotEstab_trait': a vector of "species_trait value" that established after successful reproduction/dispersal
        tempTotSeed_trait <- array(dim=length(tempTotSeed))
        tempTotSeed_trait[tempTotSeed != 0] <- pool[tempTotSeed[tempTotSeed != 0]]
        
        
        #### Recalculate the C_{ijkt} values for the reproductive seed winners at all local sites
        #### We have to do this because is the previous code, we only calculated C_{ijkt} for sites with multiple landing
        #### When there are no species, skip the calculation and set its C_{ijkt} as 0
        ## Generate 'FINAL_FEEDBACK': a vector saving feedback strength (i.e. S_{ijt}) for "sites where there is a winning immigrant
        ## Generate 'FINAL_het_compete': a vector saving C_{ijkt} for all local sites
        ## Generate 'focal_seed_cell': a 'site_ID' for non-zero immigrants that need its C_{ijkt} calculated
        FINAL_het_seed_compete <- rep(0, length(tempTotSeed_trait))
        FINAL_seed_FEEDBACK <- array(dim=length(which(tempTotSeed != 0)))
        
        for(i in 1:length(which(tempTotSeed != 0))){
          focal_seed_cell <- which(tempTotSeed != 0)[i]
          if (soil_feed[focal_seed_cell, t-1, j] == 0){
            FINAL_seed_FEEDBACK[i] <- 1
          }
          
          if(soil_feed[focal_seed_cell, t-1, j] != 0){
            FINAL_seed_FEEDBACK[i] <- demogDyn[tempTotSeed[focal_seed_cell], soil_feed[focal_seed_cell, t-1, j], (age_array[focal_seed_cell, t-1, j] + 1)]    
          }
        }
        
        FINAL_het_seed_compete[which(tempTotSeed != 0)] <- (1 - abs(patch_cond[which(tempTotSeed != 0)] - tempTotSeed_trait[which(tempTotSeed != 0)])) * FINAL_seed_FEEDBACK    
            
      } # End of if-statement for any reproduction
      #############################
      
      
      # ##############################################################
      # #### Deactivated in this 'simultaneous immigration-reproduction version' (i.e. there is no such thing as 'patch_AfterRep', there is only 'patch_AfterRecruit')
      # #### Save the renewed 'patch' (i.e. after reproduction/dispersal) to 'patch_AfterRep'
      # #### At this step, 'patch' is the "species_ID distribution" after reproduction (i.e. 'patch' = 'patch_AfterRep')
      # patch_AfterRep <- patch
      # ##############################################################
      
      ###########################################################
      #### END OF REPRODUCTION AND DISPERSAL PROCESS !! #########
      ###########################################################
      
      
      
      
      ###########################################################
      #### DETERMINE FINAL WINNER !! ############################
      ###########################################################
      
      #### Previous immigration/reproduction code determines the winner's ID (i.e. tempTotEstab & tempTotSeed) and C_{ijkt} for each process (i.e. FINAL_het_compete & FINAL_het_seed_compete)
      #### This part of the code compares these two winners and select an ultimate winner 
      ## Generate 'FINAL_Recruit': a vector with size "PATCH_SIZE" that saves the winner's "species_ID" for each local site
      FINAL_Recruit <- array(dim=PATCH_SIZE)
      
      
      #### When no immigration/reproduction, no extra recruitment, 'FINAL_Recruit' == 'patch'(unchanged) == 'patch_BeforeImm' == 'full_patch_time[, t-1, ]'
      if(length(imm_spp) == 0 & length(reproducers) == 0){
        FINAL_Recruit <- patch
      }
      
      #### When no immigration, 'FINAL_Recruit' is whatever won the reproduction process (i.e. 'tempTotSeed')
      if(length(imm_spp) == 0 & length(reproducers) != 0){
        FINAL_Recruit <- tempTotSeed
      }
      
      #### When no reproduction, 'FINAL_Recruit' is whatever won the immigration process (i.e. tempTotEstab)
      if(length(imm_spp) != 0 & length(reproducers) == 0){
        FINAL_Recruit <- tempTotEstab
      }
      
      #### The above series of if-statements is necessary since if immigration (repreduction) did not happen, tempTotEstab (tempTotSeed) would not exist, and the following code would not apply 
      #### If both immigration/reproduction occurred, 'FINAL_Recruit' is whatever is the better competitor of the two 
      if(length(imm_spp) != 0 & length(reproducers) != 0){
        for(i in 1:PATCH_SIZE){
          FINAL_Recruit[i] <- ifelse((FINAL_het_seed_compete[i] - FINAL_het_compete[i]) > 0, tempTotSeed[i], tempTotEstab[i]) 
        }
      }
      
      
      #### Save the renewed 'patch' (i.e. after reproduction/dispersal) to 'patch_AfterRecruit'
      # #### At this step, 'patch' is the "species_ID distribution" after immigration/reproduction (i.e. 'patch' = 'patch_AfterRecruit')
      patch[patch == 0] <- FINAL_Recruit[patch == 0]
      patch_AfterRecruit <- patch
      
      ###########################################################
      #### END OFDETERMINE FINAL WINNER !! ######################
      ###########################################################
      
      
      
      
      ###########################################################
      #### UPDATING AGE !! ######################################
      ###########################################################
      
      #### This part of the code updates individual age information (i.e. age_array)
      #### Note that this is different from updating 'soil_feed'
      #### For this version of the code, 'age_array' is updated after immigration (i.e. before reproduction), and 'soil_feed' is updated after reproduction, but both are done before death
      
      #### Sites are categorized into: (1) New colonizers of originally empty sites, (2) remain occupying survivors, and (3) remain empty sites
      #### In general, it's 'soil_feed' that is critical is determining whether there is feedback or not, the numeric value of 'age_array' is not that critical 
      #### This is because even if 'age_array' == 0, there is still feedback (i.e. LV-1 feedback) since when indexing we use ('age_array' + 1)
      #### When the value of 'age_array' == 1, then this feedback would be LV-2 feedback
      
      #### Here, I compare 'patch', instead of 'patch_AfterImm' as that in Peter Zee's original code, to 'full_patch_time[, t-1, j]'
      #### The reason is that by doing so, I can move this age updating part freely, and it would determine plant age after the desired life cycle step 
      #### If right before immigration, 'patch' == 'patch_BeforeImm': age affects feedback at both immigration and reproduction processes of this time step
      #### If right after immigration/reproduction, 'patch' == 'patch_AfterRecruit': age affects feedback only at the next time step
      #### If right after death, 'patch' == 'patch_AfterDeath': age affects feedback only at the next time step

      #### Note: Imagine a case where an empty site was colonized during immigration, reproduced seed, but then wiped out due to death events
      #### Note: First thing to think about is whether this individual would create feedback or not (I THINK IT SHOULD)
      #### Note: If so, this individual didn't survived more than one year and thus should create LV-1 feedback (correspondingly, this code correctly counts it as 'empty' and 'age_array' will remain zero)
      
      
      #### Count new colonizers
      #### (1) If 'patch' minus 'full_patch_time[, t-1, j]' equals 'patch_AfterRecruit', that means 'full_patch_time[, t-1, j]' was zero 
      #### (2) Within these, some have zeros for both 'patch' and 'full_patch_time[, t-1, j]', those don't count
      new_state <- which((patch - full_patch_time[, t-1, j]) == patch)
      new_arrivals <- setdiff(new_state, which(patch == 0))
      
      
      #### Count survivors
      #### (1) If 'patch' minus 'full_patch_time[, t-1, j]' equals zero, that means the patch was occupied before and survived (i.e. cannot be colonized by immigration)
      #### (2) Within these, some have zeros for both 'patch' and 'full_patch_time[, t-1, j]', those don't count
      same_state <- which(patch == full_patch_time[, t-1, j])
      survivors <- setdiff(same_state, which(patch == 0))
      
      
      #### Count empty space
      #### Zero at 'full_patch_time[, t-1, j]' AND zero at both 'patch' (the former would have more zeros since it was before immigration)
      empty <- intersect(which((full_patch_time[, t-1, j] == 0)), which((patch == 0)))
      
      
      #### Update age_array
      #### FIXED BUG ON Oct. 17, 2018
      #### Recall that 'age_array' is a zero 3-dimension arrary of plant age for each site (1st-dim.) through time (2nd-dim.) for each independent patch (3rd-dim.)
      #### For the positions of 'survivors', ADD ONE (from 'age_array[, t-1, j]') to the positions and save it to 'age_array[, t, j]'
      #### For the positions of 'new arrivals', SET AS ONE and save it to 'age_array[, t, j]'
      #### No need to change 'empty' to zero since 'age_array' was setup as zero from the beginning
      age_array[new_arrivals, t, j] <- 1
      age_array[survivors, t, j] <- age_array[survivors, t-1, j] + 1
      
      ###########################################################
      #### END OF UPDATING AGE !! ###############################
      ###########################################################
      
      
      
      
      ###########################################################
      #### UPDATING FEEDBACK-GENERATING MEMBERSHIP !! ###########
      ###########################################################
      
      #### For this version of the code, 'soil_feed' is updated after reproduction but before death
      #### Here, I use 'patch', instead of 'patch_AfterRep' as that in Peter Zee's original code, so that I can move this part freely to update 'soil_feed' after the desired life cycle step 
      #### The timing of this line determines what stage does the indiviual need to pass in order to generate feedback
      #### If right before immigration, 'patch' == 'patch_BeforeImm': only FULL SURVIVORS of last time step can generate feedbacks
      #### If right after immigration/reproduction, 'patch' == 'patch_AfterRecruit': successful established individuals after immigration and reproduction can cause feedback, even if it was wipe out later in the time step
      #### If right after death, 'patch' == 'patch_AfterDeath': only FULL SURVIVORS of this time step can generate feedbacks (effectively the same as 'patch_BeforeImm')
    
      
      #### Renew feedback-generating membership after reproduction
      #### At this step, 'patch' is the "species_ID distribution" after immigration/reproduction (i.e. 'patch' = 'patch_AfterRecruit')
      soil_feed[, t, j] <- patch
      
      ###########################################################
      #### END OF UPDATING FEEDBACK-GENERATING MEMBERSHIP !! ####
      ###########################################################
      
      
  
      
      ###########################################################
      #### MORTALITY PROCESS !! #################################
      ###########################################################
      
      #### Choose sites where we kill occupying individuals (set its "species_ID" within patch to zero) via comparing this random number to mortality probability (i.e. DEATH_RATE)
      ## Generate 'deathDRAW': a vector of random numbers (length == #_of_site) from unifor distribution 
      deathDRAW <- runif(PATCH_SIZE)
      patch[which(deathDRAW < DEATH_RATE)] <- 0
      
      #### Save the renewed 'patch' (i.e. after mortality) to 'patch_AfterDeath'
      #### At this step, 'patch' is the "species_ID distribution" after death events (i.e. 'patch' = 'patch_AfterDeath')
      patch_AfterDeath <- patch

      ###########################################################
      #### END OF MORTALITY PROCESS !! ##########################
      ###########################################################
      
      
      
      
      ###########################################################
      #### CALCULATING SIMULATION INDEX/OUTPUT !! ###############
      ###########################################################
      
      #### That was the last step of simulation for one time step
      #### Save the final 'patch' (after immigration, reproduction, and death) to 'full_patch_time[, t, j]'
      #### In addition, calculate how many sites remained empty after the whole simulation after at this time step
      full_patch_time[, t, j] <- patch
      empty_sites[j, t] <- sum(full_patch_time[, t, j] == 0)
      
      
      #### Save "species_ID distribution" (i.e. patch) at different stages into a huge 4-dimension array
      #### This records "species_ID change" for each site (1st-dim.), at different stages (2nd-dim.), for each independent patch (3rd-dim.) through time (4th-dim.)
      #### There are in total four stages: Before immigration, after immigration/reproduction, and after mortality
      life_cycle_status[, , j, t] <- cbind(patch_BeforeImm,
                                           patch_AfterRecruit, 
                                           patch_AfterDeath)
      
 
      ###########################################################
      #### END OF CALCULATING SIMULATION INDEX/OUTPUT !! ########
      ###########################################################
      
      
    } # End of PATCH_NUM for loop (j) 
    
    # #### Processing check
    # if(t %% 10 == 0) {print(t); print(soil_feed[1:20, 1:t, j])}
    
  } # End of big TIME for loop (t)
  
  ##############################################################
  ##############################################################
  #### END OF SIMULATION !! ####################################
  ##############################################################
  ##############################################################
  
  
  
  
  
  
  ##############################################################
  ##############################################################
  #### OUTPUT SIMULATION RESULTS !! ############################
  ##############################################################
  ##############################################################
  
  #### Prepare the function output: $param_list (i.e. the callback list of parameter values)
  param_list <- list(time = time, 
                     pool.size = pool.size, 
                     patch.size = patch.size, 
                     imm.rate = imm.rate, 
                     num.founders = num.founders, 
                     wavelet.abiotic = wavelet.abiotic, 
                     feed.regime = feed.regime, 
                     disp.var = disp.var, 
                     r.prob = r.prob,
                     fecundity =fecundity,
                     death.rate = death.rate,
                     prop.full = prop.full,
                     dyn.increment = dyn.increment,
                     dyn.direction = dyn.direction, 
                     brown.var = brown.var)
  
  
  #### Prepare the function output: $input 
  #### This includes: $input$pool: the regional species pool 'trait value', either supplied through input (i.e. POOL_LIST) or randomly generated
  #### This includes: $input$env: the 'environmental quality value', either supplied through input (i.e. ENV_LIST) or randomly generated
  input <- list(pool = pool, 
                env = patch_cond, 
                feed = feed)
  
  #### Prepare the function output: $output 
  #### This includes: $output$full_patch_time: the "species_ID distribution" among sites through time
  #### This includes: $output$age_array: the "individual age distribution" among sites through time
  #### This includes: $output$soil_feed: the "species_ID distribution" among sites through time recording the player that generates feedback
  #### This includes: $output$life_cycle_status: colonization and extinction events throughout the simulation
  output <- list(full_patch_time = full_patch_time,
                 age_array = age_array,
                 soil_feed = soil_feed,
                 life_cycle_status = life_cycle_status)
  
  
  #### Send out the final output, including $param_list and $output
  results <- list(param_list = param_list, input = input, output = output)
  return(results)
  
  ##############################################################
  ##############################################################
  #### END OF OUTPUT SIMULATION RESULTS !! #####################
  ##############################################################
  ##############################################################
  
} # End of function body




