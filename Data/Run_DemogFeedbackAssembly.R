########################################################################################################################
#### MASTER CODE running a demo simulation where you can set the various parameters ####################################
#### This also will plot the richness through time and a colored map of the community species ID-location through time #
########################################################################################################################


###########################################################
#### Set working directory, library, and link to other MASTER CODE
load('wavelets.RData')
source('generateDynFeedbacks.R')
source('defineFeed.R')
source('func_DemogDynAssembly.R')
library(vegan)
library(ggplot2)


###########################################################
#### Create input data 
pool <- runif(50)
patch_cond <- WaveletTransformed_H[ceiling(runif(1)*100), , 3]
feed <- defineFeed(4, 50)$feed
slope.array <- matrix(runif(nrow(feed) * ncol(feed)), nrow(feed), ncol(feed))
bidirect.array <- matrix(sample(c(-1, 1), replace=TRUE, size=nrow(feed) * ncol(feed)), nrow(feed), ncol(feed))


#### Start the demo simulation using the 'func_DemogDynAssembly_IRcombine' function in the "func_DemogDynAssembly_IRcombine.R" file
#### Check the run time for this demo simulation using "Sys.time", and the file size using "object.size"
tic = Sys.time()
Sim <- func_DemogDynAssembly_IRcombine (time = 1500,                    # number of timesteps
                                        pool.size = 50,                 # number of species in pool
                                        patch.num = 1,                  # number of patches within the meta-community landscape
                                        patch.size = 1024,              # size of patch / number of sites
                                        imm.rate = 0.05,                # immigration rate
                                        num.founders = 10,              # number of founders
                                        wavelet.abiotic = 3,            # wavelet ((1) less to (4) more spatial structure )
                                        feed.regime = 4,                # feedback regime (1=none, 2=positive, 3=negative, 4=complex, 5=mixed)
                                        disp.var = 150,                 # dispersal distance
                                        r.prob = 1,                     # probability of producing seeds
                                        fecundity = 5,                  # number of seeds produced
                                        death.rate = 0.3,               # death rate
                                        prop.full = 0,                  # proportion of patch initially colonized
                                        dyn.increment = 0.02,           # increment of feedback change with age
                                        dyn.direction = 'mag',          # direction of demographic change ('mag', 'decay', 'brownian')
                                        brown.var = 0.01,               # var for brownian motion
                                        pool.input = TRUE,              # TRUE/FALSE value indicating whether specific 'species trait values' are supplied or not
                                        env.input = TRUE,               # TRUE/FALSE value indicating whether specific 'environmental quality values' are supplied or not
                                        feedback.input = TRUE,          # TRUE/FALSE value indicating whether specific 'demongraphically dynamic feedback matrix' are supplied or not
                                        slope.input = TRUE,             # TRUE/FALSE value indicating whether specific 'slope matrix' are supplied or not
                                        bidirection.input = TRUE,       # TRUE/FALSE value indicating whether specific 'bidirectionality matrix' are supplied or not
                                        pool.list = pool,               # specific 'species trait values'
                                        env.list = patch_cond,          # specific 'environmental quality values' 
                                        feedback.matrix = feed,         # specific 'demongraphically dynamic feedback matrix' 
                                        slope.matrix = slope.array,     # specific 'slope matrix'
                                        bidirection.matrix = bidirect.array)   # specific 'bidirectionality matrix'

toc = Sys.time()

#### Run time and file size check of the demo simulation 
format(object.size(Sim), units="Mb")
toc - tic


###########################################################
#### Save the $output$comm of the demo simulation, this data frame has "location as rows" and "time as columns"
comm <- Sim$output$full_patch_time[, , 1]

#### Calculate richness through time
richness <- apply(comm, 2, function(x){length(unique(x[x!=0]))})         

##### Plot community richness through time
par(mar=c(2, 2, 1, 1), fig=c(0, 1, 0.8, 1)) 
plot(richness, type='l', ylim=c(0,Sim$param_list$pool.size), ylab="Richness")

#### Plot a colored map of the community species ID-location through time 
par(mar=c(2, 2, 1, 1), fig=c(0, 1, 0, 0.8), new=TRUE)
col.vec <- c('white', rainbow(Sim$param_list$pool.size))

plot(rep(1, Sim$param_list$patch.size),  (1:Sim$param_list$patch.size), pch=15, cex=0.2, col=col.vec[comm[, 1]+1], xlim=c(0, (Sim$param_list$time)+1), ann=FALSE, axes=FALSE)
for (i in 2:Sim$param_list$time){
  points(rep(i, Sim$param_list$patch.size), (1:Sim$param_list$patch.size), pch=15, cex=0.2, col=col.vec[comm[, i]+1])
}

abline(v=which(diff(richness)!=0), lwd = 0.2) # What time step is there a change in richness
axis(3, labels=F)
box()

