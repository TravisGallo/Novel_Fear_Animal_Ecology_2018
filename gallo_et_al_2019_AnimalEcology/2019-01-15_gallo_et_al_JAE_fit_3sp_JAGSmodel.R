#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## FITTING 3-SPECIES OCCUPANCY MODEL, BINOMIAL VIGILANCE MODEL, AND OVERLAP ANALYSIS ##
## CODE ASSOCIATED WITH GALLO ET AL 2019 ANIMAL ECOLOGY##
## URBANIZATION ALTERS PREDATOR AVOIDANCE BEHAVIOR ##
## BY TRAVIS GALLO and MASON FIDINO ##
## Urban Wildlife Institute ##
## Last Updated 2019-01-15 ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Load utility
source("2019-01-15_gallo_et_al_JAE_utility_script.R")
# Load needed packages and modules
package_load(c("reshape2", "runjags", "rjags", "parallel", "data.table","scales",
               "overlap","stringr","dplyr", "MCMCpack","denstrip"))
load.module('glm')

## Read in data
# Read in the species names
species_names <- read.table("./Data/2019-01-15_gallo_et_al_JAE_species_used_in_sp10_sp13_analysis.txt", header = TRUE)

# Read in the site names
site_names <- read.table("./Data/2019-01-15_gallo_et_al_JAE_sites_used_in_sp10_sp13_analysis.txt", header = TRUE)
# We do not use data from this site
site_names <- site_names[-23,] 


## Set up data for three species co-occurence model and vigilance binomial model
# Function that sets up data for JAGs run
vig_JAGS_setup <- function(index, model) {
  
  # Read in z-matrix for the two species
  z <- df_2_array(read.table("./Data/2019-01-15_gallo_et_al_JAE_z_matrix_sp10_sp13.txt", 
    header = TRUE, sep = "\t"))[index,,-c(1,2)]
  
  # Build y-array and j-matrix
  y_array <- df_2_array(read.table("./Data/2019-01-15_gallo_et_al_JAE_y_matrix_sp10_sp13.txt", 
    header = TRUE, sep = "\t"))[index,,-c(1,2)]
  j_mat <- read.table("./Data/2019-01-15_gallo_et_al_JAE_j_matrix_sp10_sp13.txt", 
    header = TRUE, sep = "\t")[,-c(1,2)]
  
  # Read in vigilance data for deer
  vig_deer <- df_2_array(read.csv("./Data/2019-01-15_gallo_et_al_JAE_deer_vigilance.csv", header = TRUE))
  # Heads up photos
  hups_deer <- vig_deer[1,,]
  # Total photos
  phot_deer <- vig_deer[2,,]
  #  phot_deer cannot have NA values as it is not a random variable
  #  supplied to JAGS. 
  phot_deer[is.na(phot_deer)] <- 0.001
  
  # Read in vigilance data for rabbit
  vig_rab <- df_2_array(read.csv("./Data/2019-01-15_gallo_et_al_JAE_rabbit_vigilance.csv", header = TRUE))
  # Head up photos
  hups_rab <- vig_rab[1,,]
  # Total photos
  phot_rab <- vig_rab[2,,]
  #
  phot_rab[is.na(phot_rab)] <- 0.001
  
  # Removing any site with only 1 season of data
  tg <- which(rowSums(is.na(z[1,,]))>9)
  # Remove these sites from all data sets
  z <- z[,-tg,]
  y <- y_array[,-tg,]
  j <- j_mat[-tg,]
  hups_deer <- hups_deer[-tg,]
  phot_deer <- phot_deer[-tg,]
  hups_rab <- hups_rab[-tg,]
  phot_rab <- phot_rab[-tg,]
  
  # Convert head up data to long form
  hups_loc_deer <- which(!is.na(hups_deer)==TRUE, arr.ind=TRUE)
  hups_long_deer <- rep(0,nrow(hups_loc_deer))
  for(i in 1:length(hups_long_deer)){
    hups_long_deer[i] <- hups_deer[hups_loc_deer[i,1],hups_loc_deer[i,2]]
  }
  hups_loc_rab <- which(!is.na(hups_rab)==TRUE, arr.ind=TRUE)
  hups_long_rab <- rep(0,nrow(hups_loc_rab))
  for(i in 1:length(hups_long_rab)){
    hups_long_rab[i] <- hups_rab[hups_loc_rab[i,1],hups_loc_rab[i,2]]
  }
  
  # Read in site-level habitat covariates
  covdat <- read.csv("./Data/2019-01-15_gallo_et_al_JAE_urban_covs.csv", header = TRUE)
  # Remove sites with only one season of data
  hab <- as.data.frame(covdat[-tg,3:5])
  # Calculate principle component of habitat covariates to be an index of urbanization
  pc <- prcomp(hab, scale. = TRUE)
  
  # Set up data for JAGS
  data_list <- list(
    y = y, # detection / non-detection data
    nyear = dim(z)[3], # number of seasons sampled
    nsite = dim(z)[2], # number of sites sampled
    nspec = dim(z)[1], # number of species
    J = as.matrix(j),  # number of days sampled per site and season
    pcov = pc$x[,1],   # urb covariate for persistance
    lpcov = pc$x[,1],  # urb covariate for detection
    vcov = pc$x[,1],   # urb covariate for vigilance
    inxscov = as.matrix(data.frame(a = 1, b = pc$x[,1])), # design matrix for inxs
    ncov_phi = 1, # number of slope terms for persistence
    phot_deer = phot_deer, # number of deer photos
    phot_rab = phot_rab,   # number of rabbit photos
    hup_deer = hups_long_deer, # number of head up deer photos
    hup_rab = hups_long_rab, # number of head up rabbit photos
    hups_loc_deer=as.matrix(hups_loc_deer), # connects long format data to site
    n_event_deer=nrow(hups_loc_deer), # used to loop through vigilance analysis
    hups_loc_rab=as.matrix(hups_loc_rab), # connects long format data to site
    n_event_rab=nrow(hups_loc_rab) # used to loop through vigilance analysis
    )
  
  # Set up initial values
  inits <- function(chain){
    gen_list <- function(chain = chain){
      list( 
        z = z,
        psi0 = rnorm(3, 0, 1),
        psi_urb = rnorm(3, 0, 1),
        pmu = rnorm(3, 0, 1),
        lp = rnorm(3, 0, 1),
        gmu = rnorm(3, 0, 1),
        gyr = matrix(rnorm(10*3, 0, 1), ncol= 3, nrow = 10),
        tau_yr = rgamma(3, 1, 1),
        v0 = rnorm(2, 0, 1),
        vcoy0 = rnorm(2, 0, 1), 
      	vcoy_urb = rnorm(2, 0, 1),
        vig_urb = rnorm(2),
      	ppar = rnorm(4),
      	gpar = rnorm(4),
      	dpar = rnorm(4),
        .RNG.name = switch(chain,
          "1" = "base::Wichmann-Hill",
          "2" = "base::Marsaglia-Multicarry",
          "3" = "base::Super-Duper",
          "4" = "base::Mersenne-Twister",
          "5" = "base::Wichmann-Hill",
          "6" = "base::Marsaglia-Multicarry",
          "7" = "base::Super-Duper",
          "8" = "base::Mersenne-Twister"),
        .RNG.seed = sample(1:1e+06, 1)
      )
    }
    return(switch(chain,           
      "1" = gen_list(chain),
      "2" = gen_list(chain),
      "3" = gen_list(chain),
      "4" = gen_list(chain),
      "5" = gen_list(chain),
      "6" = gen_list(chain),
      "7" = gen_list(chain),
      "8" = gen_list(chain)
    )
    )
  }
  
  # Parameters to monitor
  tmon <- c("v0", "vig_urb", "vcoy0", "vcoy_urb", "SIF_deer", "SIF_deer_high", 
  	"SIF_deer_low","SIF_rab",  "SIF_rab_high", "SIF_rab_low", "a", "da", "psi0",
    "psi_urb", "pmu", "lp", "gmu", "gyr", "tau_yr", "gb", "pb", "lpb", "z")
  
  # Run model
  mout <- run.jags( model = model, 
    monitor = tmon , 
    data = data_list ,  
    inits = inits , 
    n.chains = 2 ,
    adapt = 10000,
    burnin = 10000, 
    sample = 50000,
    thin = 2,
    summarise = FALSE ,
    plots = FALSE,
    method = "parallel")
  
  # Summarize model
  mod_sum <- add.summary(mout)
  
  # MCMC chains of all parameters that we tracked
  mmat <- as.matrix(as.mcmc.list(mout), chains = TRUE)
  
  return(list(mmat=mmat, model_results=mout, mod_sum=mod_sum))
  
}

# Data index for species combos
three <- c(1,2,7)

# Run 3 species model
mod_3sp <- vig_JAGS_setup(three, "2019-01-15_gallo_et_al_JAE_3sp_interaction_JAGSmodel.R")

# Summarize model results
mod_sum <- t(apply(mod_3sp$mmat[,1:141], 2, quantile, probs = c(0.025,0.5,0.975)))
taus <- t(apply(1/sqrt(mod_3sp$mmat[,133:135]), 2, quantile, probs = c(0.025,0.5,0.975)))


## Set up data for vigilance visualization

# Function to backtransform logit parameters below
ex <- function(x) 1/(1+exp(-x))

# Function to predict vigilance for Figure 2 visual
pred_vig <- function (mmat, index){
  
  # Bring back urbanization covariate
  covdat <- read.csv("./Data/2019-01-15_gallo_et_al_JAE_urban_covs.csv", header = TRUE)
  hab <- as.data.frame(covdat[-23,3:5])
  pc <- prcomp(hab, scale. = TRUE)
  # Calculate minimum and maximum to predict between
  x <- range(pc$x)
  # Sequence of covariate values to predict across
  co <- seq(x[1], x[2], 0.1)
  
  # Without coyote
  # Create matrix with intercept as 1 and our covariate sequence as following column
  xmat <- cbind(1, co)
  # Matrix multiplication with new covariate sequence and coefficient values
  y <- xmat %*% t(mmat[,c(paste0("v0[",index,"]"), paste0("vig_urb[",index,"]"))])
  # Summarize
  y2 <- apply(y, 1, quantile, probs = c(0.025, 0.5, 0.975))
  # Backtransform
  yp <- ex(y2)
  
  # With coyote - same process as above now we include a 1 that indicates coyotes present
  w <- xmat %*% t(mmat[,c(paste0("vcoy0[",index,"]"), paste0("vcoy_urb[",index,"]"))])
  w2 <- apply(w, 1, quantile, probs = c(0.025, 0.5, 0.975))
  # Backtransform
  wp <- ex(w2)
  
  # Predict coyote occupancy
  gam <- xmat %*% t(mmat[,c(paste0("gmu[1]"),paste0("gb[1]"))])
  col <- ex(gam) # backtransform logit parameter
  phi <- xmat %*% t(mmat[,c(paste0("pmu[1]"),paste0("pb[1]"))])
  ext <- 1-ex(phi) # backtransform logit parameter
  occu <- col/(col+ext)
  # Summarize
  occu2 <- apply(occu, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
  return(list(co=co,yp=yp,wp=wp,occu=occu2))
}

# Run function
pred_deer_vig <- pred_vig(mod_3sp$mmat, 1)
pred_rabbit_vig <- pred_vig(mod_3sp$mmat, 2)

# Set up data to plot vigilance data points with best fit lines in Figure 2c
vig_dat <- function(vig_sp){
  
  # Load z matrix to manage data below
  z <- df_2_array(read.table("./Data/2019-01-15_gallo_et_al_JAE_z_matrix_sp10_sp13.txt", 
                             header = TRUE, sep = "\t"))[c(1,2,7),,-c(1,2)]
  # Sites to remove
  tg <- which(rowSums(is.na(z[1,,]))>9)
  
  # Coyote occupancy estimated from model (z)
  coyote_z <- mod_3sp$mmat[,grep("^z\\[1", colnames(mod_3sp$mmat))]
  # Calculate mean occupancy rate
  coyt_occu <- apply(coyote_z, 2, mean)
  coyt_occu_mat <- matrix(coyt_occu, nrow=102, ncol=11, byrow=FALSE)
  # Indicated sites and seasons that have a coyote occupancy > 0.75
  coyt_occu_mat[which(coyt_occu_mat > 0.74)] <- 1
  coyt_occu_mat[which(coyt_occu_mat < 1)] <- 0
  
  # Proportion of heads up photos at each site and season
  vig <- df_2_array(read.csv(paste("./Data/2019-01-15_gallo_et_al_JAE_",vig_sp,"_vigilance.csv", 
                                   sep=""), header = TRUE))
  
  # Urbanization covariate - same as above
  covdat <- read.csv("./Data/2019-01-15_gallo_et_al_JAE_urban_covs.csv", header = TRUE)
  covdat <- covdat[-tg,]
  hab <- as.data.frame(covdat[,3:5])
  pc <- prcomp(hab, scale. = TRUE)
  # Bin urbanization covariate
  urb_bin <- cut(pc$x[,1], seq(-5,5,0.5))
  levels(urb_bin) <- seq(-4.75,4.75, 0.5)
  urb_bin <- as.numeric(as.character(urb_bin))
  
  # Create dataframe to work with
  vig_dat <- data.frame(hup = as.numeric(vig[1,-tg,]), tot = as.numeric(vig[2,-tg,]),
                        urb = rep(urb_bin, 11), coyt = as.numeric(coyt_occu_mat))
  # Group by coyotes present and the urbanization bin
  vig_results <- vig_dat %>% 
    group_by(urb, coyt) %>% 
    summarize(hup=sum(hup, na.rm=TRUE), tot = sum(tot, na.rm = TRUE)) %>% 
    as.data.frame
  # Calculate proportions of head up photos
  vig_results$prop <- vig_results$hup/vig_results$tot
  # Calculate the complement (head down)
  vig_results$q <- 1 - vig_results$prop
  # Calculate SE (this is only included in supplemental material)
  vig_results$se <- with(vig_results, {
    sqrt((hup*prop*q)/tot)})
  # Calculate upper and lower bounds
  vig_results$hi <- vig_results$prop + vig_results$se
  vig_results$hi[vig_results$hi > 1] <- 1
  vig_results$lo <- vig_results$prop - vig_results$se
  vig_results$lo[vig_results$lo < 0] <- 0
  
  return(vig_results)
  
}

# Extract results
vig_rab <- vig_dat("Rabbit")
vig_deer <- vig_dat("Deer")


## Overlap Analysis

# Function that will calculated overlap of daily activity for each prey species
overlap_extract <- function (mmat, spB) {  
  
  # z matrix estimated from occupancy model to use in overlap analysis
  coyote_z <- mmat[,grep("^z\\[1", colnames(mmat))]
  
  # Set up photo data to run through coyote z-matrix
  # Create sampling season indicator
  seasons <- c("FA10","WI11","SP11","SU11","FA11","WI12","SP12","SU12","FA12","WI13","SP13")
  
  # Load and summarize photo data
  data_pics <- read.csv("./Data/2019-01-15_gallo_et_al_JAE_activity.data.csv")
  
  # Only use data from 2010-2013
  data_pics$yearA <- str_sub(data_pics$PictureDate,-2,-1)
  data_pics <- data_pics[-which(data_pics$year>13|data_pics$year<10|data_pics$year==" "),]
  
  # Convert minutes to radians
  data_pics$radians <- (data_pics$Minutes/1440)*2*pi
  data_pics$Species <- trimws(data_pics$Species)
  data_pics$Species <- tolower(data_pics$Species)
  # Make entry for spB be lowercase to match data
  spB <- tolower(spB)
  
  # Extract prey only data
  if(spB == "rabbit") {
    data_prey <- data_pics[which(data_pics$Species==spB),]
    # remove inconsistent data from Austin Park
    data_prey <- data_prey[-which(data_prey$SurveyID == "R03-AUP1-WI11"),]
  } else {
    data_prey <- data_pics[which(data_pics$Species==spB),]
  }
  
  # Loop through to sort data and create activity density for prey species when coyote
  # were present and absent
  # Read back in site names
  site_names <- read.table("./Data/2019-01-15_gallo_et_al_JAE_sites_used_in_sp10_sp13_analysis.txt", 
                           header = TRUE)
  # We do not use data from this site
  site_names <- site_names[-23,]
  
  # Create objects to store data
  overlap_est <- rep(0,nrow(coyote_z))
  overlap_dens <- vector("list",nrow(coyote_z))
  for(i in 1:nrow(coyote_z)){
    
    # Sites where coyotes were present or absent
    z_t <- matrix(coyote_z[i,],ncol=11, nrow=102)
    rownames(z_t) <- site_names
    colnames(z_t) <- seasons
    
    # Indicates when and where a coyote was present
    ind <- which(z_t==1,arr.ind=TRUE)
    # Index by sites and seasons where coyotes were present
    # Put them together seperated by a dash to match our site/season naming convention
    SurveyID <- paste(rownames(z_t)[ind[,1]], colnames(z_t)[ind[,2]], sep="-")
    
    # Process photo data
    # Sites and seasons where prey and coyote were captured together
    has_both_pics <- data_prey[which(data_prey$SurveyID %in% SurveyID),]
    if(sum(is.na(has_both_pics$radians)) > 0){
      has_both_pics <- has_both_pics[-which(is.na(has_both_pics$radians)==TRUE),]
    }
    # Sites and season where a coyote was not captured
    alone_pics <- data_prey[-which(data_prey$SurveyID %in% SurveyID),]
    if(sum(is.na(alone_pics$radians)) > 0){
      alone_pics <- alone_pics[-which(is.na(alone_pics$radians)==TRUE),]
    }
    
    # For each scenerio calculate a daily activity density
    # Calculate a proportion of overlap (d-hat = 4) at each MCMC step
    overlap_est[i] <-overlapEst(alone_pics[,7], has_both_pics[,7])[2]
    # Modified overlapPlot function to only return density values and not a plot
    overlap_dens[[i]] <- overlap_mod(alone_pics[,7], has_both_pics[,7])
    # Need the minimum of the two to create polygon when plotting
    overlap_dens[[i]]$densityAB <- apply(overlap_dens[[i]][,2:3], 1, min)
  }
  
  # Extract coyote activity in radians
  coyt_activity <- data_pics[which(data_pics$Species=="coyote"),]
  coyt_rads=coyt_activity[,7]
  coyt_rads=coyt_rads[!is.na(coyt_rads)]
  
  return(list(overlap_est=overlap_est, 
              overlap_dens=overlap_dens, 
              coyt_rads=coyt_rads))
  
}

# Run analysis for each prey species
deer_overlap <- overlap_extract(mod_3sp$mmat, "Deer")
rabbit_overlap <- overlap_extract(mod_3sp$mmat, "Rabbit")

# Extract coyote data to include in plots
data_pics <- read.csv("./Data/2019-01-15_gallo_et_al_JAE_activity.data.csv")
# Only use data from 2010-2013
data_pics$yearA <- str_sub(data_pics$PictureDate,-2,-1)
data_pics <- data_pics[-which(data_pics$year>13|data_pics$year<10|data_pics$year==" "),]
# Convert minutes to radians - will be used for visualization
data_pics$radians <- (data_pics$Minutes/1440)*2*pi
data_pics$Species <- trimws(data_pics$Species)
coyt_activity <- data_pics[which(data_pics$Species=="Coyote"),]
coyt_rads=coyt_activity[,7]
coyt_rads=coyt_rads[!is.na(coyt_rads)]


