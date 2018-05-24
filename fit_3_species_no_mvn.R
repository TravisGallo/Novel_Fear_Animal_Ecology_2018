# Load utility
source("coyote_deer_utility.R")
# Load needed packages and modules
package_load(c("reshape2", "runjags", "rjags", "parallel", "data.table","scales",
  "overlap","stringr","dplyr", "MCMCpack","denstrip"))
load.module('glm')


# Read in the species names
species_names <- read.table("species_used_in_sp10_sp13_analysis_6_1_17.txt", header = TRUE)

# Read in the site names
site_names <- read.table("sites_used_in_sp10_sp13_analysis_6_1_17.txt", header = TRUE)
site_names <- site_names[-23,] 


# Function that sets up data for JAGs run
vig_JAGS_setup <- function(index, model) {
  
  # Read in z-matrix for the two species
  z <- df_2_array(read.table("z_matrix_sp10_sp13_6_1_17.txt", 
    header = TRUE, sep = "\t"))[index,,-c(1,2)]
  

  # Build y-array and j-matrix
  y_array <- df_2_array(read.table("y_matrix_sp10_sp13_6_1_17.txt", 
    header = TRUE, sep = "\t"))[index,,-c(1,2)]
  j_mat <- read.table("j_matrix_sp10_sp13_6_1_17.txt", 
    header = TRUE, sep = "\t")[,-c(1,2)]
  
  # read in vigilance data
  vig_deer <- df_2_array(read.csv("deer_vigilance.csv", header = TRUE))
  hups_deer <- vig_deer[1,,]
  
  phot_deer <- vig_deer[2,,]
  phot_deer[is.na(phot_deer)] <- 0.001
  
  phot_bin_deer <- phot_deer
  phot_bin_deer[phot_bin_deer > 1] <- 1
  # read in vigilance data
  vig_rab <- df_2_array(read.csv("rabbit_vigilance.csv", header = TRUE))
  hups_rab <- vig_rab[1,,]
  
  phot_rab <- vig_rab[2,,]
  phot_rab[is.na(phot_rab)] <- 0.001
  
  phot_bin_rab <- phot_rab
  phot_bin_rab[phot_bin_rab > 1] <- 1
  
  #removing any site with only 1 season of data
  tg <- which(rowSums(is.na(z[1,,]))>9)
  
  z <- z[,-tg,]
  y <- y_array[,-tg,]
  j <- j_mat[-tg,]
  hups_deer <- hups_deer[-tg,]
  phot_deer <- phot_deer[-tg,]
  phot_bin_deer <- phot_bin_deer[-tg,]
  hups_rab <- hups_rab[-tg,]
  phot_rab <- phot_rab[-tg,]
  phot_bin_rab <- phot_bin_rab[-tg,]
  
  # turn heads up data to long form
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
  
  # create pc of habitat covs
  covdat <- read.csv("urban_covs.csv", header = TRUE)
  hab <- as.data.frame(covdat[-tg,3:5])
  pc <- prcomp(hab, scale. = TRUE)
  
  # set up data for JAGS
  #y[2,,] <- NA
  #z[2,,] <- NA
  data_list <- list(
    y = y,
    nyear = dim(z)[3], 
    nsite = dim(z)[2], 
    nspec = dim(z)[1],
    J = as.matrix(j),
    pcov = pc$x[,1],
    lpcov = pc$x[,1],
    vcov = pc$x[,1],
    inxscov = as.matrix(data.frame(a = 1, b = pc$x[,1])),
    ncov_phi = 1,
    phot_deer = phot_deer,
    phot_rab = phot_rab,
    hup_deer = hups_long_deer,
    hup_rab = hups_long_rab,
    hups_loc_deer=as.matrix(hups_loc_deer),
    n_event_deer=nrow(hups_loc_deer),
    hups_loc_rab=as.matrix(hups_loc_rab),
    n_event_rab=nrow(hups_loc_rab)
    )
  
  # set up initial values
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
  
  tmon <- c("v0", "vig_urb", "vcoy0", "vcoy_urb", "SIF_deer", "SIF_deer_high", 
  	"SIF_deer_low","SIF_rab",  "SIF_rab_high", "SIF_rab_low", "a", "da", "psi0",
    "psi_urb", "pmu", "lp", "gmu", "gyr", "gb", "pb", "lpb", "z")
  
  
  #set.seed(1234)
  
  mout <- run.jags( model = model, 
    monitor = tmon , 
    data = data_list ,  
    inits = inits , 
    n.chains = detectCores()-7 ,
    adapt = 500,#10000,
    burnin = 500,#10000 , 
    sample = ceiling(1000/(detectCores()-7)) ,
    thin = 2,
    summarise = FALSE ,
    plots = FALSE,
    method = "parallel")
  
  # summarize model
  mod_sum <- add.summary(mout)
  # pull out mcmc chains of all parameters that we tracked
  mmat <- as.matrix(as.mcmc.list(mout), chains = TRUE)
  
  return(list(mmat=mmat, model_results=mout, mod_sum=mod_sum))
  
}

# Create data indexs for species combos

three <- c(1,2,7)

# run models
hunk <- vig_JAGS_setup(three, "three_species.R")

coyt_deer_mod <- vig_JAGS_setup(coyt_deer,"deer", "test_model2.R")
coyt_rabbit_mod <- vig_JAGS_setup(coyt_rabbit,"rabbit", "test_model2.R")


m1_coy <- t(apply(coyt_deer_mod[[1]], 2, quantile, probs = c(0.025,0.5,0.975)))
m2_coy <- t(apply(coyt_rabbit_mod[[1]], 2, quantile, probs = c(0.025,0.5,0.975)))
cbind(m1_coy, m2_coy)

m1 <- cor(coyt_deer_mod[[1]])
# check logit scale parameters
coyt_deer_mod$mod_sum
coyt_rabbit_mod$mod_sum



# summarize vigilance data
deer_sum <- apply(coyt_deer_mod$mmat, 2, quantile, probs=c(0.025, 0.5,0.975))
colnames(coyt_deer_mod$mmat)
deer_sum[,1:82]
rabbit_sum <- apply(coyt_rabbit_mod$mmat, 2, quantile, probs=c(0.025, 0.5,0.975))

# set up data for prediction plots

pred_vig <- function (mmat){
  
  # bring back urbanization covariate
  covdat <- read.csv("urban_covs.csv", header = TRUE)
  hab <- as.data.frame(covdat[-23,3:5])
  pc <- prcomp(hab, scale. = TRUE)
  x <- range(pc$x)
  co <- seq(x[1], x[2], 0.1)
  
  # without coyote
  xmat <- cbind(1, co)
  y <- xmat %*% t(mmat[,c(2,3)])
  y2 <- apply(y, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
  # backtransform
  ex <- function(x) 1/(1+exp(-x))
  
  yp <- ex(y2)
  
  # with coyote
  xmat2 <- cbind(xmat, 1)
  w <- xmat2 %*% t(mmat[,c(2,3,4)])
  w2 <- apply(w, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
  # backtransform
  wp <- ex(w2)
  
  # predict occupancy
  gam <- xmat %*% t(mmat[,c(36,58)])
  col <- ex(gam)
  phi <- xmat %*% t(mmat[,c(32,60)])
  ext <- 1-ex(phi)
  occu <- col/(col+ext)
  occu2 <- apply(occu, 1, quantile, probs = c(0.025, 0.5, 0.975))
  
  return(list(co=co,yp=yp,wp=wp,occu=occu2))
}

pred_deer_vig <- pred_vig(coyt_deer_mod$mmat)
pred_rabbit_vig <- pred_vig(coyt_rabbit_mod$mmat)

######################
## Overlap Analysis ##
######################

# function that will extract overlap information for each prey species
overlap_extract <- function (mmat, spB) {  
  
  # pull out z matrix to use in overlap analysis
  zmat <- mmat[,grep("z", colnames(mmat))]
  coyote_z <- zmat[,seq(1,dim(zmat)[2]-1,2)]
  
  # set up photo data to run through coyote z-matrix
  # create season indicator
  seasons <- c("FA10","WI11","SP11","SU11","FA11","WI12","SP12","SU12","FA12","WI13","SP13")
  # load and summarize photo data
  data_pics <- read.csv("activity.data.csv")
  # only use data from 2010-2013
  data_pics$yearA <- str_sub(data_pics$PictureDate,-2,-1)
  data_pics <- data_pics[-which(data_pics$year>13|data_pics$year<10|data_pics$year==" "),]
  # convert minutes to radians
  data_pics$radians <- (data_pics$Minutes/1440)*2*pi
  data_pics$Species <- trimws(data_pics$Species)
  data_pics$Species <- tolower(data_pics$Species)
  spB <- tolower(spB)
  
  # pull out prey only data
  if(spB == "rabbit") {
    data_prey <- data_pics[which(data_pics$Species==spB),]
    # removing jacked up data from Austin Park - see exploratory analysis below
    data_prey <- data_prey[-which(data_prey$SurveyID == "R03-AUP1-WI11"),]
  } else {
    data_prey <- data_pics[which(data_pics$Species==spB),]
  }
  
  # loop through to sort data
  # create activity densities for deer with and without coyotoes
  site_names <- read.table("sites_used_in_sp10_sp13_analysis_6_1_17.txt", header = TRUE)
  site_names <- site_names[-23,]
  overlap_est <- rep(0,nrow(coyote_z))
  overlap_dens <- vector("list",nrow(coyote_z))
  for(i in 1:nrow(coyote_z)){
    # sites where coyotes were present or absent
    z_t <- matrix(coyote_z[i,],ncol=11, nrow=102)
    rownames(z_t) <- site_names
    colnames(z_t) <- seasons
    ind <- which(z_t==1,arr.ind=TRUE)
    SurveyID <- paste(rownames(z_t)[ind[,1]], colnames(z_t)[ind[,2]], sep="-")
    # process photo data
    has_both_pics <- data_prey[which(data_prey$SurveyID %in% SurveyID),]
    if(sum(is.na(has_both_pics$radians)) > 0){
      has_both_pics <- has_both_pics[-which(is.na(has_both_pics$radians)==TRUE),]
    }
    alone_pics <- data_prey[-which(data_prey$SurveyID %in% SurveyID),]
    if(sum(is.na(alone_pics$radians)) > 0){
      alone_pics <- alone_pics[-which(is.na(alone_pics$radians)==TRUE),]
    }
    # overlap estimation
    overlap_est[i] <-overlapEst(alone_pics[,7], has_both_pics[,7])[2]
    overlap_dens[[i]] <- overlap_mod(alone_pics[,7], has_both_pics[,7])
    overlap_dens[[i]]$densityAB <- apply(overlap_dens[[i]][,2:3], 1, min)
  }
  
  return(list(overlap_est=overlap_est, overlap_dens=overlap_dens))
  
}

deer_overlap <- overlap_extract(coyt_deer_mod$mmat, "Deer")
rabbit_overlap <- overlap_extract(coyt_rabbit_mod$mmat, "Rabbit")

# extract coyote data to include in plots
data_pics <- read.csv("activity.data.csv")
# only use data from 2010-2013
data_pics$yearA <- str_sub(data_pics$PictureDate,-2,-1)
data_pics <- data_pics[-which(data_pics$year>13|data_pics$year<10|data_pics$year==" "),]
# convert minutes to radians
data_pics$radians <- (data_pics$Minutes/1440)*2*pi
data_pics$Species <- trimws(data_pics$Species)
coyt_activity <- data_pics[which(data_pics$Species=="Coyote"),]
coyt_rads=coyt_activity[,7]
coyt_rads=coyt_rads[!is.na(coyt_rads)]

###########################################################
### Exploratory stuff to find a weird thing in the data ###
######## resulted in if-else statement above ##############
###########################################################


# there seems to be a phenomenon in the data - lets explore
# grab difference between run 1 and 2

# ran these after changing i to 1 and 2
#has_both_pics1 <- data_rabbit[which(data_rabbit$SurveyID %in% SurveyID),]
#has_both_pics2 <- data_rabbit[which(data_rabbit$SurveyID %in% SurveyID),]

# find the difference
#diff_sites <- as.factor(unique(c(as.character(has_both_pics1$SurveyID[!has_both_pics1$SurveyID %in% has_both_pics2$SurveyID]),
#as.character(has_both_pics2$SurveyID[!has_both_pics2$SurveyID %in% has_both_pics1$SurveyID]))))

# explore those sites
#weird <- data_rabbit[which(data_rabbit$SurveyID %in% diff_sites & between(data_rabbit$radians, 2,4)),]
#weird_sites <- unique(weird$SurveyID)

#AUS <- weird[which(weird$SurveyID == "R03-AUP1-WI11"),]
#unique(AUS$PictureDate)
#range(AUS$Minutes)

# seems to be this Austin Park winter 2011. 
# Camera malfunction and a lot of blank pictures identified as rabbit. We will remove.


