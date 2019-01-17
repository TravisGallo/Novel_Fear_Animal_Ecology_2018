package_load<-function(packages = c("reshape2"), quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# change the dataframe back to an array
df_2_array <- function(my_df = NULL){
  
  require(reshape2)
  my_array <- acast(my_df, species~site~season, value.var = "count")
  dimnames(my_array) <- NULL
  return(my_array)
}


betaABfromMeanSD = function( mean , sd ) {
  if ( mean <=0 | mean >= 1) stop("must have 0 < mean < 1")
  if ( sd <= 0 ) stop("sd must be > 0")
  kappa = mean*(1-mean)/sd^2 - 1
  if ( kappa <= 0 ) stop("invalid combination of mean and sd")
  a = mean * kappa
  b = ( 1.0 - mean ) * kappa
  return( list( a=a , b=b ) )
}

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      psinit = rbeta(2, 1, 1),
      pmu = rnorm(2, 0, 1),
      lp = rnorm(2, 0, 1),
      gmu = rnorm(2, 0, 1),
      gyr = matrix(rnorm(10*2, 0, 1), ncol= 2, nrow = 10),
      tau_yr = rgamma(2, 1, 1),
      gb = matrix(rnorm(2), ncol = 2, nrow = 1),
      pb = matrix(rnorm(2), ncol = 2, nrow = 1),
      lpb = matrix(rnorm(2), ncol = 2, nrow = 1),
      v0 = rnorm(1, 0, 1),
      v1 = rnorm(1, 0, 1), 
      vcoy0 = rnorm(1, 0, 1), 
      vcoy1 = rnorm(1, 0, 1),
      a = array(c(0,rnorm(1),0,rnorm(1),0,rnorm(1),0, rnorm(1)),dim=c(2,2,2)),
      da = matrix(c(0,rnorm(1),0, rnorm(1)), ncol=2, nrow=2),
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

overlap_mod <- function (A, B, xscale = 24, xcenter = c("noon", "midnight"), 
                         linetype = c(1,2), linecol = c("black", "blue"),
                         linewidth = c(1,1), olapcol = "lightgrey", rug = FALSE, 
                         extend = NULL, n.grid = 128, kmax = 3, adjust = 1, ...) 
{
  isMidnt <- match.arg(xcenter) == "midnight"
  bwA <- getBandWidth(A, kmax = kmax)/adjust
  bwB <- getBandWidth(B, kmax = kmax)/adjust
  if (is.na(bwA) || is.na(bwB)) 
    stop("Bandwidth estimation failed.")
  xsc <- if (is.na(xscale)) 
    1
  else xscale/(2 * pi)
  if (is.null(extend)) {
    xxRad <- seq(0, 2 * pi, length = n.grid)
  }
  else {
    xxRad <- seq(-pi/4, 9 * pi/4, length = n.grid)
  }
  if (isMidnt) 
    xxRad <- xxRad - pi
  xx <- xxRad * xsc
  densA <- densityFit(A, xxRad, bwA)/xsc
  densB <- densityFit(B, xxRad, bwB)/xsc
  densOL <- pmin(densA, densB)
  
  return(data.frame(x = xx, densityA = densA, densityB = densB))
}