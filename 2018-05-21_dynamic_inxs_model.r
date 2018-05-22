model{
####################################
### PRIORS #########################
####################################


# persistance mvn URB parameters
# using scaled inverse wishart
# 1 = deer/urbanization without coyote covariate
# 2 = deer/urbanization with coyote covariate
per_urb[1:2] ~ dmnorm(mus, tau_per[,])
tau_per[1:2,1:2] ~ dwish(W[,], 3)
sigma_perm[1:2,1:2] <- inverse(tau_per[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_per[k,k.prime] <- sigma_perm[k,k.prime]/
      sqrt(sigma_perm[k,k] * sigma_perm[k.prime,k.prime])
  }
  sigma_per[k] <- sqrt(sigma_perm[k,k])
}
##
# inxs persistence
## species x parameter x process array
a[1,1,2] <- 0
a[1,2,2] <- 0
a[2,1,2] ~ dt(0, 2.5, 1) # intercept
a[2,2,2] <- per_urb[2] # urb, from mvn
##
# persistence covariates
##
pb[1] ~ dt(0,2.5,1) # coyote
pb[2] <- per_urb[1] # deer, from mvn
#
# colonization mvn urb parameters
# using scaled inverse wishart
col_urb[1:2] ~ dmnorm(mus, tau_col[,])
tau_col[1:2,1:2] ~ dwish(W[,], 3)
sigma_coln[1:2,1:2] <- inverse(tau_col[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_col[k,k.prime] <- sigma_coln[k,k.prime]/
      sqrt(sigma_coln[k,k] * sigma_coln[k.prime,k.prime])
  }
  sigma_col[k] <- sqrt(sigma_coln[k,k])
}
##
# inxs colonization
##
a[1,1,1] <- 0
a[1,2,1] <- 0
a[2,1,1] ~ dt(0, 2.5, 1) # intercept
a[2,2,1] <- col_urb[2] # urb, from mvn
##
# colonization covariates
##

  gb[1] ~ dt(0,2.5,1) # coyote
  gb[2] <- col_urb[1] # deer, from mvn

# detection mvn urb parameters
# using scaled inverse wishart
det_urb[1:2] ~ dmnorm(mus, tau_det[,])
tau_det[1:2,1:2] ~ dwish(W[,], 3)
sigma_detn[1:2,1:2] <- inverse(tau_det[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_det[k,k.prime] <- sigma_detn[k,k.prime]/
      sqrt(sigma_detn[k,k] * sigma_detn[k.prime,k.prime])
  }
  sigma_det[k] <-sqrt(sigma_detn[k,k])
}
##
# urb coyote on deer prior detection state
##
da[1,1] <- 0
da[1,2] <- 0
da[2,1] ~ dt(0, 2.5, 1) # intercept
da[2,2] <- det_urb[2] # urb, from mvn
##
# covariates on detection
##
lpb[1] ~ dt(0,2.5,1) # coyote
lpb[2] <- det_urb[1] # deer from mvn normal
# vigilance mvn urb parameters
# using scaled inverse wishart
vig_urb ~ dt(0, 2.5, 1)

##
# intercepts on vigilance
##
v0 ~ dt(0, 2.5, 1)
# coyote presence on deer vigilance
##
vcoy0 ~ dt(0, 2.5, 1)
for(i in 1:2){
##
# priors for intercepts
##
psi0[i] ~ dt(0, 2.5, 1) # initial occupancy
psi_urb[i] ~ dt(0, 2.5, 1)
pmu[i] ~ dt(0, 2.5,1)               # persistence
lp[i] ~ dt(0, 2.5, 1)               # detection
gmu[i] ~ dnorm(0, 0.3)              # colonization mu
##
# temporal random effect on colonization
##
for(t in 1:(nyear-1)){
	gyr[t,i] ~  dnorm(gmu[i], tau_yr[i])
}
tau_yr[i] ~ dgamma(1,1)
}

####################################
### LATENT STATE ###################
####################################
for(i in 1:nspec){
	for(k in 1:nsite){
	##
	# Initial occupancy, Intercept only.
	##
	logit(psi[i,k,1]) <- psi0[i] + psi_urb[i]*pcov[k]
	z[i,k,1] ~ dbern(psi[i,k,1])
	for(t in 2:nyear){
	##
	# Mixture of colonization and persistence
	##
	logit(psi[i,k,t]) <- (z[i,k,t-1] * (pmu[i] + (inprod(a[i,,2], inxscov[k,]) * z[1,k,t-1]) + pb[i]*pcov[k]))+
				    ((1 - z[i,k,t-1]) * (gyr[t-1,i] + (inprod(a[i,,1], inxscov[k,]) * z[1,k,t-1]) + gb[i]*pcov[k]))
	##
	# Likelihood
	##
	z[i,k,t] ~ dbern(psi[i,k,t])
		}
	}
}
####################################
### VIGILANCE ANALYSIS #############
####################################
  for(k in 1:n_event){
      logit(vig_prob[k]) <- v0 + vig_urb * vcov[hups_loc[k,1]] + vcoy0 * z[1,hups_loc[k,1],hups_loc[k,2]]
      hup[k] ~ dbin(vig_prob[k], phot[hups_loc[k,1],hups_loc[k,2]])
  }
####################################
### OBSERVATIONAL PROCESS ##########
####################################
for(i in 1:nspec){
	for(k in 1:nsite){
		for(t in 1:nyear){
		  logit(dprob[i,k,t]) <- lp[i] + lpb[i]*lpcov[k] + (inprod(da[i,], inxscov[k,]) * z[1,k,t])
		mu[i,k,t] <- z[i,k,t] * dprob[i,k,t]
		y[i,k,t] ~ dbin(mu[i,k,t], J[k,t])
		}
	}
}
####################################
### DERIVED PARAMETERS #############
####################################
#  PSIA = Coyote
#  PSIB = Deer
# PSIAB = Deer conditional on coyote
# convert logit scale parameters to 
# probabilities
psinit[1] <- ilogit(psi0[1])
psinit[2] <- ilogit(psi0[2])
pmup[1] <- ilogit(pmu[1])
pmup[2] <- ilogit(pmu[2])
cdp <- ilogit(a[2,1,2] + pmu[2])
for(k in 1:nspec){
	for(t in 1:(nyear-1)){
	gyrp[t,k] <- ilogit(gyr[t,k])
	cdg[t,k] <- ilogit(a[2,1,1] + gyr[t,k])
	}
}
PSIA[1] <- psinit[1] * (pmup[1]) + (1 - psinit[1]) * gyrp[1,1]
PSIBa[1] <- psinit[2] * (pmup[2]) + (1 - psinit[2]) * gyrp[1,2]
PSIBA[1] <- psinit[2] * (cdp) + (1 - psinit[2]) * (cdg[1,2])
PSIB[1] <- PSIA[1] * PSIBA[1] + (1 - PSIA[1]) * PSIBa[1]
PSIAB[1] <- PSIA[1] * PSIBA[1]
SIF[1] <- PSIAB[1]/(PSIA[1]*PSIB[1])
	for(t in 2:(nyear-1)){
	PSIA[t] <- PSIA[t-1] * (pmup[1]) + ((1 - PSIA[t-1]) * gyrp[t,1])
	PSIBa[t] <- PSIBa[t-1] * (pmup[2]) + ((1 - PSIBa[t-1]) * gyrp[t,2])
	PSIBA[t] <- PSIBA[t-1] * (cdp) + ((1 - PSIBA[t-1]) * (cdg[t,2]))
	PSIB[t] <- PSIA[t] * PSIBA[t] + (1 - PSIA[t]) * PSIBa[t]
	PSIAB[t] <- PSIA[t] * PSIBA[t]
	SIF[t] <- PSIAB[t]/(PSIA[t]*PSIB[t])
	}

# making these negative because it is the parameter times -1
pmup_low[1] <- ilogit(pmu[1] - pb[1])
pmup_low[2] <- ilogit(pmu[2] - pb[1])
cdp_low <- ilogit(a[2,1,2] + pmu[2] - pb[1] - a[2,2,2])
for(k in 1:nspec){
  for(t in 1:(nyear-1)){
    gyrp_low[t,k] <- ilogit(gyr[t,k] - gb[k])
    cdg_low[t,k] <- ilogit(a[2,1,1] + gyr[t,k] - gb[k] - a[2,2,1])
  }
}
PSIA_low[1] <- psinit[1] *  (pmup_low[1]) + (1 - psinit[1]) * gyrp_low[1,1]
PSIBa_low[1] <- psinit[2] * (pmup_low[2]) + (1 - psinit[2]) * gyrp_low[1,2]
PSIBA_low[1] <- psinit[2] * (cdp_low) + (1 - psinit[2]) * (cdg_low[1,2])
PSIB_low[1] <-  PSIA_low[1] * PSIBA_low[1] + (1 - PSIA_low[1]) * PSIBa_low[1]
PSIAB_low[1] <- PSIA_low[1] * PSIBA_low[1]
SIF_low[1] <- PSIAB_low[1]/(PSIA_low[1]*PSIB_low[1])
for(t in 2:(nyear-1)){
   PSIA_low[t] <-  PSIA_low[t-1] * (pmup_low[1]) + 
     ((1 - PSIA_low[t-1]) * gyrp_low[t,1])
  PSIBa_low[t] <- PSIBa_low[t-1] * (pmup_low[2]) + ((1 - PSIBa_low[t-1]) * 
      gyrp_low[t,2])
  PSIBA_low[t] <- PSIBA_low[t-1] * (cdp_low) + ((1 - PSIBA_low[t-1]) * (cdg_low[t,2]))
   PSIB_low[t] <-  PSIA_low[t] * PSIBA_low[t] + (1 - PSIA_low[t]) * PSIBa_low[t]
  PSIAB_low[t] <-  PSIA_low[t] * PSIBA_low[t]
    SIF_low[t] <- PSIAB_low[t]/(PSIA_low[t]*PSIB_low[t])
}



# now doing it for high

pmup_high[1] <- ilogit(pmu[1] + pb[1])
pmup_high[2] <- ilogit(pmu[2] + pb[1])
cdp_high <- ilogit(a[2,1,2] + pmu[2] + pb[1] + a[2,2,2])
for(k in 1:nspec){
  for(t in 1:(nyear-1)){
    gyrp_high[t,k] <- ilogit(gyr[t,k] + gb[k])
    cdg_high[t,k] <- ilogit(a[2,1,1] + gyr[t,k] + gb[k] + a[2,2,1])
  }
}
PSIA_high[1] <- psinit[1] *  (pmup_high[1]) + (1 - psinit[1]) * gyrp_high[1,1]
PSIBa_high[1] <- psinit[2] * (pmup_high[2]) + (1 - psinit[2]) * gyrp_high[1,2]
PSIBA_high[1] <- psinit[2] * (cdp_high) + (1 - psinit[2]) * (cdg_high[1,2])
PSIB_high[1] <-  PSIA_high[1] * PSIBA_high[1] + (1 - PSIA_high[1]) * PSIBa_high[1]
PSIAB_high[1] <- PSIA_high[1] * PSIBA_high[1]
SIF_high[1] <- PSIAB_high[1]/(PSIA_high[1]*PSIB_high[1])
for(t in 2:(nyear-1)){
  PSIA_high[t] <-  PSIA_high[t-1] * (pmup_high[1]) + 
    ((1 - PSIA_high[t-1]) * gyrp_high[t,1])
  PSIBa_high[t] <- PSIBa_high[t-1] * (pmup_high[2]) + ((1 - PSIBa_high[t-1]) * 
      gyrp_high[t,2])
  PSIBA_high[t] <- PSIBA_high[t-1] * (cdp_high) + 
    ((1 - PSIBA_high[t-1]) * (cdg_high[t,2]))
  PSIB_high[t] <-  PSIA_high[t] * PSIBA_high[t] + (1 - PSIA_high[t]) * PSIBa_high[t]
  PSIAB_high[t] <-  PSIA_high[t] * PSIBA_high[t]
  SIF_high[t] <- PSIAB_high[t]/(PSIA_high[t]*PSIB_high[t])
}

}
