model{
####################################
### PRIORS #########################
####################################
# persistance mvn URB parameters
# using scaled inverse wishart
# 1 = deer/urbanization without coyote covariate
# 2 = deer/urbanization with coyote covariate
per_urb_A[1:2] ~ dmnorm(mus, tau_per_A[,])
tau_per_A[1:2,1:2] ~ dwish(W[,], 3)
sigma_perm_A[1:2,1:2] <- inverse(tau_per_A[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_per_A[k,k.prime] <- sigma_perm_A[k,k.prime]/
      sqrt(sigma_perm_A[k,k] * sigma_perm_A[k.prime,k.prime])
  }
  sigma_per_A[k] <- sqrt(sigma_perm_A[k,k])
}
# 1 = rabbit/urbanization without coyote covariate
# 2 = rabbit/urbanization with coyote covariate
per_urb_B[1:2] ~ dmnorm(mus, tau_per_B[,])
tau_per_B[1:2,1:2] ~ dwish(W[,], 3)
sigma_perm_B[1:2,1:2] <- inverse(tau_per_B[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_per_B[k,k.prime] <- sigma_perm_B[k,k.prime]/
      sqrt(sigma_perm_B[k,k] * sigma_perm_B[k.prime,k.prime])
  }
  sigma_per_B[k] <- sqrt(sigma_perm_B[k,k])
}
##
# inxs persistence
## species x parameter x process array
a[1,1,2] <- 0 # nothing for coyote
a[1,2,2] <- 0 # nothing for coyote
a[2,1,2] ~ dt(0, 2.5, 1) # intercept deer
a[2,2,2] <- per_urb_A[2] # urb deer, from mvn
a[3,1,2] ~ dt(0, 2.5, 1) # intercept rabbit
a[3,2,2] <- per_urb_B[2] # urb rabbit, from mvn
##
# persistence covariates
##
pb[1] ~ dt(0,2.5,1) # coyote
pb[2] <- per_urb_A[1] # deer, from mvn
pb[3] <- per_urb_B[1] # rabbit, from mvn
#
# colonization mvn urb parameters
# using inverse wishart
col_urb_A[1:2] ~ dmnorm(mus, tau_col_A[,])
tau_col_A[1:2,1:2] ~ dwish(W[,], 3)
sigma_coln_A[1:2,1:2] <- inverse(tau_col_A[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_col_A[k,k.prime] <- sigma_coln_A[k,k.prime]/
      sqrt(sigma_coln_A[k,k] * sigma_coln_A[k.prime,k.prime])
  }
  sigma_col_A[k] <- sqrt(sigma_coln_A[k,k])
}
# colonization mvn urb parameters
# using inverse wishart
col_urb_B[1:2] ~ dmnorm(mus, tau_col_B[,])
tau_col_B[1:2,1:2] ~ dwish(W[,], 3)
sigma_coln_B[1:2,1:2] <- inverse(tau_col_B[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_col_B[k,k.prime] <- sigma_coln_B[k,k.prime]/
      sqrt(sigma_coln_B[k,k] * sigma_coln_B[k.prime,k.prime])
  }
  sigma_col_B[k] <- sqrt(sigma_coln_B[k,k])
}
##
# inxs colonization
##
a[1,1,1] <- 0 # nothing for coyote
a[1,2,1] <- 0 # nothing for coyote
a[2,1,1] ~ dt(0, 2.5, 1) # deer intercept
a[2,2,1] <- col_urb_A[2] # deer urb, from mvn
a[3,1,1] ~ dt(0, 2.5, 1) # rabbit intercept
a[3,2,1] <- col_urb_B[2] # rabbit urb, from mvn
##
# colonization covariates
##

  gb[1] ~ dt(0,2.5,1) # coyote
  gb[2] <- col_urb_A[1] # deer, from mvn
  gb[3] <- col_urb_B[1] # rabbit, from mvn

# detection mvn urb parameters deer
# using inverse wishart
det_urb_A[1:2] ~ dmnorm(mus, tau_det_A[,])
tau_det_A[1:2,1:2] ~ dwish(W[,], 3)
sigma_detn_A[1:2,1:2] <- inverse(tau_det_A[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_det_A[k,k.prime] <- sigma_detn_A[k,k.prime]/
      sqrt(sigma_detn_A[k,k] * sigma_detn_A[k.prime,k.prime])
  }
  sigma_det_A[k] <-sqrt(sigma_detn_A[k,k])
}
# detection
# using inverse wishart
det_urb_B[1:2] ~ dmnorm(mus, tau_det_B[,])
tau_det_B[1:2,1:2] ~ dwish(W[,], 3)
sigma_detn_B[1:2,1:2] <- inverse(tau_det_B[,])
for(k in 1:2){
  for(k.prime in 1:2){
    rho_det_B[k,k.prime] <- sigma_detn_B[k,k.prime]/
      sqrt(sigma_detn_B[k,k] * sigma_detn_B[k.prime,k.prime])
  }
  sigma_det_B[k] <-sqrt(sigma_detn_B[k,k])
}
##
# urb coyote on deer prior detection state
##
da[1,1] <- 0 # nothing for coyote
da[1,2] <- 0 # nothing for coyote
da[2,1] ~ dt(0, 2.5, 1) # deer intercept
da[2,2] <- det_urb_A[2] # deer urb, from mvn
da[3,1] ~ dt(0, 2.5, 1) # rabbit intercept
da[3,2] <- det_urb_B[2] # rabbit urb, from mvn
##
# covariates on detection
##
lpb[1] ~ dt(0,2.5,1) # coyote
lpb[2] <- det_urb_A[1] # deer from mvn
lpb[3] <- det_urb_B[2] # rabbit from mvn
# vigilance 
for(i in 1:2){ # one for each prey species
##
# intercepts on vigilance
##
v0[i] ~ dt(0, 2.5, 1)
# urbanization
vig_urb[i] ~ dt(0, 2.5, 1)
# coyote presence on prey vigilance
##
vcoy0[i] ~ dt(0, 2.5, 1)
}
# other parameter for occupancy linear predictor
for(i in 1:3){
##
# priors for intercepts
##
psi0[i] ~ dt(0, 2.5, 1) # initial occupancy
psi_urb[i] ~ dt(0, 2.5, 1) # urb on initial occupancy
pmu[i] ~ dt(0, 2.5,1)               # persistence intercept
lp[i] ~ dt(0, 2.5, 1)               # detection intercept
gmu[i] ~ dnorm(0, 0.3)              # colonization mu for ranef
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
	logit(psi[i,k,t]) <-(z[i,k,t-1] * (pmu[i] + (inprod(a[i,,2], inxscov[k,]) * z[1,k,t-1]) + pb[i]*pcov[k]))+
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
  for(k in 1:n_event_deer){
      logit(vig_prob_deer[k]) <- v0[1] + vig_urb[1] * vcov[hups_loc_deer[k,1]] + 
        vcoy0[1] * z[1,hups_loc_deer[k,1],hups_loc_deer[k,2]]
      hup_deer[k] ~ dbin(vig_prob_deer[k], phot_deer[hups_loc_deer[k,1],hups_loc_deer[k,2]])
  }
for(k2 in 1:n_event_rab){
      logit(vig_prob_rab[k2]) <- v0[2] + vig_urb[2] * vcov[hups_loc_rab[k2,1]] + 
        vcoy0[2] * z[1,hups_loc_rab[k2,1],hups_loc_rab[k2,2]]
      hup_rab[k2] ~ dbin(vig_prob_rab[k2], phot_rab[hups_loc_rab[k2,1],hups_loc_rab[k2,2]])
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

#  PSIA = Coyote
#  PSIB = Deer
# PSIAB = Deer conditional on coyote
# convert logit scale parameters to 
# probabilities
psinit[1] <- ilogit(psi0[1])
psinit[2] <- ilogit(psi0[2])
psinit[3] <- ilogit(psi0[3])
pmup[1] <- ilogit(pmu[1])
pmup[2] <- ilogit(pmu[2])
pmup[3] <- ilogit(pmu[3])
cdp[1] <- ilogit(a[2,1,2] + pmu[2])
cdp[2] <- ilogit(a[3,1,2] + pmu[3])

for(t in 1:(nyear-1)){
  for(k in 1:nspec){
    gyrp[t,k] <- ilogit(gyr[t,k])
  }
  cdg[t,1] <- ilogit(a[2,1,1] + gyr[t,2])
  cdg[t,2] <- ilogit(a[3,1,1] + gyr[t,3])
}
PSIA[1] <- psinit[1] * (pmup[1]) + (1 - psinit[1]) * gyrp[1,1]
PSIBa[1] <- psinit[2] * (pmup[2]) + (1 - psinit[2]) * gyrp[1,2]
PSICa[1] <- psinit[3] * (pmup[3]) + (1 - psinit[3]) * gyrp[1,3]
PSIBA[1] <- psinit[2] * (cdp[1]) + (1 - psinit[2]) * (cdg[1,1])
PSICA[1] <- psinit[3] * (cdp[2]) + (1 - psinit[2]) * (cdg[1,2])
PSIB[1] <- PSIA[1] * PSIBA[1] + (1 - PSIA[1]) * PSIBa[1]
PSIC[1] <- PSIA[1] * PSICA[1] + (1 - PSIA[1]) * PSICa[1]
PSIAB[1] <- PSIA[1] * PSIBA[1]
PSIAC[1] <- PSIA[1] * PSICA[1]
SIF_deer[1] <- PSIAB[1]/(PSIA[1]*PSIB[1])
SIF_rab[1] <- PSIAC[1]/(PSIA[1]*PSIC[1])
for(t in 2:(nyear-1)){
  PSIA[t] <- PSIA[t-1] * (pmup[1]) + ((1 - PSIA[t-1]) * gyrp[t,1])
  PSIBa[t] <- PSIBa[t-1] * (pmup[2]) + ((1 - PSIBa[t-1]) * gyrp[t,2])
  PSICa[t] <- PSICa[t-1] * (pmup[3]) + ((1 - PSICa[t-1]) * gyrp[t,3])
  PSIBA[t] <- PSIBA[t-1] * (cdp[1]) + ((1 - PSIBA[t-1]) * (cdg[t,1]))
  PSICA[t] <- PSICA[t-1] * (cdp[2]) + ((1 - PSICA[t-1]) * (cdg[t,2]))
  PSIB[t] <- PSIA[t] * PSIBA[t] + (1 - PSIA[t]) * PSIBa[t]
  PSIC[t] <- PSIA[t] * PSICA[t] + (1 - PSIA[t]) * PSICa[t]
  PSIAB[t] <- PSIA[t] * PSIBA[t]
  PSIAC[t] <- PSIA[t] * PSICA[t]
  SIF_deer[t] <- PSIAB[t]/(PSIA[t]*PSIB[t])
  SIF_rab[t] <- PSIAC[t]/(PSIA[t]*PSIC[t])
}

# making these negative because it is the parameter times -1
pmup_low[1] <- ilogit(pmu[1] - pb[1])
pmup_low[2] <- ilogit(pmu[2] - pb[2])
pmup_low[3] <- ilogit(pmu[3] - pb[3])
cdp_low[1] <- ilogit(a[2,1,2] + pmu[2] - pb[2] - a[2,2,2])
cdp_low[2] <- ilogit(a[3,1,2] + pmu[3] - pb[3] - a[3,2,2])

for(t in 1:(nyear-1)){
    for(k in 1:3){
    gyrp_low[t,k] <- ilogit(gyr[t,k] - gb[k])
    }
  cdg_low[t,1] <- ilogit(a[2,1,1] + gyr[t,2] - gb[2] - a[2,2,1])
  cdg_low[t,2] <- ilogit(a[3,1,1] + gyr[t,3] - gb[3] - a[3,2,1])
}
PSIA_low[1] <-  psinit[1] *  (pmup_low[1]) + (1 - psinit[1]) * gyrp_low[1,1]
PSIBa_low[1] <- psinit[2] * (pmup_low[2]) + (1 - psinit[2]) * gyrp_low[1,2]
PSICa_low[1] <- psinit[3] * (pmup_low[3]) + (1 - psinit[3]) * gyrp_low[1,3]
PSIBA_low[1] <- psinit[2] * (cdp_low[1]) + (1 - psinit[2]) * (cdg_low[1,1])
PSICA_low[1] <- psinit[3] * (cdp_low[2]) + (1 - psinit[3]) * (cdg_low[1,2])
PSIB_low[1] <-  PSIA_low[1] * PSIBA_low[1] + (1 - PSIA_low[1]) * PSIBa_low[1]
PSIC_low[1] <-  PSIA_low[1] * PSICA_low[1] + (1 - PSIA_low[1]) * PSICa_low[1]
PSIAB_low[1] <- PSIA_low[1] * PSIBA_low[1]
PSIAC_low[1] <- PSIA_low[1] * PSICA_low[1]
SIF_deer_low[1] <- PSIAB_low[1]/(PSIA_low[1]*PSIB_low[1])
SIF_rab_low[1]  <- PSIAC_low[1]/(PSIA_low[1]*PSIC_low[1])
for(t in 2:(nyear-1)){
  PSIA_low[t] <-  PSIA_low[t-1] * (pmup_low[1]) + ((1 - PSIA_low[t-1]) * gyrp_low[t,1])
  PSIBa_low[t] <- PSIBa_low[t-1] * (pmup_low[2]) + ((1 - PSIBa_low[t-1]) * gyrp_low[t,2])
  PSICa_low[t] <- PSICa_low[t-1] * (pmup_low[3]) + ((1 - PSICa_low[t-1]) * gyrp_low[t,3])
  PSIBA_low[t] <- PSIBA_low[t-1] * (cdp_low[1]) + ((1 - PSIBA_low[t-1]) * (cdg_low[t,1]))
  PSICA_low[t] <- PSICA_low[t-1] * (cdp_low[2]) + ((1 - PSICA_low[t-1]) * (cdg_low[t,2]))
  PSIB_low[t] <-  PSIA_low[t] * PSIBA_low[t] + (1 - PSIA_low[t]) * PSIBa_low[t]
  PSIC_low[t] <-  PSIA_low[t] * PSICA_low[t] + (1 - PSIA_low[t]) * PSICa_low[t]
  PSIAB_low[t] <-  PSIA_low[t] * PSIBA_low[t]
  PSIAC_low[t] <-  PSIA_low[t] * PSICA_low[t]
  SIF_deer_low[t] <- PSIAB_low[t]/(PSIA_low[t]*PSIB_low[t])
  SIF_rab_low[t] <- PSIAC_low[t]/(PSIA_low[t]*PSIC_low[t])
}



# now doing it for high
#
pmup_high[1] <- ilogit(pmu[1] + pb[1])
pmup_high[2] <- ilogit(pmu[2] + pb[2])
pmup_high[3] <- ilogit(pmu[3] + pb[3])
cdp_high[1] <- ilogit(a[2,1,2] + pmu[2] + pb[2] + a[2,2,2])
cdp_high[2] <- ilogit(a[3,1,2] + pmu[3] + pb[3] + a[3,2,2])

for(t in 1:(nyear-1)){
  for(k in 1:3){
    gyrp_high[t,k] <- ilogit(gyr[t,k] + gb[k])
  }
  cdg_high[t,1] <- ilogit(a[2,1,1] + gyr[t,2] + gb[2] + a[2,2,1])
  cdg_high[t,2] <- ilogit(a[3,1,1] + gyr[t,3] + gb[3] + a[3,2,1])
}
PSIA_high[1] <-  psinit[1] * (pmup_high[1]) + (1 - psinit[1]) * gyrp_high[1,1]
PSIBa_high[1] <- psinit[2] * (pmup_high[2]) + (1 - psinit[2]) * gyrp_high[1,2]
PSICa_high[1] <- psinit[3] * (pmup_high[3]) + (1 - psinit[3]) * gyrp_high[1,3]
PSIBA_high[1] <- psinit[2] * (cdp_high[1]) + (1 - psinit[2]) * (cdg_high[1,1])
PSICA_high[1] <- psinit[3] * (cdp_high[2]) + (1 - psinit[3]) * (cdg_high[1,2])
PSIB_high[1] <-  PSIA_high[1] * PSIBA_high[1] + (1 - PSIA_high[1]) * PSIBa_high[1]
PSIC_high[1] <-  PSIA_high[1] * PSICA_high[1] + (1 - PSIA_high[1]) * PSICa_high[1]
PSIAB_high[1] <- PSIA_high[1] * PSIBA_high[1]
PSIAC_high[1] <- PSIA_high[1] * PSICA_high[1]
SIF_deer_high[1] <- PSIAB_high[1]/(PSIA_high[1]*PSIB_high[1])
SIF_rab_high[1]  <- PSIAC_high[1]/(PSIA_high[1]*PSIC_high[1])
for(t in 2:(nyear-1)){
  PSIA_high[t] <-  PSIA_high[t-1] * (pmup_high[1]) + ((1 - PSIA_high[t-1]) * gyrp_high[t,1])
  PSIBa_high[t] <- PSIBa_high[t-1] * (pmup_high[2]) + ((1 - PSIBa_high[t-1]) * gyrp_high[t,2])
  PSICa_high[t] <- PSICa_high[t-1] * (pmup_high[3]) + ((1 - PSICa_high[t-1]) * gyrp_high[t,3])
  PSIBA_high[t] <- PSIBA_high[t-1] * (cdp_high[1]) + ((1 - PSIBA_high[t-1]) * (cdg_high[t,1]))
  PSICA_high[t] <- PSICA_high[t-1] * (cdp_high[2]) + ((1 - PSICA_high[t-1]) * (cdg_high[t,2]))
  PSIB_high[t] <-  PSIA_high[t] * PSIBA_high[t] + (1 - PSIA_high[t]) * PSIBa_high[t]
  PSIC_high[t] <-  PSIA_high[t] * PSICA_high[t] + (1 - PSIA_high[t]) * PSICa_high[t]
  PSIAB_high[t] <-  PSIA_high[t] * PSIBA_high[t]
  PSIAC_high[t] <-  PSIA_high[t] * PSICA_high[t]
  SIF_deer_high[t] <- PSIAB_high[t]/(PSIA_high[t]*PSIB_high[t])
  SIF_rab_high[t] <- PSIAC_high[t]/(PSIA_high[t]*PSIC_high[t])
}


}
