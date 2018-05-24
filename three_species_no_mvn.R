model{
####################################
### PRIORS #########################
####################################
# Interaction parameters
for(par in 1:4){
	ppar[par] ~ dt(0, 2.5, 1)
	gpar[par] ~ dt(0, 2.5, 1)
	dpar[par] ~ dt(0, 2.5, 1)
}
##
# inxs persistence
## species x parameter x process array
a[1,1,2] <- 0 # nothing for coyote
a[1,2,2] <- 0 # nothing for coyote
a[2,1,2] <- ppar[1] # intercept deer
a[2,2,2] <- ppar[2] # urb deer, from mvn
a[3,1,2] <- ppar[3] # intercept rabbit
a[3,2,2] <- ppar[4] # urb rabbit, from mvn
##
# persistence covariates
##
pb[1] ~ dt(0,2.5,1) # coyote
pb[2] ~ dt(0,2.5,1) # deer, from mvn
pb[3] ~ dt(0,2.5,1) # rabbit, from mvn
#
##
# inxs colonization
##
a[1,1,1] <- 0 # nothing for coyote
a[1,2,1] <- 0 # nothing for coyote
a[2,1,1] <- gpar[1] # deer intercept
a[2,2,1] <- gpar[2] # deer urb, from mvn
a[3,1,1] <- gpar[3] # rabbit intercept
a[3,2,1] <- gpar[4] # rabbit urb, from mvn
##
# colonization covariates
##
  gb[1] ~ dt(0, 2.5, 1) # coyote
  gb[2] ~ dt(0, 2.5, 1) # deer, from mvn
  gb[3] ~ dt(0, 2.5, 1) # rabbit, from mvn
# detection
##
# urb coyote on deer prior detection state
##
da[1,1] <- 0 # nothing for coyote
da[1,2] <- 0 # nothing for coyote
da[2,1] <- dpar[1] # deer intercept
da[2,2] <- dpar[2] # deer urb, from mvn
da[3,1] <- dpar[3] # rabbit intercept
da[3,2] <- dpar[4] # rabbit urb, from mvn
##
# covariates on detection
##
lpb[1] ~ dt(0,2.5,1) # coyote
lpb[2] ~ dt(0, 2.5, 1) # deer from mvn
lpb[3] ~ dt(0, 2.5, 1) # rabbit from mvn
# vigilance 
for(i in 1:2){ # one for each prey species
##
# intercepts on vigilance
##
v0[i] ~ dt(0, 2.5, 1)
# urbanization
vig_urb[i] ~ dt(0, 2.5, 1)
vcoy_urb[i] ~ dt(0, 2.5, 1)
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
	for(k in 1:nsite){
		for(i in 1:nspec){
	##
	# Initial occupancy, Intercept only.
	##
	logit(psi[i,k,1]) <- psi0[i] + psi_urb[i]*pcov[k]
	z[i,k,1] ~ dbern(psi[i,k,1])
		}
	for(t in 2:nyear){
	##
	# Mixture of colonization and persistence
	##
		# Linear predictor for coyote
	logit(psi[1,k,t]) <- (z[1,k,t-1] * (pmu[1] + pb[1] * pcov[k])) +
		(1 - z[1,k,t-1]) * (gyr[t-1,1] + gb[1] * pcov[k])
		# Linear predictor for deer
	logit(psi[2,k,t]) <- 
	 (z[2,k,t-1]*(1-z[1,k,t-1])*(pmu[2] + pb[2] * pcov[k])) + # phi w/o coyote
	 (z[2,k,t-1]*z[1,k,t-1]*(pmu[2] + inprod(a[2,,2], inxscov[k,]))) + #phi w/ coyote
	 ((1 - z[2,k,t-1])*(1-z[1,k,t-1])*(gyr[t-1,2] + gb[2] * pcov[k])) + #gam w/o coyote
	 ((1 - z[2,k,t-1])*z[1,k,t-1]*(gyr[t-1,2]+inprod(a[2,,1],inxscov[k,]))) # gam w/ coy
	# Linear predictor for rabbit
	logit(psi[3,k,t]) <- 
		(z[3,k,t-1]*(1-z[1,k,t-1])*(pmu[3] + pb[3] * pcov[k])) + # phi w/o coyote
		(z[3,k,t-1]*z[1,k,t-1]*(pmu[3] + inprod(a[3,,2], inxscov[k,]))) + #phi w/ coyote
		((1 - z[3,k,t-1])*(1-z[1,k,t-1])*(gyr[t-1,3] + gb[3] * pcov[k])) + #gam w/o coyote
		((1 - z[3,k,t-1])*z[1,k,t-1]*(gyr[t-1,3]+inprod(a[3,,1],inxscov[k,]))) # gam w/ coy
	##
	# Likelihood
	##
	for(i in 1:nspec){
	z[i,k,t] ~ dbern(psi[i,k,t])
		}
	}
}
####################################
### VIGILANCE ANALYSIS #############
####################################
  for(k in 1:n_event_deer){
      logit(vig_prob_deer[k]) <- v0[1] + # baseline
      	(vig_urb[1] * vcov[hups_loc_deer[k,1]] * # if no coyote
      	(1 - z[1,hups_loc_deer[k,1],hups_loc_deer[k,2]])) + 
        vcoy0[1] * z[1,hups_loc_deer[k,1],hups_loc_deer[k,2]] + # if coyote
      	(vcoy_urb[1] * vcov[hups_loc_deer[k,1]] *               #if coyote
      	z[1,hups_loc_deer[k,1],hups_loc_deer[k,2]])
      hup_deer[k] ~ dbin(vig_prob_deer[k], 
      	                 phot_deer[hups_loc_deer[k,1],hups_loc_deer[k,2]])
  }
for(k2 in 1:n_event_rab){
      logit(vig_prob_rab[k2]) <- v0[2] + # baseline
      	(vig_urb[2] * vcov[hups_loc_rab[k2,1]] * # if no coyote
      	(1 - z[1,hups_loc_rab[k2,1],hups_loc_rab[k2,2]])) + 
        vcoy0[2] * z[1,hups_loc_rab[k2,1],hups_loc_rab[k2,2]] + # if coyote
      	(vcoy_urb[2] * vcov[hups_loc_rab[k2,1]] *               # if coyote
      	 z[1,hups_loc_rab[k2,1],hups_loc_rab[k2,2]]) 
      hup_rab[k2] ~ dbin(vig_prob_rab[k2], 
      	                 phot_rab[hups_loc_rab[k2,1],hups_loc_rab[k2,2]])
  }

####################################
### OBSERVATIONAL PROCESS ##########
####################################

	for(k in 1:nsite){
		for(t in 1:nyear){
		  logit(dprob[1,k,t]) <- lp[1] +  lpb[1]*lpcov[k]
		  logit(dprob[2,k,t]) <- lp[2] + (lpb[2]*lpcov[k] *  (1 - z[1,k,t])) + 
		  	                       (inprod(da[2,], inxscov[k,]) * z[1,k,t])
		  logit(dprob[3,k,t]) <- lp[3] + (lpb[3]*lpcov[k] *  (1 - z[1,k,t])) + 
		  	                       (inprod(da[3,], inxscov[k,]) * z[1,k,t])
		  for(i in 1:nspec){
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
  cdg_low[t,1] <- ilogit(a[2,1,1] + gyr[t,2] - a[2,2,1])
  cdg_low[t,2] <- ilogit(a[3,1,1] + gyr[t,3] - a[3,2,1])
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
cdp_high[1] <- ilogit(a[2,1,2] + pmu[2] + a[2,2,2])
cdp_high[2] <- ilogit(a[3,1,2] + pmu[3] + a[3,2,2])

for(t in 1:(nyear-1)){
  for(k in 1:3){
    gyrp_high[t,k] <- ilogit(gyr[t,k] + gb[k])
  }
  cdg_high[t,1] <- ilogit(a[2,1,1] + gyr[t,2] + a[2,2,1])
  cdg_high[t,2] <- ilogit(a[3,1,1] + gyr[t,3] + a[3,2,1])
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
