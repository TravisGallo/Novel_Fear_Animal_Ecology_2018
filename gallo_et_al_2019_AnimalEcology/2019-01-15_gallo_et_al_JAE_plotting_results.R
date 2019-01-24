# Run analysis
source("2019-01-15_gallo_et_al_JAE_fit_3sp_JAGSmodel.R")


## Plot species interaction factors
# mmat is results from model output
plot_SIF <- function(mmat, species, y_label=TRUE, legend=FALSE){
  
  # Extract species interaction factors at high, medium, and low urbanization
  if(species == "White-tailed deer"){
    sif <- mmat[,10:19]
    sif_low <- mmat[,30:39]
    sif_hi <- mmat[,20:29]
  } else {
    sif <- mmat[,40:49]
    sif_low <- mmat[,60:69]
    sif_hi <- mmat[,50:59]
  }
  
  # Create empty plot
  plot(1:10, xlim=c(1,10), ylim=c(0,3), col="white", xaxt="n", main=species, yaxt="n", 
       xlab="", ylab="", bty='l', frame=FALSE, cex.main=3, font.lab=2)
  # x-axis
  axis(1,seq(1,10,1),labels=seq(2,11,1), mgp=c(3,0.35,0))
  mtext("Season",1, line=1.5, font=2, cex=1.1)
  # y-axis
  axis(2,seq(0,3,1), las=1)
  
  # y-axis label
  if(y_label){
    mtext("Species interaction factor",2, line=2.2, font=2, cex=1.1)
  }
  
  # Line at 1 which indicates independend spatial distributions
  abline(h=1, lty=2, lwd=1, col="grey")
  
  # Add point estimates for sif_med for each season
  for(i in 1:ncol(sif)){
  points(i,median(sif[,i]), pch=16, col="#889265")
  }
  # Add error bars
  for(i in 1:ncol(sif)){
    segments(i,quantile(sif[,i], probs=0.025),i,quantile(sif[,i], probs=0.975), 
             lwd=1, col="#889265")
  }
  
  # Add point estimates for sif_low for each season
  for(i in 1:ncol(sif)){
    points(i-0.25,median(sif_low[,i]), pch=15, col='#caf270')
  }
  # Add error bars
  for(i in 1:ncol(sif_low)){
    segments(i-0.25,quantile(sif_low[,i], probs=0.025),i-0.25,quantile(sif_low[,i], probs=0.975), 
             lwd=1, col='#caf270')
  }
  
  # Add point estimates for sif_hi for each season
  for(i in 1:ncol(sif_hi)){
    points(i+0.25,median(sif_hi[,i]), pch=17, col='#453b52')
  }
  # Add error bars
  for(i in 1:ncol(sif_hi)){
    segments(i+0.25,quantile(sif_hi[,i], probs=0.025),i+0.25,quantile(sif_hi[,i], probs=0.975), 
             lwd=1, col='#453b52')
  }
  
  # Legend
  if(legend){
    legend(7.75,3,c("Less urban","Medium urban","More urban" ), 
           col=c("#caf270","#889265","#453b52"), pch=c(15,16,17), 
           pt.cex=1.25,lwd=c(2.5,2.5,2.5), lty=c(1,1,1), bty="n", xpd=TRUE, cex=0.85)
  }
}


## Plot Vigilance

# yp, wp, co, come from pred_vig function
plot_vig <- function(yp, wp, species, co, occu, vig_results, y_label=TRUE, legend=FALSE){
  
  # Best fit line for vigilance without coyotes
  plot(yp[2,]~co, type = 'l', lwd=1.5, xlim=c(-3,3), xlab= "", xaxt="n", ylim=c(0,1), 
      ylab ="", yaxt="n", frame=FALSE)
  
  # Plot coyote occupancy along bottom
  denstrip(co,occu[2,], at = 0.045, width=0.1, colmax="black", colmin="white", 
           scale=max(occu[2,]), xpd=TRUE)
  
  # Axis after density strip so they go over the top of the denstrip
  # x-axis
  axis(1,seq(-3,3,1),line=0.05, mgp=c(3,0.35,0))
  mtext("Urbanization index",1, line=1.5, font=2, cex=1.1)
  # y-axis
  axis(2,seq(0,1,0.2), las=1)
  
  # y-axis label
  if(y_label){
    mtext("Probability of vigilance",2, line=2, font=2, cex=1.1)
  }
  
  # 95% credible intervals
  x1 <- pred_deer_vig$co
  x2 <- rev(pred_deer_vig$co)
  y1 <- yp[1,]
  y2 <- rev(yp[3,])
  polygon(c(x1, x2), c(y1, y2), col = alpha("black", .60), border = NA)
  
  # Best fit line for vigilance with coyotes
  lines(wp[2,]~co, col="#F6274C", lwd = 1.5)
  # 95% credible intervals
  y1 <- wp[1,]
  y2 <- rev(wp[3,])
  polygon(c(x1, x2), c(y1, y2), col = alpha("#F6274C", .60), border = NA)
  
  # Mean data w/o coyotes
  points(vig_results[which(vig_results$coyt == 0 & vig_results$prop > 0 & vig_results$prop < 1),"urb"] + 0.05, 
         vig_results[which(vig_results$coyt == 0 & vig_results$prop > 0 & vig_results$prop < 1),"prop"],
         pch=16, cex=0.9)
  # Mean data points w/ coyotes
  points(vig_results[which(vig_results$coyt == 1 & vig_results$prop > 0 & vig_results$prop < 1),"urb"] - 0.05, 
         vig_results[which(vig_results$coyt == 1 & vig_results$prop > 0 & vig_results$prop < 1),"prop"],
         pch=16, col="#F6274C", cex=0.9)
  
  # Legend
  if(legend){
    denstrip.legend(2,0.95,colmax="black", colmin="white", horiz=TRUE, main="", 
                    width= 0.025, len=2, nticks=1)
    legend(0.8,1.2,c("Coyotes present", "Coyotes absent"), col=c("#F6274C","black"), 
           lwd=3, lty=1, bty="n", xpd=TRUE, cex=0.85)
  }
}


## Plot overlap estimates

# Density and estimate come from overlap_extract function
plot_overlap <- function(density, estimate, species, coyt_den, y_label=TRUE, legend=FALSE){
  
  # Red line is with coyotes
  # Black line is without coyotes
  
  # Plot the first MCMC iteration of daily activity when coyotes are absent
  plot(density[[1]][,1],density[[1]][,2], type="l", lwd=1, xaxt="n", yaxt="n", xlab="", 
       ylim=c(0,0.13), ylab="", bty='l', frame=FALSE, col=alpha("black",0.01))
  
  # x-axis
  axis(1,seq(0,24,3), labels=c("0:00","3:00","6:00","9:00","12:00","15:00","18:00","21:00","24:00"), 
       mgp=c(3,0.35,0))
  mtext("Time of day",1, line=1.5, font=2, cex=1.1)
  # y-axis
  axis(2,seq(0.00,0.13,0.01), las=1)
  
  # y-axis label
  if(y_label){
    mtext("Density of activity",2, line=2.2, font=2, cex=1.1)
  }
  
  # Dispaly the density estmate and 95% credible intervals
  text(12,0.13, substitute(paste("", hat(Delta)[4], "=", med, " (95% BCI ", lo,"-",hi,")", sep=""), 
                           list(med=sprintf("%.2f", round(median(estimate),digits=2)),
                                lo=sprintf("%.2f",round(quantile(estimate, probs=0.025), digits=2)), 
                                hi=sprintf("%.2f",round(quantile(estimate,probs=0.975), digits=2)))), 
       cex=1.5, xpd=TRUE)
  
  # Draw the overalap polygon for first MCMC iteration
  polygon(c(0,density[[1]][,1],24),c(0,density[[1]][,4],0),
          col=alpha("grey89",0.01), border=NA)
  
  # Draw first MCMC iteraction of daily activity when coyotes are present
  lines(density[[1]][,1],density[[1]][,3], col=alpha("#F6274C", 0.01), lwd=1)
  
  # Loop through to draw all of the overlap polygons for each MCMC interation
  for (i in 2:length(density)){
    polygon(c(0,density[[i]][,1],24),c(0,density[[i]][,4],0),col="grey89", border=NA)
  }
  # Loop through to draw all lines for each MCMC interation
  for (i in 2:length(density)){
    # Without coyote
    lines(density[[i]][,1],density[[i]][,2], col=alpha("black", 0.01), lwd=1)
    # With coyote
    lines(density[[i]][,1],density[[i]][,3], col=alpha("#F6274C", 0.01), lwd=1)
  }
  
  # Plot overall activity of coyote
  densityPlot(coyt_den, lty=3, lwd=2, col="darkgrey", add=TRUE, extend=NULL)
  
  # Legend
  if(legend){
    legend(17.2,0.13,c("Coyotes present", "Coyotes absent", "Coyote activity","Activity overlap"), 
           col=c("#F6274C","black","darkgrey","grey89"), pch=c(NA,NA,NA,15), pt.cex=2, 
           lwd=c(3,3,1,NA), lty=c(1,1,3,0), bty="n", xpd=TRUE, cex=0.85)
  }
}

## Plot all results into a 3x2 panel

# Create layout
m <- matrix(seq(1,6,1),nrow=3, ncol=2, byrow=TRUE)

# Call device
tiff("gallo_et_al_JAE_Figure2.tiff",width=173, height=180, units="mm", res=600, compression="lzw")
# Call layout
layout(m)
# Set parameters
par(mar=c(4,3.6,2,0) + 0.1)
par(ps=7)

# Plot species interaction factors
plot_SIF(mod_3sp$mmat,"White-tailed deer", y_label=TRUE, legend=FALSE)
plot_SIF(mod_3sp$mmat, "Eastern cottontail", y_label=FALSE, legend=TRUE)
# Plot overlap results
plot_overlap(deer_overlap$overlap_dens, deer_overlap$overlap_est, 
                  "White-tailed deer", deer_overlap$coyt_rads, y_label=TRUE, legend=FALSE)
plot_overlap(rabbit_overlap$overlap_dens, rabbit_overlap$overlap_est, 
                    "Eastern cottontail", deer_overlap$coyt_rads, y_label=FALSE, legend=TRUE)
# Plot vigilance
plot_vig(pred_deer_vig$yp, pred_deer_vig$wp, "White-tailed Deer",pred_deer_vig$co, 
                 pred_deer_vig$occu, vig_deer, y_label=TRUE, legend=FALSE)
plot_vig(pred_rabbit_vig$yp, pred_rabbit_vig$wp, "Eastern Cottontail",pred_rabbit_vig$co, 
                 pred_rabbit_vig$occu, vig_rab, y_label=FALSE, legend=TRUE)

dev.off()
