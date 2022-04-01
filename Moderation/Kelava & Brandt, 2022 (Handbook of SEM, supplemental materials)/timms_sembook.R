# timms usa 8th grader: example for sembook

dat4 <- read.table("sample_timms2015.dat",header=T)
N <- nrow(dat4)


# table for descriptives
library(psych)


describe(dat4)

# read functions for standardization
source("standardization.R")


#################################################
#################################################
# Mplus via Mplus
#################################################
#################################################
library(MplusAutomation)

# save data for Mplus without colnames
write.table(dat4,"mplus.dat",col.names = F,row.names = F,na="999")

# run model in mplus
runModels("lms.inp")

# get output
lmsout <- readModels("lms.out")
# all parameter estimates (unstandardized)
lmsout$parameters$unstandardized
# regression coefficients
lmsreg <- lmsout$parameters$unstandardized[c(10:12),]

# standardization based on standardization.R
# get the latent variables variances and covariances
covxi.lms <- lmsout$parameters$unstandardized[c(26,27,13,37),3]
# standardization function
lmsreg$std <- std.reg.lms(lmsreg$est,covxi.lms)
  
# output lms
lmsreg


#################################################
#################################################
# NSEMM via Mplus
#################################################
#################################################
library(MplusAutomation)

# save data for Mplus without colnames
write.table(dat4,"mplus.dat",col.names = F,row.names = F,na="999")

# run model in mplus
runModels("nsemm.inp")

# get output
nsemmout <- readModels("nsemm.out")
# all parameter estimates (unstandardized)
nsemmout$parameters$unstandardized
#only regression coefficients
nsemmreg <- nsemmout$parameters$unstandardized[c(47:49),-7]

# standardization based on standardization.R
# get the latent variables variances and covariances 
covxi.nsemm <- nsemmout$parameters$unstandardized[c(26,27,13,
                                                    63,64,50,
                                                    74),3]
# get the latent predictor means E[xi]
exi.nsemm   <- nsemmout$parameters$unstandardized[c(14,15,51,52),3]
# get mean of latent class variable (to compute mixture proportions)
muc.nsemm   <- nsemmout$parameters$unstandardized[75,3]

# standardization function
nsemmreg$std <- std.reg.nsemm(nsemmreg$est,covxi.nsemm,exi.nsemm,muc.nsemm)

nsemmreg

#################################################
#################################################
# UPI via lavaan
# We only use 3 PIs here, there could be more (up to 9); then residual covariances need to be included
#################################################
#################################################
library(lavaan)

# create product indicators "exp x int = ei"
dat.upi <- dat4
dat.upi$ei1 <- dat.upi$exp1*dat.upi$intrinsic1
dat.upi$ei2 <- dat.upi$exp2*dat.upi$intrinsic2
dat.upi$ei3 <- dat.upi$exp3*dat.upi$intrinsic3

# lavaan model 
model.upi <- '
# measurement models
xi1 =~ exp1 + exp2 + exp3
xi2 =~ intrinsic1 + intrinsic2 + intrinsic3
eta =~ math1 + math2 + math3
# measurement model for latent product term
xi1xi2 =~ ei1 + ei2 + ei3
# structural model
eta ~  xi1 + xi2 + xi1xi2
'

# run model with MLR
upi1 <- sem(model.upi,dat.upi,estimator="MLR")
# get output
upiout <- summary(upi1)

# regression coefficients
upireg <- upiout$PE[13:15,]

# standardization
# get the latent variables variances and covariances 
covxi.upi <- upiout$PE[c(28,29,31,32,33,34,30),5]
# standardization under assumption of multivariate normality (i.e. model-implied variance of product term)
upireg$std1 <- std.reg.lms(upireg$est,covxi.upi[c(1,2,4,7)])
# standardization based on lavaan-estimates for latent product term
upireg$std2 <- std.reg.upi(upireg$est,covxi.upi)


#################################################
#################################################
# MIIV-2SLS via MIIVsem package
#################################################
#################################################
library(MIIVsem)

# 1) MIIVsem needs same data as upi (i.e. same product indicators) (above, lines 95-98)
# 2) MIIVsem can use the same lavaan model specification as upi (above, lines 100-110)

# this provides additional information about each estimation equation and the MIIVs
miivs(model.upi)

# run model
sls1 <- miive(model.upi,dat.upi)
sls1
#sargan test is significant for several indicators (xi1, xi1xi2). this indicates model misfit (exclusion restriction violated)

# get the regression coefficients (estimates, se, t-, and p-values)
slsreg <- data.frame(matrix(NA,3,4))
slsreg[,1] <- sls1$coefficients[c("eta~xi1","eta~xi2","eta~xi1xi2")]
slsreg[,2] <- sqrt(diag(sls1$coefCov))[c("eta~xi1","eta~xi2","eta~xi1xi2")]
slsreg[,3] <- slsreg[,1]/slsreg[,2]
slsreg[,4] <- 2*(1-pt(abs(slsreg[,3]),df=N-1))

colnames(slsreg) <- c("est","se","t","p")
rownames(slsreg) <- paste0("b",2:4)

slsreg

# standardization based on manifest data
# get the relevant variances and covariances 
covxi.sls <- c(var(dat.upi$exp1),var(dat.upi$intrinsic1),var(dat.upi$ei1),
               cov(dat.upi$exp1,dat.upi$intrinsic1),cov(dat.upi$exp1,dat.upi$ei1),cov(dat.upi$intrinsic1,dat.upi$ei1),
               var(dat4$math1))
# standardization under assumption of multivariate normality (i.e. model-implied variance of product term)
slsreg$std1 <- std.reg.lms(slsreg$est,covxi.sls[c(1,2,4,7)])
# standardization based on descriptive statistics
slsreg$std2 <- std.reg.upi(slsreg$est,covxi.sls)

slsreg

#################################################
#################################################
# 2SMM via lavaan and sample code
#################################################
#################################################
library(lavaan)

# stage 1 cfa model (order of parameters and variables in this file is relevant for stage 2 function)
model.cfa <- readLines("cfa_lavaan.lav")

# function for 2smm (both stages)
# the function uses the corrected standard errors from Wall and Amemyia 2003 (but correction is super-minimal)
source("2smm_delta_lavaan_interaction.R")

# run 2smm 
results.2smm <- run2smm(dat4)
# two objects: regression coefficients and latent variables covariance matrix for standardization
results.2smm$regression

# regression coefficients
smmreg <- data.frame(results.2smm$regression[2:4,])

# standardization based on latent variables covariance matrix
covxi.smm <- results.2smm$covxi
smmreg$std2 <- std.reg.smm(smmreg$est,covxi.smm)

smmreg

#################################################
#################################################
# Bayesian implementation via stan
#################################################
#################################################
library(rstan)
# some options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# data needs to be a list, "ones" is a technical support vector
data2 <- list(N=N,y=dat4[,1:3],x=dat4[,1:6+3],ones=rep(1,N))
# parameters to be monitored
params <- c("b1","lx","ly","tx","ty","sigmaeta","phi","sigmax","sigmay")

# read model code from external file
rt1 <- stanc("standard_bayes1.stan")
# compile
sm1 <- stan_model(stanc_ret = rt1, verbose=FALSE)
# run stan (mcmc sampling) with standard 4k iterations and 4 chains
fit1 <- sampling(sm1, data=data2)

# optional: save output
#saveRDS(fit1,"fit1_standard_bayes.Rdata")

# optional read output from saved location above
#fit1 <- readRDS("fit1_standard_bayes.Rdata")

# convergence checks
# rhat<1.1
stan_rhat(fit1)
# trace plots for regression coefficients for chain mixing
stan_trace(fit1,"b1")
# densitiy plots to check posterior (univariate) distributions
stan_dens(fit1,"b1",separate_chains =T)

# output
print(fit1,pars=params)
print(fit1,pars="b1")

# get relevant parameters
la2  <- as.data.frame(fit1)
params2 <- c(paste0("b1[",2:4,"]"),"phi[1,1]","phi[2,2]","phi[2,1]","sigmaeta")

# mean, sd, 95% credible interval
m.est  <- apply(la2[params2],2,mean)
sd.est <- apply(la2[params2],2,sd)
low.est <- apply(la2[params2],2,sort)[.025*dim(la2)[1],]
hig.est <- apply(la2[params2],2,sort)[.975*dim(la2)[1],]

# regression coefficients
bayesreg <- data.frame(m.est,sd.est,low.est,hig.est)[1:3,]

# standardization
covxi.bay <- m.est[-c(1:3)]
covxi.bay[4] <- covxi.bay[4]^2
bayesreg$std <- std.reg.lms(bayesreg[,1],covxi.bay)

round(bayesreg,3)



####################################
# SUMMARY
####################################
lmsreg
nsemmreg
upireg
slsreg
smmreg
bayesreg



################################################################################################################################
################################################################################################################################
# plots
################################################################################################################################
################################################################################################################################

################################################################
################################################################
# 3d graphic
################################################################
################################################################
# standardized regression coefficients #note: unstandardized version is commented out
b1 <- lmsreg[,7] #lmsreg[,3]
# standard deviations of latent factors
xr <- 1 #sqrt(covxi.lms[1])
yr <- 1 #sqrt(covxi.lms[2])

# interaction function
f2<- function(zn1,zn2){b1[1]*zn1+b1[2]*zn2+b1[3]*zn1*zn2}

# 3d plot
curve_3d <- function(f2, x_range=c(-3,3)*xr,y_range=c(-3,3)*yr, 
                     col=1:6){ 
  if (!require(rgl) ) {stop("load rgl")}
  xvec <- seq(x_range[1], x_range[2], len=15)
  yvec <- seq(y_range[1], y_range[2], len=15)
  fz   <- outer(xvec, yvec, FUN=f2)
  open3d()
  persp3d(xvec, yvec, fz,col="grey",xlab="Expectancy",ylab="Intrinsic Value",zlab="Math",smooth=F) }

# opens interactive plot
curve_3d(f2)

# saves a copy as png
rgl.snapshot("3dplot.png")

################################################################
################################################################
# alternative illustration
################################################################
################################################################
# setup grid
Nd <- 10
x <- seq(-3*xr, 3*xr, length.out = Nd)
y <- seq(-3*yr, 3*yr, length.out = Nd)
z <- outer(x, y, f2)

# open jpeg file 
jpeg("3dplot2.jpg",width=1000,height=1000)

# make the plot
persp(x = x,
      y = y,
      z = z, xlim = range(x), ylim = range(y),
      zlim = range(z, na.rm = TRUE),
      xlab = "Expectancy", ylab = "Intrinsic Value", zlab = "Math",
      #main = NULL, sub = NULL,
      theta =-40, phi = 35, r = sqrt(3), d = 1,
      scale = TRUE, expand = 1,
      col = "white", border = NULL, ltheta = -135, lphi = 0,
      shade = NA, box = TRUE, axes = TRUE, nticks = 5,
      ticktype = "simple",lwd=2,cex.lab=2)

#close the jpeg file
dev.off()


################################################################
################################################################
# simple slopes
################################################################
################################################################

# conditional regression function
slope.un <- function(x,b1,val){
  y <- b1[1]*val+b1[2]*x+b1[3]*x*val
  y
}

################################################################
#unstandardized
################################################################
b1 <- lmsreg[,3]
xr <- sqrt(covxi.lms[1])
yr <- sqrt(covxi.lms[2])

jpeg("simple_slope_lms_ustd1.jpg",width=480,height=480)
plot(0,0,col="white",xlim=c(-3,3)*yr,ylim=c(-2,1.3),ylab="Predicted Math",xlab="Intrinsic",axes=F)
whereline <- c(-2,-1,0,1,2)
for(k in 1:5){
  curve(slope.un(x=x,val=whereline[k]*xr,b1=b1),add=T,lty=k)
}

legend("topleft",paste0(c("-2SD","-1SD","+average","+1SD","+2SD")," Expectancy"),lty=1:5,bty="n")
axis(1);axis(2)
dev.off()

jpeg("simple_slope_lms_ustd2.jpg",width=480,height=480)
plot(0,0,col="white",xlim=c(-3,3)*xr,ylim=c(-2.5,2.5),ylab="Predicted Math",xlab="Expectancy",axes=F)
whereline <- c(-2,-1,0,1,2)
for(k in 1:5){
  curve(slope.un(x=x,val=whereline[k]*yr,b1=b1[c(2,1,3)]),add=T,lty=k)
}

legend("topleft",paste0(c("-2SD","-1SD","+average","+1SD","+2SD")," Intrinsic"),lty=1:5,bty="n")
axis(1);axis(2)
dev.off()

################################################################
#standardized
################################################################
b1 <- lmsreg$std
xr <- 1
yr <- 1

jpeg("simple_slope_lms_std1.jpg",width=1000,height=1000)
plot(0,0,col="white",xlim=c(-3,3)*yr,ylim=c(-3,3),ylab="Predicted Math",xlab="Intrinsic",axes=F)
whereline <- c(-2,-1,0,1,2)
for(k in 1:5){
  curve(slope.un(x=x,val=whereline[k]*xr,b1=b1),add=T,lty=k)
}

legend("topleft",paste0("Expectancy=",-2:2),lty=1:5,bty="n")
axis(1);axis(2)
dev.off()


jpeg("simple_slope_lms_std2.jpg",width=1000,height=1000)
plot(0,0,col="white",xlim=c(-3,3)*xr,ylim=c(-3,2),ylab="Predicted Math",xlab="Expectancy",axes=F)
whereline <- c(-2,-1,0,1,2)
for(k in 1:5){
  curve(slope.un(x=x,val=whereline[k]*yr,b1=b1[c(2,1,3)]),add=T,lty=k)
}

legend("topleft",paste0("Intrinsic=",c("-2","-1","0","+1","+2")),lty=1:5,bty="n")
axis(1);axis(2)
dev.off()



















