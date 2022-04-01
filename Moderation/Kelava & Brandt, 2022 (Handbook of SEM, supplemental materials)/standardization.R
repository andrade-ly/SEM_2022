# standardization functions


#################################################
# standardization function for lms
# model.based explained variance (standard normal variables with zero mean)
#################################################
std.reg.lms <- function(b1,covxi){
  #covxi<- covxi.lms
  covxi1 <- diag(3)
  covxi1[1,1] <- covxi[1]
  covxi1[2,2] <- covxi[2]
  covxi1[2,1] <- covxi1[1,2] <- covxi[3]
  covxi1[3,3] <- covxi[1]*covxi[2]+covxi[3]^2
  veta <- c(t(b1)%*%covxi1%*%b1+covxi[4])
  
  # standardized coefficients
  b2 <- b1
  b2[1:2] <- b1[1:2]*sqrt(covxi[1:2]/veta)
  b2[3]   <- b1[3]*sqrt(covxi[1]*covxi[2]/veta)
  b2
}

#################################################
# standardization function for upi based on estimates for 
# variances and covariances for latent product term
#################################################
std.reg.upi <- function(b1,covxi){
  
  covxi1 <- diag(3)
  covxi1[1,1] <- covxi[1]
  covxi1[2,2] <- covxi[2]
  covxi1[3,3] <- covxi[3]
  covxi1[2,1] <- covxi1[1,2] <- covxi[4]
  covxi1[3,1] <- covxi1[1,3] <- covxi[5]
  covxi1[3,2] <- covxi1[2,3] <- covxi[6]
  
  veta <- c(t(b1)%*%covxi1%*%b1+covxi[7])
  
  # standardized coefficients
  b2 <- b1
  b2[1:2] <- b1[1:2]*sqrt(covxi[1:2]/veta)
  b2[3]   <- b1[3]*sqrt(covxi[1]*covxi[2]/veta)
  b2
}

#################################################
# standardization function for 2smm based on factor score estimates
# (there is no residual variance for eta)
#################################################

std.reg.smm <- function(b1,covxi){
  
  veta <- covxi[4,4]
  covxi1 <- diag(covxi)[1:2]
  # standardized coefficients
  b2 <- b1
  b2[1:2] <- b1[1:2]*sqrt(covxi1[1:2]/veta)
  b2[3]   <- b1[3]*sqrt(covxi[1]*covxi1[2]/veta)
  b2
}




#################################################
# standardization function for nsemm (see Brandt et al., 2016)
#################################################
source("phimix_fun.R")

std.reg.nsemm <- function(b1,covxi,exi,muc){
  #covxi <- covxi.nsemm
  #exi <- exi.nsemm
  #muc <- muc.nsemm
  #b1 <- nsemmreg[,3]
  
  ga     <- c(b1,0,0) #quadratic effects are zero
  psi    <-	covxi[7]
  phi.g1 <-	matrix(covxi[c(1,3,3,2)],2,2) #Var(xi|g=1)
  phi.g2 <- matrix(covxi[c(4,6,6,5)],2,2) #Var(xi|g=1)
  
  kap.g1 <- c(exi[1:2]) #E(xi|g=1)
  kap.g2 <- c(exi[3:4]) #E(xi|g=2)
  
  
  # mean and variance of mixture distribution
  P      <- c(exp(muc),exp(0))/(exp(muc)+exp(0))  #Proportions
  kapmix <- P[1]*kap.g1+P[2]*kap.g2               #E(xi)
  phimix <- phimix.fun(kap.g1,kap.g2,phi.g1,phi.g2,P) #Var(xi)
  phi00  <- ga%*%phimix%*%ga+psi                  #Var(eta)
  
  ga.st <- c((ga[1]+ga[3]*kapmix[2]+2*ga[4]*kapmix[1])*sqrt(phimix[1,1]/phi00), #ga1
             (ga[2]+ga[3]*kapmix[1]+2*ga[5]*kapmix[2])*sqrt(phimix[2,2]/phi00), #ga2
              ga[3]*sqrt(phimix[1,1]*phimix[2,2]/phi00), #om12
              ga[4]*sqrt(phimix[1,1]*phimix[1,1]/phi00), #om11
              ga[5]*sqrt(phimix[2,2]*phimix[2,2]/phi00)) #om22
  
  ga.st[-c(4:5)]
}
