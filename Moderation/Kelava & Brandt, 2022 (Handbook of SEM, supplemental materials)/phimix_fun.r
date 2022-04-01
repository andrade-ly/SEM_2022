#################################################################
phimix.fun <- function(kap.g1,kap.g2,phi.g1,phi.g2,P){
#phimix.fun <- function(mu1g,mu2g,phi11g,phi12g,phi22g,P){
  #function calculates the covariance matrix of a mixture variable xi=(xi1,xi2)'
  #ARGUMENTS:
  #kap: matrix of class specific means 
  #phi1,phi2: 2 class specific covariance matrices (2x2) 
  #P: vector of probabilities for class membership (posterior probs)
  #RETURNS: the covariance matrix phi (5x5)
    
  # Krimskrams
  #mu1g <- c(0,1)   # das ist x1 für g==1,2
  #mu2g <- c(-1,2)  # das ist x2 für g==1,2
  #phi11g <- c(1,1)
  #phi12g <- c(.5,-.5)
  #phi22g <- c(1,1)
  #P <- c(.4,.6)
  #kapmix <- matrix(kapmix.fun(rbind(mu1g,mu2g),P),2,1)

  ######################################
  #Umsortierung
  ######################################
  mu1g <- c(kap.g1[1],kap.g2[1])
  mu2g <- c(kap.g1[2],kap.g2[2])
  
  phi11g <- c(phi.g1[1,1],phi.g2[1,1])
  phi12g <- c(phi.g1[1,2],phi.g2[1,2])
  phi22g <- c(phi.g1[2,2],phi.g2[2,2])
  
  ######################################
  #1st mixed moments (mu)
  ######################################
  mu1 <- sum(P*mu1g)
  mu2 <- sum(P*mu2g)
  
  #apply(dat,2,mean) korrekt
  ######################################
  #2nd moments
  ######################################
  
  #####################
  # 1. group specific central moments (nug)
  #####################
  nu11g <- phi11g
  nu12g <- phi12g
  nu22g <- phi22g
  
  #####################
  # 2. group specific non-central moments (mug)
  #####################
  mu11g <- mu1g^2+nu11g
  mu12g <- mu1g*mu2g+nu12g
  mu22g <- mu2g^2+nu22g
  
  #c(mean(xg1[,1]^2),mean(xg2[,1]^2))
  #c(mean(xg1[,1]*xg1[,2]),mean(xg2[,1]*xg2[,2]))
  #c(mean(xg1[,2]^2),mean(xg2[,2]^2)) #korrekt
  #####################
  #3. mixed non-central moments (mu)
  #####################
  mu11 <- sum(P*mu11g)
  mu12 <- sum(P*mu12g)
  mu22 <- sum(P*mu22g)
  
  #mean(dat[,1]^2)
  #mean(dat[,1]*dat[,2])
  #mean(dat[,2]^2) #korrekt
  #####################
  #2nd mixed central moments (nu)
  #####################
  nu11 <- mu11-mu1^2
  nu12 <- mu12-mu1*mu2
  nu22 <- mu22-mu2^2
  
  phi11 <- nu11
  phi12 <- nu12
  phi22 <- nu22
  
  #cov(dat[,1:2]) #korrekt
  ######################################
  #3rd moments
  ######################################
  
  #####################
  #1. group specific central moments (nug)
  #####################
  nu111g <- nu112g <- nu122g <- nu222g <- c(0,0)
  #nu111g[1] <- nu111g[2] <- 0
  #nu112g[1] <- nu112g[2] <- 0
  #nu122g[1] <- nu122g[2] <- 0
  #nu222g[1] <- nu222g[2] <- 0
  
  #zmoms3(xg1) #korrekt
  #zmoms3(xg2) #korrekt
  #####################
  #2. group specific non-central moments (mug)
  #####################
  mu111g <- mu1g^3+3*mu1g*phi11g
  mu222g <- mu2g^3+3*mu2g*phi22g
  mu112g <- mu1g^2*mu2g+mu2g*phi11g+2*mu1g*phi12g
  mu122g <- mu1g*mu2g^2+mu1g*phi22g+2*mu2g*phi12g
  
  #moms3(x1) #korrekt
  #moms3(x2) #korrekt
  #####################
  #3. mixed non-central moments (mu) (2 components)
  #####################
  mu111 <- sum(P*mu111g)
  mu222 <- sum(P*mu222g)
  mu112 <- sum(P*mu112g)
  mu122 <- sum(P*mu122g)
  
  #moms3(dat[,1:2]) #korrekt
  #####################
  #4. mixed central moments (nu) (2 components)
  #####################
  nu111 <- mu111-mu1^3-3*mu1*phi11
  nu222 <- mu222-mu2^3-3*mu2*phi22
  nu112 <- mu112-mu1^2*mu2-mu2*phi11-2*mu1*phi12
  nu122 <- mu122-mu1*mu2^2-2*mu2*phi12-mu1*phi22
  
  #zmoms3(dat[,1:2]) #korrekt
  ######################################  
  #4th moments
  ######################################
  
  #####################
  #1. group specific central moments (nug)
  #####################
  nu1111g <- 3*nu11g^2
  nu1112g <- 3*nu11g*nu12g
  nu1122g <- nu11g*nu22g+2*nu12g^2
  nu1222g <- 3*nu22g*nu12g
  nu2222g <- 3*nu22g^2
  
  #zmoms4(xg1) #korrekt
  #zmoms4(xg2) #korrekt
  #####################
  #2. group specific non-central moments (mug)
  #####################
  mu1111g <- 3*phi11g^2+mu1g^4+6*mu1g^2*phi11g
  mu1112g <- 3*phi11g*phi12g+mu1g^3*mu2g+3*mu2g*mu1g*phi11g+3*mu1g^2*phi12g 
  mu1122g <- phi11g*phi22g+2*phi12g^2+mu1g^2*mu2g^2+mu2g^2*phi11g+4*mu2g*mu1g*phi12g+mu1g^2*phi22g
  mu1222g <- 3*phi22g*phi12g+mu1g*mu2g^3+3*mu2g^2*phi12g+3*mu2g*mu1g*phi22g
  mu2222g <- 3*phi22g^2+mu2g^4+6*mu2g^2*phi22g
  
  #moms4(xg1) #korrekt
  #moms4(xg2) #korrekt
  #####################
  #3. mixed non-central moments (mu) (2 components)
  #####################
  mu1111 <- sum(P*mu1111g)
  mu1112 <- sum(P*mu1112g)
  mu1122 <- sum(P*mu1122g)
  mu1222 <- sum(P*mu1222g)
  mu2222 <- sum(P*mu2222g)
  
  #moms4(dat[,1:2]) #korrekt
  #####################
  #4. mixed central moments (nu) (2 components)
  #####################
  nu1111 <- mu1111-4*mu111*mu1+3*mu1^4+6*mu1^2*phi11
  nu1112 <- mu1112-3*mu112*mu1+3*mu1^3*mu2+3*mu1*mu2*phi11+3*mu1^2*phi12-mu111*mu2
  nu1122 <- mu1122-2*mu112*mu2+3*mu1^2*mu2^2+mu2^2*phi11-2*mu122*mu1+4*mu1*mu2*phi12+mu1^2*phi22
  nu1222 <- mu1222-3*mu122*mu2+3*mu2^3*mu1+3*mu1*mu2*phi22+3*mu2^2*phi12-mu222*mu1
  nu2222 <- mu2222-4*mu222*mu2+3*mu2^4+6*mu2^2*phi22
  
  #zmoms4(dat[,1:2]) #korrekt
  ##########################  
  #covariance matrix
  ##########################  
  phimix <- diag(5)
  phimix[1,1] <- phi11
  phimix[1,2] <- phimix[2,1] <- phi12
  phimix[2,2] <- phi22
  
  #variances of product terms
  phimix[3,3] <- -2*mu1*mu2*phi12-phi12^2+mu1122-mu1^2*mu2^2
  phimix[4,4] <- -2*mu1^2*phi11-phi11^2+mu1111-mu1^4
  phimix[5,5] <- -2*mu2^2*phi22-phi22^2+mu2222-mu2^4
  
  #covariances of product terms
  phimix[1,3] <-phimix[3,1] <- -mu1*phi12+mu112-mu1^2*mu2
  phimix[1,4] <-phimix[4,1] <- -mu1*phi11+mu111-mu1^3
  phimix[1,5] <-phimix[5,1] <- -mu1*phi22-mu1*mu2^2+mu122
  phimix[2,3] <-phimix[3,2] <- -mu2*phi12-mu1*mu2^2+mu122
  phimix[2,4] <-phimix[4,2] <- -mu2*phi11-mu1^2*mu2+mu112
  phimix[2,5] <-phimix[5,2] <- -mu2*phi22-mu2^3+mu222
  phimix[3,4] <-phimix[4,3] <- -mu1^2*phi12-mu1*mu2*phi11-phi12*phi11+mu1112-mu1^3*mu2
  phimix[3,5] <-phimix[5,3] <- -mu2^2*phi12-mu1*mu2*phi22-phi12*phi22+mu1222-mu2^3*mu1
  phimix[4,5] <-phimix[5,4] <- -phi11*phi22+mu1122-mu1^2*mu2^2-mu2^2*phi11-mu1^2*phi22
  
  #phimix;cov(dat);round(phimix-cov(dat),3)
  #ALT
  #phimix[3,3] <- mu1^2*phi22+mu2^2*phi11+2*mu1*mu2*phi12+phi11*phi22+phi12^2+
  #  mu1*nu122+2*mu2*nu112+nu1122
  #phimix[4,4] <- 4*mu1^2*phi11+2*phi11^2+4*mu1*nu111+nu1111
  #phimix[5,5] <- 4*mu2^2*phi22+2*phi22^2+4*mu2*nu222+nu2222
  #phimix[1,3] <-phimix[3,1] <- mu1*phi12+mu2*phi11+nu112
  #phimix[1,4] <-phimix[4,1] <- 2*mu1*phi11+nu111
  #phimix[1,5] <-phimix[5,1] <- 2*mu2*phi12+nu122
  #phimix[2,3] <-phimix[3,2] <- mu1*phi22+mu2*phi12+nu122
  #phimix[2,4] <-phimix[4,2] <- 2*mu1*phi12+nu112
  #phimix[2,5] <-phimix[5,2] <- 2*mu2*phi22+nu122
  #phimix[3,4] <-phimix[4,3] <- 2*mu1^2*phi12+2*mu1*mu2*phi11+2*phi11*phi12+nu1112+
  #  3*mu1*nu112+mu2*nu111
  #phimix[3,5] <-phimix[5,3] <- 2*mu2^2*phi12+2*mu1*mu2*phi22+2*phi22*phi12+nu1222+
  #  3*mu2*nu122+mu1*nu222
  #phimix[4,5] <-phimix[5,4] <- 4*mu1*mu2*phi12+2*phi12^2+nu1122+2*mu1*nu122+2*mu2*nu112
  
  return(phimix)
}
