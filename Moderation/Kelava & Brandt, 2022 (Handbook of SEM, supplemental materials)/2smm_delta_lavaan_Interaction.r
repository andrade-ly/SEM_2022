################################################################################
# analyze 2smm
################################################################################

getfscores <- function(cfa,y1){
  #input: cfa results of lavaan
  #ouput: factorscores
  
  # getparameters
  beta0hat <- matrix(coef(cfa)[c(paste0("math",c(2,3),"~1"),
                                 paste0("exp",c(2,3),"~1"),
                                 paste0("intrinsic",c(2,3),"~1"))],6,1)
  
  beta1hat <- matrix(0,6,3)
  beta1hat[1:2,1] <- coef(cfa)[paste0("eta=~math",c(2,3))]
  beta1hat[3:4,2] <- coef(cfa)[paste0("xi1=~exp",c(2,3))]
  beta1hat[5:6,3] <- coef(cfa)[paste0("xi2=~intrinsic",c(2,3))]
  
  #lx[9:13,1:5] <- diag(5)
  
  psihat <- diag(9)
  diag(psihat) <- coef(cfa)[c(paste0("math",c(2,3),"~~math",c(2,3)),
                              paste0("exp",c(2,3),"~~exp",c(2,3)),
                              paste0("intrinsic",c(2,3),"~~intrinsic",c(2,3)),
                              paste0("math",1,"~~math",1),
                              paste0("exp",1,"~~exp",1),
                              paste0("intrinsic",1,"~~intrinsic",1))]
  
  b1.mod1 <- cbind(matrix(c(0),3,6),diag(3))
  b1.mod2 <- rbind(diag(6),-1*t(beta1hat))
  
  # eq(6)
  Gamma <- b1.mod1 %*% psihat %*% b1.mod2 %*% solve(t(b1.mod2) %*% psihat %*% (b1.mod2))
  sigmaee <- cbind(-Gamma,(diag(3)+ Gamma %*% beta1hat)) %*% psihat %*%  t(b1.mod1)
  
  # eq(5)
  fs<-matrix(c(0),N,3)
  for(i in 1:N){
    fs[i,] <- cbind(-Gamma,(diag(3)+ Gamma %*% beta1hat))%*%(t(as.matrix(y1[i,]))-rbind(beta0hat,matrix(0,3,1)))
  }
  
  colnames(fs) <- c("math","exp","intrinsic")
  out <- list(fs,sigmaee,Gamma,psihat,beta1hat,beta0hat)
  names(out) <- c("fs","sigmaee","Gamma","psihat","beta1hat","beta0hat")
  out
}



run2smm <- function(dat){
  #dat <- dat4
  alpha<-matrix(0,1,4)
  se <-matrix(0,1,4)
  
  
  dat <- dat[,c(2,3,5,6,8,9,1,4,7)]
  n <- nrow(dat)
  
  cfa <- sem(model.cfa, data=dat, estimator="ML")
  
  fs0 <- getfscores(cfa,dat)
  fs      <- fs0$fs
  sigmaee <- fs0$sigmaee
  # for standard error corrections
  Gamma <- fs0$Gamma
  psihat <- fs0$psihat
  beta1hat <- fs0$beta1hat
  beta0hat <- fs0$beta0hat
  
  f1 <- fs[,1]
  f2 <- fs[,2]
  f3 <- fs[,3]
  
  # Annahme des simple models (Scenario 3)
  # 1 (iv)
  mu123 <-0  #e1*e2*e3   mu13
  mu223 <-0  #e2^2*e3    mu23
  mu233 <-0  #e2*e3^2    mu33
  
  nu2233<-sigmaee[2,2]*sigmaee[3,3] #e2^2*e3^2 mu14
  
  # 2 (ii)
  
  # Elemente von M^
  ai11<-rep(1,n)
  ai21<-f2
  ai31<-f3
  ai41<-f2*f3-sigmaee[2,3] 
  ai22<-f2^2-sigmaee[2,2] 
  ai32<-f2*f3-sigmaee[2,3]
  ai42<-f2^2*f3-2*sigmaee[2,3]*f2-sigmaee[2,2]*f3-mu223
  ai33<-f3^2-sigmaee[3,3] 
  ai43<-f2*f3^2-2*sigmaee[2,3]*f3-sigmaee[3,3]*f2-mu233
  ai44<-f2^2*f3^2-sigmaee[3,3]*f2^2-sigmaee[2,2]*f3^2-4*sigmaee[2,3]*f2*f3+2*sigmaee[2,2]*sigmaee[3,3]+4*sigmaee[2,3]^2-2*mu233*f2-2*mu223*f3-nu2233
  
  # Elemente von m^
  di1<-f1
  di2<-f1*f2-sigmaee[1,2]
  di3<-f1*f3-sigmaee[1,3]
  di4<-f1*f2*f3-sigmaee[1,2]*f3-sigmaee[1,3]*f2-sigmaee[2,3]*f1-mu123
  
  Mhat<-matrix(1/n*c(sum(ai11),sum(ai21),sum(ai31),sum(ai41),
                     sum(ai21),sum(ai22),sum(ai32),sum(ai42),
                     sum(ai31),sum(ai32),sum(ai33),sum(ai43),
                     sum(ai41),sum(ai42),sum(ai43),sum(ai44)),4,4)
  
  mhat<-matrix(1/n*c(sum(di1),sum(di2),sum(di3),sum(di4)),4,1)
  
  # alpha^ = c(alpha1,gamma1,gamma2,omega12)
  # 2 (iii)
  alpha[1,]<-solve(Mhat)%*%mhat
  
  ##########################################################################################################################################################
  #Hier kommen dann die SE
  ##########################################################################################################################################################
  
  omegastar<-matrix(0,4,4)
  
  for(i in 1:n){
    Mihat<- matrix(c(ai11[i],ai21[i],ai31[i],ai41[i],
                     ai21[i],ai22[i],ai32[i],ai42[i],
                     ai31[i],ai32[i],ai33[i],ai43[i],
                     ai41[i],ai42[i],ai43[i],ai44[i]),4,4)
    
    mihat<- matrix(c(di1[i],di2[i],di3[i],di4[i]),4,1)
    la   <- mihat-Mihat%*%alpha[1,]
    omegastar.i <- la%*%t(la)
    omegastar  <-omegastar+omegastar.i
  }
  
  #i<-1
  
  omegastar <- 1/n*omegastar
  
  ################################################################################
  # Einschub aus Korrektur SE
  ################################################################################
  
  Chat <- matrix(0,4,48)
  
  ###########################################################################################################################
  # dl[i](a)/df1[i][i] ######################################################################################################
  # Die Vektoren stammen aus Maple ##########################################################################################
  ###########################################################################################################################
  y1 <- dat$math1
  x1 <- dat$exp1
  x4 <- dat$intrinsic1
  
  for(i in 1:n){#i<-1
    
    dla.df1 <- matrix(c(1,f2[i],f3[i],f2[i]*f3[i]-sigmaee[2,3]),4,1)
    
    dla.df2 <- matrix(c(-alpha[1,2]-f3[i]*alpha[1,4],f1[i]-alpha[1,1]-2*f2[i]*alpha[1,2]-f3[i]*alpha[1,3]-(2*f2[i]*f3[i]-2*sigmaee[2,3])
                        *alpha[1,4],-f3[i]*alpha[1,2]-(f3[i]^2-sigmaee[3,3])*alpha[1,4],f1[i]*f3[i]-sigmaee[1,3]-f3[i]*alpha[1,1]-(2*f2[i]*f3[i]-2*sigmaee[2,3])
                        *alpha[1,2]-(f3[i]^2-sigmaee[3,3])*alpha[1,3]-(2*f2[i]*f3[i]^2-2*sigmaee[3,3]*f2[i]-4*sigmaee[2,3]*f3[i]-2*mu233)*alpha[1,4]),4,1)
    
    dla.df3 <- matrix(c(-alpha[1,3]-f2[i]*alpha[1,4],-f2[i]*alpha[1,3]-(f2[i]^2-sigmaee[2,2])*alpha[1,4],f1[i]-alpha[1,1]-f2[i]*alpha[1,2]-2*
                          f3[i]*alpha[1,3]-(2*f2[i]*f3[i]-2*sigmaee[2,3])*alpha[1,4],f1[i]*f2[i]-sigmaee[1,2]-f2[i]*alpha[1,1]-(f2[i]^2-sigmaee[2,2])*alpha[1,2]-
                          (2*f2[i]*f3[i]-2*sigmaee[2,3])*alpha[1,3]-(2*f2[i]^2*f3[i]-2*sigmaee[2,2]*f3[i]-4*sigmaee[2,3]*f2[i]-2*mu223)*alpha[1,4]),4,1)
    
    ###########################################################################################################################
    ###########################################################################################################################
    
    dla.df <- cbind(dla.df1,dla.df2,dla.df3)   # 4x3 Matrix
    
    ###########################################################################################################################
    # dl[i](a)/dsigmaee #######################################################################################################
    ###########################################################################################################################
    
    dla.ds11 <- matrix(c(0,0,0,0),4,1)
    
    dla.ds12 <- matrix(c(0,-1,0,-f3[i]),4,1)
    
    dla.ds13 <- matrix(c(0,0,-1,-f2[i]),4,1)
    
    dla.ds22 <- matrix(c(0,alpha[1,2]+f3[i]*alpha[1,4],0,f3[i]*alpha[1,2]+alpha[1,4]*f3[i]^2-2*alpha[1,4]*sigmaee[3,3]),4,1)
    
    dla.ds23 <- matrix(c(alpha[1,4],
                         alpha[1,3]+2*f2[i]*alpha[1,4],
                         alpha[1,2]+2*f3[i]*alpha[1,4],
                         -f1[i]+alpha[1,1]+2*f2[i]*alpha[1,2]+2*f3[i]*alpha[1,3]+4*alpha[1,4]*f2[i]*f3[i]-8*alpha[1,4]*sigmaee[2,3]),4,1)
    
    dla.ds33 <- matrix(c(0,0,alpha[1,3]+f2[i]*alpha[1,4],f2[i]*alpha[1,3]+alpha[1,4]*f2[i]^2-2*alpha[1,4]*sigmaee[2,2]),4,1)
    
    ###########################################################################################################################
    
    dla.ds <- cbind(dla.ds11,dla.ds12,dla.ds13,dla.ds22,dla.ds23,dla.ds33)   # 4x3 Matrix
    
    ###########################################################################################################################
    # C[i] ####################################################################################################################
    ###########################################################################################################################
    
    x.var <- matrix(c(y1[i],x1[i],x4[i]),3,1)
    
    Chati<- cbind(dla.df%*%Gamma,
                  dla.df%*%Gamma%*% t(kronecker(x.var,diag(6))),
                  matrix(0,4,18),
                  dla.ds)
    
    Chat <- Chat+Chati
  }
  
  Chat <- Chat/n
  
  ###########################################################################################################################
  # Upsilon #################################################################################################################
  ###########################################################################################################################
  
  ###########################################################################################################################
  # Var(beta1) ##############################################################################################################
  ###########################################################################################################################
  
  vb1 <- matrix(0,18,18)
  
  vb1[1,1]   <- vcov(cfa)[1,1]
  vb1[2,2]   <- vcov(cfa)[2,2]
  vb1[9,9]   <- vcov(cfa)[3,3]
  vb1[10,10] <- vcov(cfa)[4,4]
  vb1[17,17] <- vcov(cfa)[5,5]
  vb1[18,18] <- vcov(cfa)[6,6]
  
  vb1[2,1]   <- vb1[1,2]   <- vcov(cfa)[1,2]
  vb1[9,1]   <- vb1[1,9]   <- vcov(cfa)[1,3]
  vb1[10,1]  <- vb1[1,10]  <- vcov(cfa)[1,4]
  vb1[17,1]  <- vb1[1,17]  <- vcov(cfa)[1,5]
  vb1[18,1]  <- vb1[1,18]  <- vcov(cfa)[1,6]
  
  vb1[9,2]   <- vb1[2,9]   <- vcov(cfa)[2,3]
  vb1[10,2]  <- vb1[2,10]  <- vcov(cfa)[2,4]
  vb1[17,2]  <- vb1[2,17]  <- vcov(cfa)[2,5]
  vb1[18,2]  <- vb1[2,18]  <- vcov(cfa)[2,6]
  
  vb1[10,9]  <- vb1[9,10]  <- vcov(cfa)[3,4]
  vb1[17,9]  <- vb1[9,17]  <- vcov(cfa)[3,5]
  vb1[18,9]  <- vb1[9,18]  <- vcov(cfa)[3,6]
  
  vb1[17,10] <- vb1[10,17] <- vcov(cfa)[4,5]
  vb1[18,10] <- vb1[10,18] <- vcov(cfa)[4,6]
  
  vb1[18,17] <- vb1[17,18] <- vcov(cfa)[5,6]
  
  # cov(beta0hat,beta1hat)      ist eigentlich auch 0
  
  cb1b0 <- matrix(0,6,18)
  
  cb1b0[,1:2]  <-  vcov(cfa)[1:6,7:8]
  cb1b0[,9:10] <-  vcov(cfa)[1:6,9:10]
  cb1b0[,17:18]<-  vcov(cfa)[1:6,11:12]
  
  # cov(psihat,beta1hat)        ist eigentlich auch 0
  cpsb1 <- matrix(0,9,18)
  
  cpsb1[,1:2]  <-  vcov(cfa)[13:21,7:8]
  cpsb1[,9:10] <-  vcov(cfa)[13:21,9:10]
  cpsb1[,17:18]<-  vcov(cfa)[13:21,11:12]
  
  ###########################################################################################################################
  # cov.theta ###############################################################################################################
  ###########################################################################################################################
  
  var.theta <- matrix(0,33,33)
  
  var.theta[1:6,1:6]     <-vcov(cfa)[1:6,1:6]           # cov(beta0hat)
  var.theta[25:33,25:33] <-vcov(cfa)[13:21,13:21]       # cov(psihat)
  var.theta[1:6,25:33]   <-vcov(cfa)[1:6,13:21]         # cov(beta0hat,psihat)  sind eigentlich 0
  var.theta[25:33,1:6]   <-vcov(cfa)[13:21,1:6]         # cov(psihat,beta0hat)
  var.theta[7:24,7:24]   <-vb1                          # cov(beta1hat)
  var.theta[1:6,7:24]    <-cb1b0                        # cov(beta1hat,beta0hat)
  var.theta[7:24,1:6]    <-t(cb1b0)                     # cov(beta0hat,beta1hat)
  
  var.theta[25:33,7:24]  <-cpsb1                        # cov(psihat,beta1hat)
  var.theta[7:24,25:33]  <-t(cpsb1)                     # cov(psihat,beta1hat)
  
  
  
  
  ###########################################################################################################################
  # sigmaee_vv^-1  ##########################################################################################################
  ###########################################################################################################################
  ps11 <- psihat[1:6,1:6]
  ps22 <- psihat[7:9,7:9]
  sigmavv_1 <- solve(ps11+beta1hat%*%ps22%*%t(beta1hat))
  
  # Aufgeteilte Psi-Matrix
  
  # Spalten von beta1hat
  bv1 <- beta1hat[,1]
  bv2 <- beta1hat[,2]
  bv3 <- beta1hat[,3]
  #cbind(bv1,bv2,bv3)
  
  ###########################################################################################################################
  # die vech-Funktion  ######################################################################################################
  ###########################################################################################################################
  
  vec <- function(A){c(A[1:3,1],A[2:3,2],A[3,3])}
  
  ###########################################################################################################################
  # dsigma/dpsi  ############################################################################################################
  ###########################################################################################################################
  
  # eq (24a)
  E1<-E2<-E3<-E4<-E5<-E6<-matrix(0,6,6)
  E1[1,1]<- 1
  E2[2,2]<- 1
  E3[3,3]<- 1
  E4[4,4]<- 1
  E5[5,5]<- 1
  E6[6,6]<- 1
  dsigma.dps1.1 <- ps22%*%t(beta1hat)%*%sigmavv_1%*%E1 %*% sigmavv_1%*%beta1hat%*%ps22 #sigma/ty1 = sigma[1,1]
  dsigma.dps1.2 <- ps22%*%t(beta1hat)%*%sigmavv_1%*%E2 %*% sigmavv_1%*%beta1hat%*%ps22
  dsigma.dps1.3 <- ps22%*%t(beta1hat)%*%sigmavv_1%*%E3 %*% sigmavv_1%*%beta1hat%*%ps22
  dsigma.dps1.4 <- ps22%*%t(beta1hat)%*%sigmavv_1%*%E4 %*% sigmavv_1%*%beta1hat%*%ps22
  dsigma.dps1.5 <- ps22%*%t(beta1hat)%*%sigmavv_1%*%E5 %*% sigmavv_1%*%beta1hat%*%ps22
  dsigma.dps1.6 <- ps22%*%t(beta1hat)%*%sigmavv_1%*%E6 %*% sigmavv_1%*%beta1hat%*%ps22
  
  #eq (24b)             
  # das läuft nicht, da (3x3) - (3x3)(3x6)(6x6)(3x3) -... nicht geht
  # im 2ten Summanden habe ich daher E7%*%t(beta1hat)%*%sigmavv_1%*%ps22 ausgetauscht!
  E7<-E8<-E9<-matrix(0,3,3)
  E7[1,1]<- 1
  E8[2,2]<- 1
  E9[3,3]<- 1
  dsigma.dps2.1 <- E7 - E7%*%t(beta1hat)%*%sigmavv_1%*%beta1hat%*%ps22 - ps22%*%t(beta1hat)%*%sigmavv_1%*%beta1hat%*%E7 + ps22%*%t(beta1hat)%*%sigmavv_1%*%bv1%*%t(bv1)%*%sigmavv_1%*%beta1hat%*%ps22
  dsigma.dps2.2 <- E8 - E8%*%t(beta1hat)%*%sigmavv_1%*%beta1hat%*%ps22 - ps22%*%t(beta1hat)%*%sigmavv_1%*%beta1hat%*%E8 + ps22%*%t(beta1hat)%*%sigmavv_1%*%bv2%*%t(bv2)%*%sigmavv_1%*%beta1hat%*%ps22
  dsigma.dps2.3 <- E9 - E9%*%t(beta1hat)%*%sigmavv_1%*%beta1hat%*%ps22 - ps22%*%t(beta1hat)%*%sigmavv_1%*%beta1hat%*%E9 + ps22%*%t(beta1hat)%*%sigmavv_1%*%bv3%*%t(bv3)%*%sigmavv_1%*%beta1hat%*%ps22
  
  dsigma.dps<- cbind(vec(dsigma.dps1.1),vec(dsigma.dps1.2),vec(dsigma.dps1.3),vec(dsigma.dps1.4),vec(dsigma.dps1.5),vec(dsigma.dps1.6),
                     vec(dsigma.dps2.1),vec(dsigma.dps2.2),vec(dsigma.dps2.3))
  
  
  ###########################################################################################################################
  # dsigma/dbeta1 ###########################################################################################################
  ###########################################################################################################################
  
  #eq (24c) 
  #Nummerierung von E nach beta1hatmatrix            
  E11<-E21<-E31<-E41<-E51<-E61<-matrix(0,6,3)
  E12<-E22<-E32<-E42<-E52<-E62<-matrix(0,6,3)
  E13<-E23<-E33<-E43<-E53<-E63<-matrix(0,6,3)
  
  E11[1,1]<- E21[2,1]<- E31[3,1]<- E41[4,1]<- E51[5,1]<- E61[6,1]<- 1
  E12[1,2]<- E22[2,2]<- E32[3,2]<- E42[4,2]<- E52[5,2]<- E62[6,2]<- 1
  E13[1,3]<- E23[2,3]<- E33[3,3]<- E43[4,3]<- E53[5,3]<- E63[6,3]<- 1
  
  # Berechnung aller Ableitungen. Es ergeben sich symmetrisch 3x3-Matrizen
  
  dsigma.db11 <- -ps22[1,1]*t(E11)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E11*ps22[1,1]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv1%*%t(E11[,1])+E11[,1]%*%t(bv1))%*%(ps22[1,1]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db21 <- -ps22[1,1]*t(E21)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E21*ps22[1,1]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv1%*%t(E21[,1])+E21[,1]%*%t(bv1))%*%(ps22[1,1]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db31 <- -ps22[1,1]*t(E31)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E31*ps22[1,1]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv1%*%t(E31[,1])+E31[,1]%*%t(bv1))%*%(ps22[1,1]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db41 <- -ps22[1,1]*t(E41)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E41*ps22[1,1]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv1%*%t(E41[,1])+E41[,1]%*%t(bv1))%*%(ps22[1,1]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db51 <- -ps22[1,1]*t(E51)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E51*ps22[1,1]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv1%*%t(E51[,1])+E51[,1]%*%t(bv1))%*%(ps22[1,1]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db61 <- -ps22[1,1]*t(E61)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E61*ps22[1,1]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv1%*%t(E61[,1])+E61[,1]%*%t(bv1))%*%(ps22[1,1]*sigmavv_1%*%beta1hat%*%ps22)
  
  dsigma.db12 <- -ps22[2,2]*t(E12)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E12*ps22[2,2]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv2%*%t(E12[,2])+E12[,2]%*%t(bv2))%*%(ps22[2,2]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db22 <- -ps22[2,2]*t(E22)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E22*ps22[2,2]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv2%*%t(E22[,2])+E22[,2]%*%t(bv2))%*%(ps22[2,2]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db32 <- -ps22[2,2]*t(E32)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E32*ps22[2,2]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv2%*%t(E32[,2])+E32[,2]%*%t(bv2))%*%(ps22[2,2]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db42 <- -ps22[2,2]*t(E42)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E42*ps22[2,2]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv2%*%t(E42[,2])+E42[,2]%*%t(bv2))%*%(ps22[2,2]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db52 <- -ps22[2,2]*t(E52)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E52*ps22[2,2]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv2%*%t(E52[,2])+E52[,2]%*%t(bv2))%*%(ps22[2,2]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db62 <- -ps22[2,2]*t(E62)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E62*ps22[2,2]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv2%*%t(E62[,2])+E62[,2]%*%t(bv2))%*%(ps22[2,2]*sigmavv_1%*%beta1hat%*%ps22)
  
  dsigma.db13 <- -ps22[3,3]*t(E13)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E13*ps22[3,3]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv3%*%t(E13[,3])+E13[,3]%*%t(bv3))%*%(ps22[3,3]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db23 <- -ps22[3,3]*t(E23)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E23*ps22[3,3]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv3%*%t(E23[,3])+E23[,3]%*%t(bv3))%*%(ps22[3,3]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db33 <- -ps22[3,3]*t(E33)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E33*ps22[3,3]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv3%*%t(E33[,3])+E33[,3]%*%t(bv3))%*%(ps22[3,3]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db43 <- -ps22[3,3]*t(E43)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E43*ps22[3,3]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv3%*%t(E43[,3])+E43[,3]%*%t(bv3))%*%(ps22[3,3]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db53 <- -ps22[3,3]*t(E53)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E53*ps22[3,3]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv3%*%t(E53[,3])+E53[,3]%*%t(bv3))%*%(ps22[3,3]*sigmavv_1%*%beta1hat%*%ps22)
  dsigma.db63 <- -ps22[3,3]*t(E63)%*%sigmavv_1%*%beta1hat%*%ps22-ps22%*%t(beta1hat)%*%sigmavv_1%*%E63*ps22[3,3]+ps22%*%t(beta1hat)%*%sigmavv_1%*%(bv3%*%t(E63[,3])+E63[,3]%*%t(bv3))%*%(ps22[3,3]*sigmavv_1%*%beta1hat%*%ps22)
  
  # Diese Matrix ist dann dveccSigmaee/dvecbeta1'
  dsigma.db<-cbind(vec(dsigma.db11),vec(dsigma.db21),vec(dsigma.db31),
                   vec(dsigma.db41),vec(dsigma.db51),vec(dsigma.db61),
                   vec(dsigma.db12),vec(dsigma.db22),vec(dsigma.db32),
                   vec(dsigma.db42),vec(dsigma.db52),vec(dsigma.db62),
                   vec(dsigma.db13),vec(dsigma.db23),vec(dsigma.db33),
                   vec(dsigma.db43),vec(dsigma.db53),vec(dsigma.db63))
  
  ###########################################################################################################################
  # zusammenfassen der Matrizen  ############################################################################################
  ###########################################################################################################################
  
  dsigma.theta  <- matrix(0,48,33)
  dsigma.theta[1:6,1:6]    <- 1
  dsigma.theta[7:24,7:24]  <- 1
  dsigma.theta[43:48,7:24] <-dsigma.db
  dsigma.theta[43:48,25:33]<-dsigma.dps
  
  Upsilon <- dsigma.theta %*% var.theta %*% t(dsigma.theta)
  
  
  ################################################################################
  ################################################################################
  ################################################################################
  
  
  omegahat <- omegastar + Chat %*% Upsilon %*% t(Chat)
  #eq (12)
  Vhat <- 1/n * solve(Mhat)%*% omegahat %*% solve(Mhat)
  #Vstar <- 1/n * solve(Mhat)%*% omegastar %*% solve(Mhat)
  
  se[1,] <- sqrt(diag(Vhat)) 
  #se2 <- sqrt(diag(Vstar)) 
  
  t.est <- alpha/se
  p.est <- 2*(1-pt(abs(t.est),df=n-1))
  
  results1 <- t(rbind(alpha,se,t.est,p.est))
  colnames(results1) <- c("est","se","t","p")
  rownames(results1) <- paste0("b",1:4)
  
  #covxi1 <- coef(cfa)[c("xi1~~xi1","xi2~~xi2","xi1~~xi2")]
  covxi <- cov(cbind(fs[,c(2,3)],fs[,2]*fs[,3],fs[,1]))
  
  
  results <- list(results1,covxi)
  names(results) <- c("regression","covxi")
  
  results
  
  
  
}



