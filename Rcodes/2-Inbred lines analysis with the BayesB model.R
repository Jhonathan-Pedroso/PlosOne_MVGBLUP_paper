############################################################################################################
###########################Genomic inbreds lines analysis using BayesB method###############################

#Function to build the incidence matrix of the allele substitution effect 
Wa.vit= function (K) {
  freq=colMeans(K)/2
  
  Fa=t(matrix(freq,ncol(K),nrow(K)))
  
  W=K-2*Fa
  rm(K); rm(Fa)
  return (W)
}

#Function to build the incidence matrix of the dominance deviation effect 
Wd.vit = function (K) {
  freq=colMeans(K)/2
  K=round(K,0)
  k <- K
  index2 <- k==2
  index1 <- k==1
  index0 <- k==0
  FD <- t(matrix(-2*(1-freq)^2,ncol(K),nrow(K)))
  k[index2] <- FD[index2]
  FD <- t(matrix(2*freq*(1-freq),ncol(K),nrow(K)))
  k[index1] <- FD[index1]
  FD <- t(matrix(-2*(freq)^2,ncol(K),nrow(K)))
  k[index0] <- FD[index0]
  Dv <- k
  rm(index2); rm(index1); rm(index0); rm(k); rm(K)
  return(Dv)
}
#################

##All the phenotypic and genotypic data processed used to obtain the results of the manuscript can be download here:
setwd("/home/jhonathan/Documentos/MVGBLUP/Data")
load("data_inbreds.RData")

install.packages("BGLR") #Installing BayesB package

nIter<-10000; burnIn<-1000 #Setting the number of iterations and burnIn 
library(BGLR) 

setwd("/home/jhonathan/Documentos/MVGBLUP/Results")
for (i in 1:5) {
  if (i==1) {
    y=as.matrix(data[,(i+2)]) #Separating phenotypic data from one of the five traits per loop
  }
  if (i==2) {
    y=as.matrix(data[,(i+2)])
  }
  if (i==3) {
    y=as.matrix(data[,(i+2)])
  }
  if (i==4) {
    y=as.matrix(data[,(i+2)])
  }
  if (i==5) {
    y=as.matrix(data[,(i+2)])
  }
  
  Z.loop=Z[y!="NaN",] # removing missing data
  
  ##Imputing markers
  require(rrBLUP)
  n=nrow(Z.loop); Z.imp=as.numeric(Z.loop); rm(Z.loop)
  Z.imp=matrix(Z.imp,n,ncol(Z))
  Z.imp=A.mat(Z.imp,min.MAF=F, max.missing=NULL,return.imputed=T)
  Z.imp=matrix(round(Z.imp$imputed,2),n,ncol(Z)) + 1  # Changing for the genotypic codification 2, 1 and 0
  
  Wa.loop=Wa.vit(Z.imp) #Building the incidence matrix of the allele substitution effect 
  Wd.loop=Wd.vit(Z.imp) #Building the incidence matrix of the dominance deviation effect
  
  y <- as.numeric(y[y!="NaN"]) # removing missing data
  
  ETA<-list(Ad=list(X=Wa.loop,model="BayesB",probIn=(100/27000)), #Setting the BayesB model
            Dom=list(X=Wd.loop,model="BayesB",probIn=(100/27000)))
  fmBB<-BGLR(y=y,ETA=ETA, nIter=nIter, burnIn=burnIn,saveAt=F)  #Running the BayesB analysis
  
  rm(Wa.loop, Wd.loop, Z.imp, ETA)
  if (i==1) {
    save(fmBB, file = "fmBB_PH.RData") #Saving the BayesB from the analysis of each trait per loop
  }
  if (i==2) {
    save(fmBB, file = "fmBB_EH.RData")
  }
  if (i==3) {
    save(fmBB, file = "fmBB_EL.RData")
  }
  if (i==4) {
    save(fmBB, file = "fmBB_ERN.RData")
  }
  if (i==5) {
    save(fmBB, file = "fmBB_KW.RData")
  }
}

