############################################################################################################
####################################Analysis with the  GBLUP models#########################################

setwd("/home/jhonathan/Documentos/MVGBLUP/Data/")
load("data_hybrids.RData") #Loading the hybrids data

# Additive relationship matrix function
A.vit= function (K) {
  freq=colMeans(K)/2
  
  Fa=t(matrix(freq,ncol(K),nrow(K)))
  
  W=K-2*Fa
  cor=2*sum(freq*(1-freq))
  Av=W%*%t(W)/cor
  return (Av)
}

#Dominance relationship matrix function
D.vit= function (K) {
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
  cor=4*sum((freq*(1-freq))^2)
  Dv <- k%*%t(k)/cor
  index2=0
  index1=0
  index0=0
  k=0
  K=0
  return(Dv)
}

#Function to perform the spectral decomposition of the relationship matrices and correction of the eigenvalues
spec.dec <-function(K) {
  E=eigen(K)
  Dg= diag(E$values)
  k=0
  for (i in 1:nrow(Dg)) {
    if (Dg[i,i]<1e-4) {
      Dg[i,i]=Dg[(i-1),(i-1)]-0.01*Dg[(i-1),(i-1)]
    }
  }
  C=matrix(E$vectors,nrow(K),ncol(K))
  K = C%*%Dg%*%t(C)
}

A=A.vit(Z.h); D=D.vit(Z.h);  #Building additive and dominance relationship matrices
A=spec.dec(A); D=spec.dec(D); # Correcting eigenvalues to obtain a positive definite matrix
Ga1=chol2inv(chol(A)) # Computing the inverse of the additive relationship matrix
Gd1=chol2inv(chol(D)) # Computing the inverse of the dominance relationship matrix

##Running the univariate GBLUP-A model in all single-cross hybrids data scenarios
for (i in 1:20) {
  
  if (i==1) {  
    y=y.h.PH.0.3 #Pushing the data of the plant height trait constructed in the 0.3 heritability scenario 
    comp=comp.PH.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==2) { 
    y=y.h.EH.0.3 #Loading the data of the ear height trait in the 0.3 heritability scenario
    comp=comp.EH.0.3[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==3) {
    y=y.h.EL.0.3 #Pushing the data of the ear lenght trait in the 0.3 heritability scenario
    comp=comp.EL.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==4) {
    y=y.h.ERN.0.3 #Pushing the data of the ear row number trait in the 0.3 heritability scenario
    comp=comp.ERN.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==5) {
    y=y.h.KW.0.3 #Pushing the data of the kernel weight trait in the 0.3 heritability scenario
    comp=comp.KW.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==6) {  
    y=y.h.PH.0.5 #Pushing the data of the plant height trait constructed in the 0.5 heritability scenario 
    comp=comp.PH.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==7) { 
    y=y.h.EH.0.5 #Loading the data of the ear height trait in the 0.5 heritability scenario
    comp=comp.EH.0.5[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==8) {
    y=y.h.EL.0.5 #Pushing the data of the ear lenght trait in the 0.5 heritability scenario
    comp=comp.EL.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==9) {
    y=y.h.ERN.0.5 #Pushing the data of the ear row number trait in the 0.5 heritability scenario
    comp=comp.ERN.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==10) {
    y=y.h.KW.0.5 #Pushing the data of the kernel weight trait in the 0.5 heritability scenario
    comp=comp.KW.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==11) {  
    y=y.h.PH.0.7 #Pushing the data of the plant height trait constructed in the 0.7 heritability scenario 
    comp=comp.PH.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==12) { 
    y=y.h.EH.0.7 #Loading the data of the ear height trait in the 0.7 heritability scenario
    comp=comp.EH.0.7[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==13) {
    y=y.h.EL.0.7 #Pushing the data of the ear lenght trait in the 0.7 heritability scenario
    comp=comp.EL.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==14) {
    y=y.h.ERN.0.7 #Pushing the data of the ear row number trait in the 0.7 heritability scenario
    comp=comp.ERN.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==15) {
    y=y.h.KW.0.7 #Pushing the data of the kernel weight trait in the 0.7 heritability scenario
    comp=comp.KW.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==16) {  
    y=y.h.PH.hist #Pushing the data of the plant height trait constructed in the historical heritability scenario 
    comp=comp.PH.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==17) { 
    y=y.h.EH.hist #Loading the data of the ear height trait in the historical heritability scenario
    comp=comp.EH.hist[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==18) {
    y=y.h.EL.hist #Pushing the data of the ear lenght trait in the historical heritability scenario
    comp=comp.EL.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==19) {
    y=y.h.ERN.hist #Pushing the data of the ear row number trait in the historical heritability scenario
    comp=comp.ERN.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==20) {
    y=y.h.KW.hist #Pushing the data of the kernel weight trait in the historical heritability scenario
    comp=comp.KW.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  
  g.p=a.p+d.p; ve=as.numeric(1);  # Computing the parametric total genetic value
  va=as.numeric(1); time.new=0; time.old=0; #Setting the first guess of the additive variance (va) and setting parameters to compute the time of the iterations (time.new; time.old)
  iter=as.numeric(0) #Setting the variable to control the number of iterations of the EM algorithm
  dif=as.numeric(0)  #Setting the variable to control the difference between the variables estimated between iterations of the EM algorithm
  n=length(y) #Computing the length of the vector containing the phenotypic data
  x=matrix(1,n,1) #Computing the vector of fixed effects
  Zds=diag(n) # Computing the matrix of the random effects
  n1=ncol(x) # Computing the number of fixed effectrs
  n2=ncol(A) # Computing the number of random effects
  xx=crossprod(x,x) #Computing the X'X which will be inserted in theright-hand side of the mixed model equations
  xz=crossprod(x,Zds) #Computing the X'Z which will be inserted in the right-hand side of the mixed model equations
  zx=crossprod(Zds,x) #Computing the Z'X which will be inserted in the right-hand side of the mixed model equations
  zz=crossprod(Zds,Zds) #Computing the Z'Z which will be inserted in the right-hand side of the mixed model equations
  xy=crossprod(x,y) #Computing the X'Y which will be inserted in the left-hand side of the mixed model equations
  zy=crossprod(Zds,y) #Computing the Z'Y which will be inserted in the left-hand side of the mixed model equations
  W1=cbind(x,Zds) # Horizontal concatenation of the [X | Z] 
  maxiter=1500 #Setting the number of iterations
  
  setwd("/home/jhonathan/Documentos/MVGBLUP/Results")
  repeat{ #Beginning the iterative process
    ve<-as.numeric(ve)
    va<-as.numeric(va)
    
    Ga=Ga1*c(ve/va) #Computing the additive variance-covariance matrix
    
    ##Obtaining the inverse of the right-hand side of the mixed model equations
    C1=cbind(xx,xz)
    C2=cbind(zx,zz+Ga)
    C=chol2inv(chol(rbind(C1,C2)))
    
    ##Obtaining the right-hand side of the mixed model equations
    so=c(xy,zy)  
    so=as.matrix(so)
    
    ##Obtaining the solution of the mixed model equations (EM - maximization step)
    sol=C%*%so 
    
    beta=as.matrix(sol[1:n1]) #Pushing the fixed effects estimated from the solution of the mixed model equations
    ad=sol[(n1+1):(n1+n2)] #Pushing the additive effects estimated from the solution of the mixed model equations
    
    ##Computing the vector of residuals
    y=as.matrix(y)
    e=y-x%*%beta-Zds%*%as.matrix(ad) 
    
    ##Obtaining the new variances components from the EM equations (EM - expection step)
    ve1=(crossprod(e,e)+sum(diag(W1%*%C%*%t(W1)))*ve)/n
    va1=(crossprod(ad,Ga1)%*%ad+sum(diag(Ga1%*%(C[(n1+1):(n1+n2),(n1+1):(n1+n2)])))*ve)/n2
    
    iter=iter+1 #Computing the current iteration step
    
    dif=max(abs(c((ve1-ve),(va1-va)))) #Computing the difference between the new and the old parameters obtained in the iterative process
    
    ##Updating variance components
    va=va1
    ve=ve1
    
    ##Computing the model evaluation statistics 
    vg=va
    g=ad
    
    #Correlations between the estimated and parametric additive effects
    cora=cor(ad,a.p)
    #Correlations between the estimated and parametric total genetic effects
    corg=cor(g,g.p)
    #narrow-sense heritability
    ha=cov(ad,y)/var(y)
    #broad-sense heritability
    hg=cov(g,y)/var(y)
    #Predicted residual error sum of squares of the additive, total genetic effects, respectively.
    PRESS.a=sum(c((as.matrix(ad)-a.p)^2))
    PRESS.g=sum(c((as.matrix(g)-g.p)^2))
    
    ##Organizing the results of the current iterative process
    results=matrix(c(va,vg,
                     cora,corg,
                     ha,hg,
                     PRESS.a,PRESS.g),4,2,byrow=T)
    
    colnames(results) <-c("Addit.", "Gen")
    rownames(results) <-c("Variance","Correlation","heritability","PRESS")
    
    
    print("-----------------------Univariate GBLUP-A------------------------------------")  
    print("-----------------------------------------------------------------------------")
    print("-----------------------------------------------------------------------------")
    
    ##Setting the function to count the iteration time
    time=proc.time() 
    time.new=as.numeric(time[1])
    time_iter=(time.new-time.old)/60
    time.old=time.new
    
    ##Setting the iterative features to be printed in each iteration
    t=t(as.matrix(c(dif,iter,time_iter)))
    colnames(t) <- c("Dif:" , "Iter:", "Time/iter:")
    print(t)
    
    print("-----------------------------------------------------------------------------")
    
    print("Statistics:")
    print(results)
    
    print("-----------------------------------------------------------------------------")
    
    if ((iter==maxiter)) break
  }
  ##Saving the results of the analysis
  #Saving the results for of the analysis of the data constructed in the 0.3 heritability scenario
  if (i==1) {
    save(results,file = "Results_GBLUP-A_UV_PH_h=0.3.RData")
  }
  if (i==2) {
    save(results,file = "Results_GBLUP-A_UV_EH_h=0.3.RData")
  }
  if (i==3) {
    save(results,file = "Results_GBLUP-A_UV_EL_h=0.3.RData")
  }
  if (i==4) {
    save(results,file = "Results_GBLUP-A_UV_ERN_h=0.3.RData")
  }
  if (i==5) {
    save(results,file = "Results_GBLUP-A_UV_KW_h=0.3.RData")
  }
  #Saving the results for of the analysis of the data constructed in the 0.5 heritability scenario
  if (i==6) {
    save(results,file = "Results_GBLUP-A_UV_PH_h=0.5.RData")
  }
  if (i==7) {
    save(results,file = "Results_GBLUP-A_UV_EH_h=0.5.RData")
  }
  if (i==8) {
    save(results,file = "Results_GBLUP-A_UV_EL_h=0.5.RData")
  }
  if (i==9) {
    save(results,file = "Results_GBLUP-A_UV_ERN_h=0.5.RData")
  }
  if (i==10) {
    save(results,file = "Results_GBLUP-A_UV_KW_h=0.5.RData")
  }
  #Saving the results for of the analysis of the data constructed in the 0.7 heritability scenario
  if (i==11) {
    save(results,file = "Results_GBLUP-A_UV_PH_h=0.7.RData")
  }
  if (i==12) {
    save(results,file = "Results_GBLUP-A_UV_EH_h=0.7.RData")
  }
  if (i==13) {
    save(results,file = "Results_GBLUP-A_UV_EL_h=0.7.RData")
  }
  if (i==14) {
    save(results,file = "Results_GBLUP-A_UV_ERN_h=0.7.RData")
  }
  if (i==15) {
    save(results,file = "Results_GBLUP-A_UV_KW_h=0.7.RData")
  }
  #Saving the results for of the analysis of the data constructed in the historical heritability scenario
  if (i==16) {
    save(results,file = "Results_GBLUP-A_UV_PH_h=hist.RData")
  }
  if (i==17) {
    save(results,file = "Results_GBLUP-A_UV_EH_h=hist.RData")
  }
  if (i==18) {
    save(results,file = "Results_GBLUP-A_UV_EL_h=hist.RData")
  }
  if (i==19) {
    save(results,file = "Results_GBLUP-A_UV_ERN_h=hist.RData")
  }
  if (i==20) {
    save(results,file = "Results_GBLUP-A_UV_KW_h=hist.RData")
  }
  rm(C,C1,C2,Ga,W1,Zds,beta,e,ha,hg,results,so,sol,t,va,va1,ve,ve1,vg,x,xx,xy,xz,y,zx,zy,zz,PRESS.a,PRESS.g,a.p,ad,comp,cora,corg,d.p,
     dif,g,g.p,iter,maxiter,n,n1,n2,time,time.new,time.old,time_iter)  
}

##Running the univariate GBLUP-AD model in all single-cross hybrids data scenarios
for (i in 1:20) {
  
  if (i==1) {  
    y=y.h.PH.0.3 #Pushing the data of the plant height trait constructed in the 0.3 heritability scenario 
    comp=comp.PH.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==2) { 
    y=y.h.EH.0.3 #Loading the data of the ear height trait in the 0.3 heritability scenario
    comp=comp.EH.0.3[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==3) {
    y=y.h.EL.0.3 #Pushing the data of the ear lenght trait in the 0.3 heritability scenario
    comp=comp.EL.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==4) {
    y=y.h.ERN.0.3 #Pushing the data of the ear row number trait in the 0.3 heritability scenario
    comp=comp.ERN.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==5) {
    y=y.h.KW.0.3 #Pushing the data of the kernel weight trait in the 0.3 heritability scenario
    comp=comp.KW.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==6) {  
    y=y.h.PH.0.5 #Pushing the data of the plant height trait constructed in the 0.5 heritability scenario 
    comp=comp.PH.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==7) { 
    y=y.h.EH.0.5 #Loading the data of the ear height trait in the 0.5 heritability scenario
    comp=comp.EH.0.5[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==8) {
    y=y.h.EL.0.5 #Pushing the data of the ear lenght trait in the 0.5 heritability scenario
    comp=comp.EL.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==9) {
    y=y.h.ERN.0.5 #Pushing the data of the ear row number trait in the 0.5 heritability scenario
    comp=comp.ERN.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==10) {
    y=y.h.KW.0.5 #Pushing the data of the kernel weight trait in the 0.5 heritability scenario
    comp=comp.KW.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==11) {  
    y=y.h.PH.0.7 #Pushing the data of the plant height trait constructed in the 0.7 heritability scenario 
    comp=comp.PH.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==12) { 
    y=y.h.EH.0.7 #Loading the data of the ear height trait in the 0.7 heritability scenario
    comp=comp.EH.0.7[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==13) {
    y=y.h.EL.0.7 #Pushing the data of the ear lenght trait in the 0.7 heritability scenario
    comp=comp.EL.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==14) {
    y=y.h.ERN.0.7 #Pushing the data of the ear row number trait in the 0.7 heritability scenario
    comp=comp.ERN.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==15) {
    y=y.h.KW.0.7 #Pushing the data of the kernel weight trait in the 0.7 heritability scenario
    comp=comp.KW.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==16) {  
    y=y.h.PH.hist #Pushing the data of the plant height trait constructed in the historical heritability scenario 
    comp=comp.PH.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the plant height trait
  }
  if (i==17) { 
    y=y.h.EH.hist #Loading the data of the ear height trait in the historical heritability scenario
    comp=comp.EH.hist[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear height trait
  }
  if (i==18) {
    y=y.h.EL.hist #Pushing the data of the ear lenght trait in the historical heritability scenario
    comp=comp.EL.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear length trait
  }
  if (i==19) {
    y=y.h.ERN.hist #Pushing the data of the ear row number trait in the historical heritability scenario
    comp=comp.ERN.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
  }
  if (i==20) {
    y=y.h.KW.hist #Pushing the data of the kernel weight trait in the historical heritability scenario
    comp=comp.KW.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait 
    a.p=comp[1:400]; d.p=comp[401:800]; #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  
  g.p=a.p+d.p; ve=as.numeric(1);  # Computing the parametric total genetic value
  va=as.numeric(1); vd=as.numeric(1); time.new=0; time.old=0; #Setting the first guess of the additive (va) and dominance (vd) variances and setting parameters to compute the time of the iterations (time.new; time.old)
  iter=as.numeric(0) #Setting the variable to control the number of iteratins of the EM algorithm
  dif=as.numeric(0)  #Setting the variable to control the difference between the variables estimated between iterations of the EM algorithm
  n=length(y) #Computing the length of the observation vector
  x=matrix(1,n,1) #Computing the vector of fixed effects
  Zds=diag(n) # Computing the matrix of the random effects
  n1=ncol(x) # Computing the number of fixed effectrs
  n2=ncol(A) # Computing the number of random effects
  xx=crossprod(x,x) #Computing the X'X which will be inserted in theright-hand side of the mixed model equations
  xz=crossprod(x,Zds) #Computing the X'Z which will be inserted in the right-hand side of the mixed model equations
  zx=crossprod(Zds,x) #Computing the Z'X which will be inserted in the right-hand side of the mixed model equations
  zz=crossprod(Zds,Zds) #Computing the Z'Z which will be inserted in the right-hand side of the mixed model equations
  xy=crossprod(x,y) #Computing the X'Y which will be inserted in the left-hand side of the mixed model equations
  zy=crossprod(Zds,y) #Computing the Z'Y which will be inserted in the left-hand side of the mixed model equations
  W1=cbind(x,Zds,Zds) # Horizontal concatenation of the [X | Z | Z] 
  maxiter=1500 #Setting the number of iterations
  
  setwd("/home/jhonathan/Documentos/MVGBLUP/Results")
  repeat {
    ve<-as.numeric(ve) 
    va<-as.numeric(va) 
    vd<-as.numeric(vd)
    
    Ga=Ga1*c(ve/va) #Computing the additive variance-covariance matrix
    Gd=Gd1*c(ve/vd) #Computing the dominance variance-covariance matrix
    
    ##Obtaining the inverse of the right-hand side of the mixed model equations
    C1=cbind(xx,xz,xz)
    C2=cbind(zx,zz+Ga,zz)
    C3=cbind(zx,zz,zz+Gd)  
    C=chol2inv(chol(rbind(C1,C2,C3)))
    
    ##Obtaining the right-hand side of the mixed model equations
    so=c(xy,zy,zy) 
    so=as.matrix(so)
    
    ##Obtaining the solution of the mixed model equations (EM - maximization step)
    sol=C%*%so
    
    beta=as.matrix(sol[1:n1])  #Pushing the fixed effects estimated from the solution of the mixed model equations
    ad=sol[(n1+1):(n1+n2)] #Pushing the additive effects estimated from the solution of the mixed model equations
    dom=sol[(n1+n2+1):(n1+n2*2)] #Pushing the dominance deviation effects estimated from the solution of the mixed model equations
    
    ##Computing the vector of residuals
    y=as.matrix(y)
    e=y-x%*%beta-Zds%*%as.matrix(ad)-Zds%*%as.matrix(dom)
    
    ##Obtaining the new variances components estimated in the iterative process (EM - expection step)
    ve1=(crossprod(e,e)+sum(diag(W1%*%C%*%t(W1)))*ve)/n
    va1=(crossprod(ad,Ga1)%*%ad+sum(diag(Ga1%*%(C[(n1+1):(n1+n2),(n1+1):(n1+n2)])))*ve)/n2
    vd1=(crossprod(dom,Gd1)%*%dom+sum(diag(Gd1%*%(C[(n1+n2+1):(n1+n2*2),(n1+n2+1):(n1+n2*2)])))*ve)/n2
    
    iter=iter+1 #Computing the current iteration step
    dif=max(abs(c((ve1-ve),(va1-va),(vd1-vd)))) #Computing the difference between the new and the old parameters obtained in the iterative process
    
    ##Updating variance components
    va=va1
    ve=ve1
    vd=vd1
    
    ##Computing the model evaluation statistics 
    vg=va+vd
    g=ad+dom
    
    #Correlations between the estimated and parametric additive effects
    cora=cor(ad,a.p) 
    #Correlations between the estimated and parametric dominance deviation effects
    cord=cor(dom,d.p)
    #Correlations between the estimated and parametric total genetic effects
    corg=cor(g,g.p)
    #narrow-sense additive heritability
    ha=cov(ad,y)/var(y)
    #narrow-sense dominance heritability
    hd=cov(dom,y)/var(y)
    #broad-sense heritability
    hg=cov(g,y)/var(y)
    #Predicted residual error sum of squares of the additive, dominance deviation and total genetic effects, respectively.
    PRESS.a=sum(c((as.matrix(ad)-a.p)^2))
    PRESS.d=sum(c((as.matrix(dom)-d.p)^2))
    PRESS.g=sum(c((as.matrix(g)-g.p)^2))
    
    ##Organizing the results of the current iterative process
    results=matrix(c(va,vd,vg,
                     cora,cord,corg,
                     ha,hd,hg,
                     PRESS.a,PRESS.d,PRESS.g),4,3,byrow=T)
    
    colnames(results) <-c("Addit.","Dom.", "Gen.")
    rownames(results) <-c("Variance","Correlation","heritability","PRESS")
    
    
    print("-----------------------Univariate GBLUP-AD-----------------------------------")  
    print("-----------------------------------------------------------------------------")
    print("-----------------------------------------------------------------------------")
    
    
    ##Setting the function to count the iteration time
    time=proc.time()
    time.new=as.numeric(time[1])
    time_iter=(time.new-time.old)/60
    time.old=time.new
    
    ##Setting the iterative features to be printed in each iteration
    t=t(as.matrix(c(dif,iter,time_iter)))
    colnames(t) <- c("Dif:" , "Iter:", "Time/iter:")
    print(t)
    
    print("-----------------------------------------------------------------------------")
    
    print("Statistics:")
    print(results)
    
    print("-----------------------------------------------------------------------------")
    
    if ((iter==maxiter)) break
  } 
  
  ##Saving the results of the analysis
  #Saving the results for of the analysis of the data constructed in the 0.3 heritability scenario
  if (i==1) {
    save(results,file = "Results_GBLUP-AD_UV_PH_h=0.3.RData")
  }
  if (i==2) {
    save(results,file = "Results_GBLUP-AD_UV_EH_h=0.3.RData")
  }
  if (i==3) {
    save(results,file = "Results_GBLUP-AD_UV_EL_h=0.3.RData")
  }
  if (i==4) {
    save(results,file = "Results_GBLUP-AD_UV_ERN_h=0.3.RData")
  }
  if (i==5) {
    save(results,file = "Results_GBLUP-AD_UV_KW_h=0.3.RData")
  }
  #Saving the results for of the analysis of the data constructed in the 0.5 heritability scenario
  if (i==6) {
    save(results,file = "Results_GBLUP-AD_UV_PH_h=0.5.RData")
  }
  if (i==7) {
    save(results,file = "Results_GBLUP-AD_UV_EH_h=0.5.RData")
  }
  if (i==8) {
    save(results,file = "Results_GBLUP-AD_UV_EL_h=0.5.RData")
  }
  if (i==9) {
    save(results,file = "Results_GBLUP-AD_UV_ERN_h=0.5.RData")
  }
  if (i==10) {
    save(results,file = "Results_GBLUP-AD_UV_KW_h=0.5.RData")
  }
  #Saving the results for of the analysis of the data constructed in the 0.7 heritability scenario
  if (i==11) {
    save(results,file = "Results_GBLUP-AD_UV_PH_h=0.7.RData")
  }
  if (i==12) {
    save(results,file = "Results_GBLUP-AD_UV_EH_h=0.7.RData")
  }
  if (i==13) {
    save(results,file = "Results_GBLUP-AD_UV_EL_h=0.7.RData")
  }
  if (i==14) {
    save(results,file = "Results_GBLUP-AD_UV_ERN_h=0.7.RData")
  }
  if (i==15) {
    save(results,file = "Results_GBLUP-AD_UV_KW_h=0.7.RData")
  }
  #Saving the results for of the analysis of the data constructed in the historical heritability scenario
  if (i==16) {
    save(results,file = "Results_GBLUP-AD_UV_PH_h=hist.RData")
  }
  if (i==17) {
    save(results,file = "Results_GBLUP-AD_UV_EH_h=hist.RData")
  }
  if (i==18) {
    save(results,file = "Results_GBLUP-AD_UV_EL_h=hist.RData")
  }
  if (i==19) {
    save(results,file = "Results_GBLUP-AD_UV_ERN_h=hist.RData")
  }
  if (i==20) {
    save(results,file = "Results_GBLUP-AD_UV_KW_h=hist.RData")
    rm(C,C1,C2,Ga,W1,Zds,beta,e,ha,hg,results,so,sol,t,va,va1,ve,ve1,vg,x,xx,xy,xz,y,zx,zy,zz,PRESS.a,PRESS.g,a.p,ad,comp,cora,corg,d.p,
       dif,g,g.p,iter,maxiter,n,n1,n2,time,time.new,time.old,time_iter,C3,Gd,hd,vd,vd1,PRESS.d,cord,dom)
  }
}

##Running the multivariate GBLUP-A model
for(i in 1:5) {
  time.new=0; time.old=0 #Setting parameters to compute the time of the iterations (time.new; time.old)
  
  if (i==1) {
    y1=y.h.PH.0.3; y2=y.h.EH.0.3; y3=y.h.EL.0.3; y4=y.h.ERN.0.3; y5=y.h.KW.0.3  #Pushing phenotypic data of all traits constructed in the 0.3 heritability scenario 
    comp=comp.PH.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.3[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==2) {
    y1=y.h.PH.0.5; y2=y.h.EH.0.5; y3=y.h.EL.0.5; y4=y.h.ERN.0.5; y5=y.h.KW.0.5  #Pushing phenotypic data of all traits constructed in the 0.5 heritability scenario 
    comp=comp.PH.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.5[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==3) {
    y1=y.h.PH.0.7; y2=y.h.EH.0.7; y3=y.h.EL.0.7; y4=y.h.ERN.0.7; y5=y.h.KW.0.7  #Pushing phenotypic data of all traits constructed in the 0.7 heritability scenario 
    comp=comp.PH.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.7[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==4) {
    y1=y.h.PH.0.3; y2=y.h.EH.0.5; y3=y.h.EL.0.3;  y4=y.h.ERN.0.7; y5=y.h.KW.0.7; #Pushing phenotypic data of all traits constructed in different heritability scenario determined randomly
    comp=comp.PH.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.5[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800]  #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait  
  }
  if (i==5) {
    y1=y.h.PH.hist; y2=y.h.EH.hist; y3=y.h.EL.hist; y4=y.h.ERN.hist; y5=y.h.KW.hist  #Pushing phenotypic data of all traits constructed in the historical heritability scenario 
    comp=comp.PH.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.hist[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  
  g.p1=a.p1+d.p1; g.p2=a.p2+d.p2; g.p3=a.p3+d.p3; g.p4=a.p4+d.p4; g.p5=a.p5+d.p5  # Computing the parametric total genetic value for all traits
  
  itermax=1500; iter=0 #Setting the number of iterations (itermax) and coming up with a variable to count the number of iterations (iter)
  
  ##Setting guesses of variance-covariance components of the residual effects to the first round of iterations
  ve11=100; ve12=0; ve13=0; ve14=0; ve15=0   
  ve21=0; ve22=100; ve23=0; ve24=0; ve25=0
  ve31=0; ve32=0; ve33=100; ve34=0; ve35=0  
  ve41=0; ve42=0; ve43=0; ve44=100; ve45=0  
  ve51=0; ve52=0; ve53=0; ve54=0; ve55=100
  
  ##Setting guesses of variance-covariance components of the additive effects to the first round of iterations
  va11=100; va12=0; va13=0; va14=0; va15=0  
  va21=0; va22=100; va23=0; va24=0; va25=0   
  va31=0; va32=0; va33=100; va34=0; va35=0 
  va41=0; va42=0; va43=0; va44=100; va45=0
  va51=0; va52=0; va53=0; va54=0; va55=100
  
  #Building the incidence matrix of the fixed effects considering all five traits
  x11=matrix(1,length(y1),1); x12=matrix(0,length(y1),1); x13=matrix(0,length(y1),1); x14=matrix(0,length(y1),1); x15=matrix(0,length(y1),1)
  x21=matrix(0,length(y2),1); x22=matrix(1,length(y2),1); x23=matrix(0,length(y2),1); x24=matrix(0,length(y2),1); x25=matrix(0,length(y2),1)
  x31=matrix(0,length(y3),1); x32=matrix(0,length(y3),1); x33=matrix(1,length(y3),1); x34=matrix(0,length(y3),1); x35=matrix(0,length(y3),1)
  x41=matrix(0,length(y4),1); x42=matrix(0,length(y4),1); x43=matrix(0,length(y4),1); x44=matrix(1,length(y4),1); x45=matrix(0,length(y4),1)
  x51=matrix(0,length(y4),1); x52=matrix(0,length(y4),1); x53=matrix(0,length(y4),1); x54=matrix(0,length(y4),1); x55=matrix(1,length(y4),1)
  X=rbind(cbind(x11,x12,x13,x14,x15),
          cbind(x21,x22,x23,x24,x25),
          cbind(x31,x32,x33,x34,x35),
          cbind(x41,x42,x43,x44,x45),
          cbind(x51,x52,x53,x54,x55))
  
  
  N1=length(y1); N2=length(y2); N3=length(y3); N4=length(y4); N5=length(y5); N=N1+N2+N3+N4+N5 #Counting the number of observations for each trait and total number of observations for all traits
  q=ncol(X) #Counting the number of fixed effects
  
  #Building the incidence matrix of the random effects considering all five traits
  z11=diag(length(y1)); z12=matrix(0,nrow(A),ncol(A)); z13=matrix(0,nrow(A),ncol(A)); z14=matrix(0,nrow(A),ncol(A)); z15=matrix(0,nrow(A),ncol(A))
  z21=matrix(0,nrow(A),ncol(A)); z22=diag(length(y2)); z23=matrix(0,nrow(A),ncol(A)); z24=matrix(0,nrow(A),ncol(A)); z25=matrix(0,nrow(A),ncol(A))
  z31=matrix(0,nrow(A),ncol(A)); z32=matrix(0,nrow(A),ncol(A)); z33=diag(length(y3)); z34=matrix(0,nrow(A),ncol(A)); z35=matrix(0,nrow(A),ncol(A))
  z41=matrix(0,nrow(A),ncol(A)); z42=matrix(0,nrow(A),ncol(A)); z43=matrix(0,nrow(A),ncol(A)); z44=diag(length(y4)); z45=matrix(0,nrow(A),ncol(A))
  z51=matrix(0,nrow(A),ncol(A)); z52=matrix(0,nrow(A),ncol(A)); z53=matrix(0,nrow(A),ncol(A)); z54=matrix(0,nrow(A),ncol(A)); z55=diag(length(y5))
  Z=rbind(cbind(z11,z12,z13,z14,z15),
          cbind(z21,z22,z23,z24,z25),
          cbind(z31,z32,z33,z34,z35),
          cbind(z41,z42,z43,z44,z45),
          cbind(z51,z52,z53,z54,z55))
  
  ##Removing additional variables used to build the fixed and random incidence matrices
  rm(x11,x12,x13,x14,x15,
     x21,x22,x23,x24,x25, 
     x31,x32,x33,x34,x35,
     x41,x42,x43,x44,x45,
     x51,x52,x53,x54,x55,
     z11,z12,z13,z14,z15,
     z21,z22,z23,z24,z25,
     z31,z32,z33,z34,z35,
     z41,z42,z43,z44,z45,
     z51,z52,z53,z54,z55)
  
  W=cbind(X,Z) #Horizontal concatenation of the [X | Z]
  
  y=rbind(as.matrix(y1),as.matrix(y2),as.matrix(y3),as.matrix(y4),as.matrix(y5)) #Concatenation of the vectors containing the phenotypic data of all traits
  
  setwd("/home/jhonathan/Documentos/MVGBLUP/Results")
  repeat {
    
    ve11<-as.numeric(ve11); ve12<-as.numeric(ve12); ve13<-as.numeric(ve13); ve14<-as.numeric(ve14); ve15<-as.numeric(ve15)
    ve21<-as.numeric(ve21); ve22<-as.numeric(ve22); ve23<-as.numeric(ve23); ve24<-as.numeric(ve24); ve25<-as.numeric(ve25)
    ve31<-as.numeric(ve31); ve32<-as.numeric(ve32); ve33<-as.numeric(ve33); ve34<-as.numeric(ve34); ve35<-as.numeric(ve35)
    ve41<-as.numeric(ve41); ve42<-as.numeric(ve42); ve43<-as.numeric(ve43); ve44<-as.numeric(ve44); ve45<-as.numeric(ve45)
    ve51<-as.numeric(ve51); ve52<-as.numeric(ve52); ve53<-as.numeric(ve53); ve54<-as.numeric(ve54); ve55<-as.numeric(ve55)
    
    va11<-as.numeric(va11); va12<-as.numeric(va12); va13<-as.numeric(va13); va14<-as.numeric(va14); va15<-as.numeric(va15)
    va21<-as.numeric(va21); va22<-as.numeric(va22); va23<-as.numeric(va23); va24<-as.numeric(va24); va25<-as.numeric(va25)
    va31<-as.numeric(va31); va32<-as.numeric(va32); va33<-as.numeric(va33); va34<-as.numeric(va34); va35<-as.numeric(va35)
    va41<-as.numeric(va41); va42<-as.numeric(va42); va43<-as.numeric(va43); va44<-as.numeric(va44); va45<-as.numeric(va45)
    va51<-as.numeric(va51); va52<-as.numeric(va52); va53<-as.numeric(va53); va54<-as.numeric(va54); va55<-as.numeric(va55)
    
    ##Buiding the variance-covariance matrix of the residual effects among all traits
    R=(matrix(c(ve11,ve12,ve13,ve14,ve15,
                ve21,ve22,ve23,ve24,ve25,
                ve31,ve32,ve33,ve34,ve35,
                ve41,ve42,ve43,ve44,ve45,
                ve51,ve52,ve53,ve54,ve55),q,q,byrow=T))
    sig.R.inv=chol2inv(chol(kronecker(R,diag(N1))))
    
    ##Building the variance-covariance additive genetic matrix among all traits
    Ga=(matrix(c(va11,va12,va13,va14,va15,
                 va21,va22,va23,va24,va25,
                 va31,va32,va33,va34,va35,
                 va41,va42,va43,va44,va45,
                 va51,va52,va53,va54,va55),q,q,byrow=T))
    G=kronecker(Ga,A); 
    
    ##Building the inverse matrix of the left-hand side of the mixed model equations
    C11=t(X)%*%sig.R.inv%*%X; C12=t(X)%*%sig.R.inv%*%Z;
    C21=t(Z)%*%sig.R.inv%*%X; C22=(t(Z)%*%sig.R.inv%*%Z+chol2inv(chol(G)));
    C=chol2inv(chol(rbind(cbind(C11,C12), cbind(C21,C22))))
    
    so=rbind(t(X)%*%sig.R.inv%*%y, t(Z)%*%sig.R.inv%*%y) #Building the right-hand side of the mixed model equations
    sol=C%*%so #Obtaining the solutions of the mixed model equations (EM - maximization step)
    
    beta=sol[1:q] #Pushing the fixed effects estimated from the solution of the mixed model equations
    ad=sol[(q+1):(q+N1+4*N2)] #Pushing the additive effects estimated from the solution of the mixed model equations
    ad1=ad[1:N1]; ad2=ad[(N1+1):(N1+N2)]; ad3=ad[(N1+N2+1):(N1+2*N2)]; ad4=ad[(N1+2*N2+1):(N1+3*N2)]; ad5=ad[(N1+3*N2+1):(N1+4*N2)] #Pushing the additive effects specific for each trait
    
    e=y-X%*%beta-Z%*%ad   ##Computing the vector of residuals
    e1=e[1:N1]; e2=e[(N1+1):(N1+N2)]; e3=e[(N1+N2+1):(N1+2*N2)]; e4=e[(N1+2*N2+1):(N1+3*N2)]; e5=e[(N1+3*N2+1):(N1+4*N2)] #Pushing the residual effects specific for each trait
    
    ##Obtaining the new variances components from the EM equations (EM - expection step)
    
    S=W%*%C%*%t(W)
    
    ve11.loop=(crossprod(e1,e1)+sum(diag(S[1:N1,1:N1])))/N1
    ve12.loop=(crossprod(e1,e2)+sum(diag(S[1:N1,(N1+1):(N1+N2)])))/N1
    ve13.loop=(crossprod(e1,e3)+sum(diag(S[1:N1,(N1+N2+1):(N1+2*N2)])))/N1
    ve14.loop=(crossprod(e1,e4)+sum(diag(S[1:N1,(N1+2*N2+1):(N1+3*N2)])))/N1
    ve15.loop=(crossprod(e1,e5)+sum(diag(S[1:N1,(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve21.loop=(crossprod(e2,e1)+sum(diag(S[(N1+1):(N1+N2),1:N1])))/N1
    ve22.loop=(crossprod(e2,e2)+sum(diag(S[(N1+1):(N1+N2),(N1+1):(N1+N2)])))/N1
    ve23.loop=(crossprod(e2,e3)+sum(diag(S[(N1+1):(N1+N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve24.loop=(crossprod(e2,e4)+sum(diag(S[(N1+1):(N1+N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve25.loop=(crossprod(e2,e5)+sum(diag(S[(N1+1):(N1+N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve31.loop=(crossprod(e3,e1)+sum(diag(S[(N1+N2+1):(N1+2*N2),1:N1])))/N1
    ve32.loop=(crossprod(e3,e2)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+1):(N1+N2)])))/N1
    ve33.loop=(crossprod(e3,e3)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve34.loop=(crossprod(e3,e4)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve35.loop=(crossprod(e3,e5)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve41.loop=(crossprod(e4,e1)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),1:N1])))/N1
    ve42.loop=(crossprod(e4,e2)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+1):(N1+N2)])))/N1
    ve43.loop=(crossprod(e4,e3)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve44.loop=(crossprod(e4,e4)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve45.loop=(crossprod(e4,e5)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve51.loop=(crossprod(e5,e1)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),1:N1])))/N1
    ve52.loop=(crossprod(e5,e2)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+1):(N1+N2)])))/N1
    ve53.loop=(crossprod(e5,e3)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve54.loop=(crossprod(e5,e4)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve55.loop=(crossprod(e5,e5)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    va11.loop=(crossprod(ad1,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1):(q+N1)])))/N1  
    va12.loop=(crossprod(ad1,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1):(q+N1+N2)])))/N1 
    va13.loop=(crossprod(ad1,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1+N2):(q+N1+2*N2)])))/N1 
    va14.loop=(crossprod(ad1,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1+2*N2):(q+N1+3*N2)])))/N1 
    va15.loop=(crossprod(ad1,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1+3*N2):(q+N1+4*N2)])))/N1 
    
    va21.loop=(crossprod(ad2,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+1):(q+N1)])))/N1 
    va22.loop=(crossprod(ad2,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+1):(q+N1+N2)])))/N1  
    va23.loop=(crossprod(ad2,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va24.loop=(crossprod(ad2,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va25.loop=(crossprod(ad2,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    va31.loop=(crossprod(ad3,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+1):(q+N1)])))/N1 
    va32.loop=(crossprod(ad3,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+1):(q+N1+N2)])))/N1  
    va33.loop=(crossprod(ad3,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va34.loop=(crossprod(ad3,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va35.loop=(crossprod(ad3,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    va41.loop=(crossprod(ad4,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+1):(q+N1)])))/N1 
    va42.loop=(crossprod(ad4,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+1):(q+N1+N2)])))/N1  
    va43.loop=(crossprod(ad4,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va44.loop=(crossprod(ad4,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va45.loop=(crossprod(ad4,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    va51.loop=(crossprod(ad5,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+1):(q+N1)])))/N1 
    va52.loop=(crossprod(ad5,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+1):(q+N1+N2)])))/N1  
    va53.loop=(crossprod(ad5,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va54.loop=(crossprod(ad5,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va55.loop=(crossprod(ad5,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    iter=iter+1 #Updating the current state of the iteration
    
    ##Computing the difference between the new and the old parameters obtained in the iterative process
    dif=max(abs(c((ve11.loop-ve11),(ve12.loop-ve12),(ve13.loop-ve13),(ve14.loop-ve14),(ve15.loop-ve15),
                  (ve21.loop-ve21),(ve22.loop-ve22),(ve23.loop-ve23),(ve24.loop-ve24),(ve25.loop-ve25),
                  (ve31.loop-ve31),(ve32.loop-ve32),(ve33.loop-ve33),(ve34.loop-ve34),(ve35.loop-ve35),        
                  (ve41.loop-ve41),(ve42.loop-ve42),(ve43.loop-ve43),(ve44.loop-ve44),(ve45.loop-ve45),        
                  (ve51.loop-ve51),(ve52.loop-ve52),(ve53.loop-ve53),(ve54.loop-ve54),(ve55.loop-ve55),            
                  (va11.loop-va11),(va12.loop-va12),(va13.loop-va13),(va14.loop-va14),(va15.loop-va15),
                  (va21.loop-va21),(va22.loop-va22),(va23.loop-va23),(va24.loop-va24),(va25.loop-va25),
                  (va31.loop-va31),(va32.loop-va32),(va33.loop-va33),(va34.loop-va34),(va35.loop-va35),
                  (va41.loop-va41),(va42.loop-va42),(va43.loop-va43),(va44.loop-va44),(va45.loop-va45),
                  (va51.loop-va51),(va52.loop-va52),(va53.loop-va53),(va54.loop-va54),(va55.loop-va55))))
    
    ##Updating variance components
    ve11=ve11.loop; ve12=ve12.loop; ve13=ve13.loop; ve14=ve14.loop; ve15=ve15.loop
    ve21=ve21.loop; ve22=ve22.loop; ve23=ve23.loop; ve24=ve24.loop; ve25=ve25.loop
    ve31=ve31.loop; ve32=ve32.loop; ve33=ve33.loop; ve34=ve34.loop; ve35=ve35.loop
    ve41=ve41.loop; ve42=ve42.loop; ve43=ve43.loop; ve44=ve44.loop; ve45=ve45.loop    
    ve51=ve51.loop; ve52=ve52.loop; ve53=ve53.loop; ve54=ve54.loop; ve55=ve55.loop
    va11=va11.loop; va12=va12.loop; va13=va13.loop; va14=va14.loop; va15=va15.loop
    va21=va21.loop; va22=va22.loop; va23=va23.loop; va24=va24.loop; va25=va25.loop
    va31=va31.loop; va32=va32.loop; va33=va33.loop; va34=va34.loop; va35=va35.loop
    va41=va41.loop; va42=va42.loop; va43=va43.loop; va44=va44.loop; va45=va45.loop
    va51=va51.loop; va52=va52.loop; va53=va53.loop; va54=va54.loop; va55=va55.loop
    
    
    ##Computing the model evaluation statistics 
    
    ##narrow-sense additive heritability of all five traits
    ha1=cov(ad1,y1)/var(y1); ha2=cov(ad2,y2)/var(y2); ha3=cov(ad3,y3)/var(y3); ha4=cov(ad4,y4)/var(y4); ha5=cov(ad5,y5)/var(y5)
    
    g1=ad1; g2=ad2; g3=ad3; g4=ad4; g5=ad5
    
    #broad-sense heritability of all five traits
    hg1=cov(g1,y1)/var(y1); hg2=cov(g2,y2)/var(y2); hg3=cov(g3,y3)/var(y3); hg4=cov(g4,y4)/var(y4); hg5=cov(g5,y5)/var(y5)
    
    #Correlations between the estimated and parametric additive effects of all five traits
    cora1=cor(ad1,a.p1); cora2=cor(ad2,a.p2); cora3=cor(ad3,a.p3); cora4=cor(ad4,a.p4); cora5=cor(ad5,a.p5)
    
    #Correlations between the estimated and parametric total genetic effects of all five traits
    corg1=cor(g1,g.p1); corg2=cor(g2,g.p2); corg3=cor(g3,g.p3); corg4=cor(g4,g.p4); corg5=cor(g5,g.p5)
    
    ##Organizing the results of the current iterative process
    results=rbind(c(ha1,ha2,ha3,ha4,ha5),
                  c(hg1,hg2,hg3,hg4,hg5),
                  c(cora1,cora2,cora3,cora4,cora5),
                  c(corg1,corg2,corg3,corg4,corg5))
    
    colnames(results) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(results) <- c('ha','hg','cora','corg')
    
    #Building variance-covariance matrix of the residual effects considering all traits in the current interative process
    vare=rbind(c(ve11,ve12,ve13,ve14,ve15),
               c(ve21,ve22,ve23,ve24,ve25),
               c(ve31,ve32,ve33,ve34,ve35), 
               c(ve41,ve42,ve43,ve44,ve45),        
               c(ve51,ve52,ve53,ve54,ve55))        
    
    colnames(vare) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(vare) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    if (i==1) {
      print("--------------------MV-GBLUP-A-(h=0.3)--Jhonathan Pedroso-------------------")  
    }
    if (i==2) {
      print("--------------------MV-GBLUP-A-(h=0.5)--Jhonathan Pedroso-------------------")  
    }
    if (i==3) {
      print("--------------------MV-GBLUP-A-(h=0.7)--Jhonathan Pedroso-------------------")  
    }
    if (i==4) {
      print("--------------------MV-GBLUP-A-(Mixed)--Jhonathan Pedroso-------------------")  
    }
    if (i==5) {
      print("---------------------MV-GBLUP-A-(Hist)--Jhonathan Pedroso-------------------")  
    }
    
    print("-----------------------------------------------------------------------------")
    print("-----------------------------------------------------------------------------")
    
    ##Setting the function to count the iteration time
    time=proc.time()
    time.new=as.numeric(time[1])
    time_iter=(time.new-time.old)/60
    time.old=time.new
    
    ##Setting the iterative features to be printed in each iteration
    t=t(as.matrix(c(dif,iter,time_iter)))
    colnames(t) <- c("Dif:" , "Iter:", "Time/iter:")
    print(t)
    
    print("-----------------------------------------------------------------------------")
    
    print("Statistics:")
    print(results)
    
    print("-----------------------------------------------------------------------------")
    
    print("Residual Variance-Covariace Matrix:",)
    print(vare) 
    
    #Building the variance-covariance additive genetic matrix among all traits
    vara=rbind(c(va11,va12,va13,va14,va15),
               c(va21,va22,va23,va24,va25),
               c(va31,va32,va33,va34,va35),
               c(va41,va42,va43,va44,va45),
               c(va51,va52,va53,va54,va55))
    
    colnames(vara) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(vara) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("-----------------------------------------------------------------------------")
    
    print("Additive Variance-Covariace Matrix:",)
    print(vara) 
    
    print("-----------------------------------------------------------------------------")
    
    varg=vara #Given that it is a additive model the additive effects is equivalent to the total genetic effects
    colnames(varg) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(varg) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("Genetic Variance-Covariace Matrix:",)
    print(varg) 
    
    print("-----------------------------------------------------------------------------")
    
    
    core.matrix=solve(diag(diag((vare))^0.5))%*%(vare)%*%solve(diag(diag((vare))^0.5)) #Obtaining the residual correlation matrix among all traits
    colnames(core.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(core.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("Residual Correlation Matrix:",)
    print(core.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    cora.matrix=solve(diag(diag((vara))^0.5))%*%(vara)%*%solve(diag(diag((vara))^0.5)) #Obtaining the additive effects correlation matrix among all traits
    colnames(cora.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(cora.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("Additive Correlation Matrix:",)
    print(cora.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    corg.matrix=solve(diag(diag((varg))^0.5))%*%(varg)%*%solve(diag(diag((varg))^0.5)) #Obtaining the total genetic effects correlation matrix among all traits
    colnames(corg.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(corg.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("Genetic Correlation Matrix:",)
    print(corg.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    ad1=as.matrix(ad1); ad2=as.matrix(ad2); ad3=as.matrix(ad3); ad4=as.matrix(ad4); ad5=as.matrix(ad5)
    g1=as.matrix(g1); g2=as.matrix(g2); g3=as.matrix(g3); g4=as.matrix(g4); g5=as.matrix(g5)
    
    #Computing predicted residual error sum of squares (PRESS) of the additive and total genetic effects, respectively.
    PRESS.a1=sum(c((ad1-a.p1)^2)); PRESS.a2=sum(c((ad2-a.p2)^2)); PRESS.a3=sum(c((ad3-a.p3)^2)); PRESS.a4=sum(c((ad4-a.p4)^2)); PRESS.a5=sum(c((ad5-a.p5)^2))
    PRESS.g1=sum(c((g1-g.p1)^2)); PRESS.g2=sum(c((g2-g.p2)^2)); PRESS.g3=sum(c((g3-g.p3)^2)); PRESS.g4=sum(c((g4-g.p4)^2)); PRESS.g5=sum(c((g5-g.p5)^2))
    
    #Array of the PRESS values
    PRESS.matrix=rbind(c(PRESS.a1,PRESS.a2,PRESS.a3,PRESS.a4,PRESS.a5),
                       c(PRESS.g1,PRESS.g2,PRESS.g3,PRESS.g4,PRESS.g5))
    
    colnames(PRESS.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(PRESS.matrix) <-c("PRESS.a","PRESS.g")
    
    print("PRESS values:",)
    print(PRESS.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    results=list(results=list(results),vare=list(vare),vara=list(vara),varg=list(varg),core.matrix=list(core.matrix),cora.matrix=list(cora.matrix),corg.matrix=list(corg.matrix),PRESS.matrix=list(PRESS.matrix),dif=list(dif),iter=list(iter))  
    
    if (i==1) {
      save(results,file = "Results_MV-GBLUP-A_h(0.3)_03-11-2014.RData")
    }
    if (i==2) {
      save(results,file = "Results_MV-GBLUP-A_h(0.5)_03-11-2014.RData")
    }
    if (i==3) {
      save(results,file = "Results_MV-GBLUP-A_h(0.7)_03-11-2014.RData")
    }
    if (i==4) {
      save(results,file = "Results_MV-GBLUP-A_h(Mixed)_03-11-2014.RData")
    }
    if (i==5) {
      save(results,file = "Results_MV-GBLUP-A_h(Hist)_03-11-2014.RData")
    }
    
    if (iter==itermax) break
  }
  ######################################  
}

##Running the multivariate GBLUP-AD model
for (i in 1:5) {
  time.new=0; time.old=0 #Setting parameters to compute the time of the iterations (time.new; time.old)
  
  if (i==1) {
    y1=y.h.PH.0.3; y2=y.h.EH.0.3; y3=y.h.EL.0.3; y4=y.h.ERN.0.3; y5=y.h.KW.0.3  #Pushing phenotypic data of all traits constructed in the 0.3 heritability scenario 
    comp=comp.PH.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.3[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==2) {
    y1=y.h.PH.0.5; y2=y.h.EH.0.5; y3=y.h.EL.0.5; y4=y.h.ERN.0.5; y5=y.h.KW.0.5  #Pushing phenotypic data of all traits constructed in the 0.5 heritability scenario 
    comp=comp.PH.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.5[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.5[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==3) {
    y1=y.h.PH.0.7; y2=y.h.EH.0.7; y3=y.h.EL.0.7; y4=y.h.ERN.0.7; y5=y.h.KW.0.7  #Pushing phenotypic data of all traits constructed in the 0.7 heritability scenario 
    comp=comp.PH.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.7[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  if (i==4) {
    y1=y.h.PH.0.3; y2=y.h.EH.0.5; y3=y.h.EL.0.3;  y4=y.h.ERN.0.7; y5=y.h.KW.0.7; #Pushing phenotypic data of all traits constructed in different heritability scenario determined randomly
    comp=comp.PH.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.0.5[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800]  #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.0.7[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.0.3[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait  
  }
  if (i==5) {
    y1=y.h.PH.hist; y2=y.h.EH.hist; y3=y.h.EL.hist; y4=y.h.ERN.hist; y5=y.h.KW.hist  #Pushing phenotypic data of all traits constructed in the historical heritability scenario 
    comp=comp.PH.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the plant height trait
    a.p1=comp[1:400]; d.p1=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the plant height trait
    comp=comp.EH.hist[-c(1,2,3,4)] #Ṕushing the parametric genetic effects related to the ear height trait
    a.p2=comp[1:400]; d.p2=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear height trait
    comp=comp.EL.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear length trait 
    a.p3=comp[1:400]; d.p3=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear length trait
    comp=comp.ERN.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the ear row number trait
    a.p4=comp[1:400]; d.p4=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the ear row number trait
    comp=comp.KW.hist[-c(1,2,3,4)] #Pushing the parametric genetic effects related to the kernel weight trait
    a.p5=comp[1:400]; d.p5=comp[401:800] #Pushing the parametric additive and dominance deviation effects related to the kernel weight trait
  }
  
  g.p1=a.p1+d.p1; g.p2=a.p2+d.p2; g.p3=a.p3+d.p3; g.p4=a.p4+d.p4; g.p5=a.p5+d.p5  # Computing the parametric total genetic value for all traits
  
  itermax=1500; iter=0 #Setting the number of iterations (itermax) and coming up with a variable to count the number of iterations (iter)
  
  ##Setting guesses of variance-covariance components of the residual effects to the first round of iterations
  ve11=100; ve12=0; ve13=0; ve14=0; ve15=0   
  ve21=0; ve22=100; ve23=0; ve24=0; ve25=0
  ve31=0; ve32=0; ve33=100; ve34=0; ve35=0  
  ve41=0; ve42=0; ve43=0; ve44=100; ve45=0  
  ve51=0; ve52=0; ve53=0; ve54=0; ve55=100
  
  ##Setting guesses of variance-covariance components of the additive effects to the first round of iterations
  va11=100; va12=0; va13=0; va14=0; va15=0  
  va21=0; va22=100; va23=0; va24=0; va25=0   
  va31=0; va32=0; va33=100; va34=0; va35=0 
  va41=0; va42=0; va43=0; va44=100; va45=0
  va51=0; va52=0; va53=0; va54=0; va55=100
  
  ##Setting guesses of variance-covariance components of the dominance deviation effects to the first round of iterations
  vd11=100; vd12=0; vd13=0; vd14=0; vd15=0  
  vd21=0; vd22=100; vd23=0; vd24=0; vd25=0   
  vd31=0; vd32=0; vd33=100; vd34=0; vd35=0 
  vd41=0; vd42=0; vd43=0; vd44=100; vd45=0
  vd51=0; vd52=0; vd53=0; vd54=0; vd55=100
  
  #Building the incidence matrix of the fixed effects considering all five traits
  x11=matrix(1,length(y1),1); x12=matrix(0,length(y1),1); x13=matrix(0,length(y1),1); x14=matrix(0,length(y1),1); x15=matrix(0,length(y1),1)
  x21=matrix(0,length(y2),1); x22=matrix(1,length(y2),1); x23=matrix(0,length(y2),1); x24=matrix(0,length(y2),1); x25=matrix(0,length(y2),1)
  x31=matrix(0,length(y3),1); x32=matrix(0,length(y3),1); x33=matrix(1,length(y3),1); x34=matrix(0,length(y3),1); x35=matrix(0,length(y3),1)
  x41=matrix(0,length(y4),1); x42=matrix(0,length(y4),1); x43=matrix(0,length(y4),1); x44=matrix(1,length(y4),1); x45=matrix(0,length(y4),1)
  x51=matrix(0,length(y4),1); x52=matrix(0,length(y4),1); x53=matrix(0,length(y4),1); x54=matrix(0,length(y4),1); x55=matrix(1,length(y4),1)
  X=rbind(cbind(x11,x12,x13,x14,x15),
          cbind(x21,x22,x23,x24,x25),
          cbind(x31,x32,x33,x34,x35),
          cbind(x41,x42,x43,x44,x45),
          cbind(x51,x52,x53,x54,x55))
  
  N1=length(y1); N2=length(y2); N3=length(y3); N4=length(y4); N5=length(y5); N=N1+N2+N3+N4+N5 #Counting the number of observations for each trait and total number of observations for all traits
  q=ncol(X) #Counting the number of fixed effects
  
  #Building the incidence matrix of the random effects considering all five traits
  z11=diag(length(y1)); z12=matrix(0,nrow(A),ncol(A)); z13=matrix(0,nrow(A),ncol(A)); z14=matrix(0,nrow(A),ncol(A)); z15=matrix(0,nrow(A),ncol(A))
  z21=matrix(0,nrow(A),ncol(A)); z22=diag(length(y2)); z23=matrix(0,nrow(A),ncol(A)); z24=matrix(0,nrow(A),ncol(A)); z25=matrix(0,nrow(A),ncol(A))
  z31=matrix(0,nrow(A),ncol(A)); z32=matrix(0,nrow(A),ncol(A)); z33=diag(length(y3)); z34=matrix(0,nrow(A),ncol(A)); z35=matrix(0,nrow(A),ncol(A))
  z41=matrix(0,nrow(A),ncol(A)); z42=matrix(0,nrow(A),ncol(A)); z43=matrix(0,nrow(A),ncol(A)); z44=diag(length(y4)); z45=matrix(0,nrow(A),ncol(A))
  z51=matrix(0,nrow(A),ncol(A)); z52=matrix(0,nrow(A),ncol(A)); z53=matrix(0,nrow(A),ncol(A)); z54=matrix(0,nrow(A),ncol(A)); z55=diag(length(y5))
  Z=rbind(cbind(z11,z12,z13,z14,z15),
          cbind(z21,z22,z23,z24,z25),
          cbind(z31,z32,z33,z34,z35),
          cbind(z41,z42,z43,z44,z45),
          cbind(z51,z52,z53,z54,z55))
  
  W=cbind(X,Z,Z)
  
  ##Removing additional variables used to build the fixed and random incidence matrices
  rm(x11,x12,x13,x14,x15,
     x21,x22,x23,x24,x25, 
     x31,x32,x33,x34,x35,
     x41,x42,x43,x44,x45,
     x51,x52,x53,x54,x55,
     z11,z12,z13,z14,z15,
     z21,z22,z23,z24,z25,
     z31,z32,z33,z34,z35,
     z41,z42,z43,z44,z45,
     z51,z52,z53,z54,z55)
  
  y=rbind(as.matrix(y1),as.matrix(y2),as.matrix(y3),as.matrix(y4),as.matrix(y5))
  
  setwd("/home/jhonathan/Documentos/MVGBLUP/Results")
  repeat {
    ve11<-as.numeric(ve11); ve12<-as.numeric(ve12); ve13<-as.numeric(ve13); ve14<-as.numeric(ve14); ve15<-as.numeric(ve15)
    ve21<-as.numeric(ve21); ve22<-as.numeric(ve22); ve23<-as.numeric(ve23); ve24<-as.numeric(ve24); ve25<-as.numeric(ve25)
    ve31<-as.numeric(ve31); ve32<-as.numeric(ve32); ve33<-as.numeric(ve33); ve34<-as.numeric(ve34); ve35<-as.numeric(ve35)
    ve41<-as.numeric(ve41); ve42<-as.numeric(ve42); ve43<-as.numeric(ve43); ve44<-as.numeric(ve44); ve45<-as.numeric(ve45)
    ve51<-as.numeric(ve51); ve52<-as.numeric(ve52); ve53<-as.numeric(ve53); ve54<-as.numeric(ve54); ve55<-as.numeric(ve55)
    
    va11<-as.numeric(va11); va12<-as.numeric(va12); va13<-as.numeric(va13); va14<-as.numeric(va14); va15<-as.numeric(va15)
    va21<-as.numeric(va21); va22<-as.numeric(va22); va23<-as.numeric(va23); va24<-as.numeric(va24); va25<-as.numeric(va25)
    va31<-as.numeric(va31); va32<-as.numeric(va32); va33<-as.numeric(va33); va34<-as.numeric(va34); va35<-as.numeric(va35)
    va41<-as.numeric(va41); va42<-as.numeric(va42); va43<-as.numeric(va43); va44<-as.numeric(va44); va45<-as.numeric(va45)
    va51<-as.numeric(va51); va52<-as.numeric(va52); va53<-as.numeric(va53); va54<-as.numeric(va54); va55<-as.numeric(va55)
    
    vd11<-as.numeric(vd11); vd12<-as.numeric(vd12); vd13<-as.numeric(vd13); vd14<-as.numeric(vd14); vd15<-as.numeric(vd15)
    vd21<-as.numeric(vd21); vd22<-as.numeric(vd22); vd23<-as.numeric(vd23); vd24<-as.numeric(vd24); vd25<-as.numeric(vd25)
    vd31<-as.numeric(vd31); vd32<-as.numeric(vd32); vd33<-as.numeric(vd33); vd34<-as.numeric(vd34); vd35<-as.numeric(vd35)
    vd41<-as.numeric(vd41); vd42<-as.numeric(vd42); vd43<-as.numeric(vd43); vd44<-as.numeric(vd44); vd45<-as.numeric(vd45)
    vd51<-as.numeric(vd51); vd52<-as.numeric(vd52); vd53<-as.numeric(vd53); vd54<-as.numeric(vd54); vd55<-as.numeric(vd55)
    
    ##Buiding the variance-covariance matrix of the residual effects among all traits  
    R=(matrix(c(ve11,ve12,ve13,ve14,ve15,
                ve21,ve22,ve23,ve24,ve25,
                ve31,ve32,ve33,ve34,ve35,
                ve41,ve42,ve43,ve44,ve45,
                ve51,ve52,ve53,ve54,ve55),q,q,byrow=T))
    sig.R.inv=chol2inv(chol(kronecker(R,diag(N1))))
    
    ##Building the variance-covariance additive genetic matrix among all traits
    Ga=(matrix(c(va11,va12,va13,va14,va15,
                 va21,va22,va23,va24,va25,
                 va31,va32,va33,va34,va35,
                 va41,va42,va43,va44,va45,
                 va51,va52,va53,va54,va55),q,q,byrow=T))
    
    ##Building the variance-covariance dominance genetic matrix among all traits
    Gd=(matrix(c(vd11,vd12,vd13,vd14,vd15,
                 vd21,vd22,vd23,vd24,vd25,
                 vd31,vd32,vd33,vd34,vd35,
                 vd41,vd42,vd43,vd44,vd45,
                 vd51,vd52,vd53,vd54,vd55),q,q,byrow=T))
    G.a=kronecker(Ga,A); 
    G.d=kronecker(Gd,D); 
    
    ##Building the inverse matrix of the left-hand side of the mixed model equations
    C11=t(X)%*%sig.R.inv%*%X; C12=t(X)%*%sig.R.inv%*%Z;C13=t(X)%*%sig.R.inv%*%Z;
    C21=t(Z)%*%sig.R.inv%*%X; C22=(t(Z)%*%sig.R.inv%*%Z+chol2inv(chol(G.a)));C23=t(Z)%*%sig.R.inv%*%Z;
    C31=t(Z)%*%sig.R.inv%*%X; C32=t(Z)%*%sig.R.inv%*%Z;C33=(t(Z)%*%sig.R.inv%*%Z+chol2inv(chol(G.d)));
    
    
    C=chol2inv(chol(rbind(cbind(C11,C12,C13),
                          cbind(C21,C22,C23),
                          cbind(C31,C32,C33))))
    
    
    so=rbind(t(X)%*%sig.R.inv%*%y, t(Z)%*%sig.R.inv%*%y, t(Z)%*%sig.R.inv%*%y) #Building the right-hand side of the mixed model equations
    sol=C%*%so #Obtaining the solutions of the mixed model equations (EM - maximization step)
    
    beta=sol[1:q] #Pushing the fixed effects estimated from the solution of the mixed model equations
    ad=sol[(q+1):(q+N1+4*N2)] #Pushing the additive effects estimated from the solution of the mixed model equations
    ad1=ad[1:N1]; ad2=ad[(N1+1):(N1+N2)]; ad3=ad[(N1+N2+1):(N1+2*N2)]; ad4=ad[(N1+2*N2+1):(N1+3*N2)]; ad5=ad[(N1+3*N2+1):(N1+4*N2)] #Pushing the additive effects specific for each trait
    
    dom=sol[(q+N1+4*N2+1):(q+N1+9*N2)] #Pushing the dominance deviation effects estimated from the solution of the mixed model equations
    dom1=dom[1:N1]; dom2=dom[(N1+1):(N1+N2)]; dom3=dom[(N1+N2+1):(N1+2*N2)]; dom4=dom[(N1+2*N2+1):(N1+3*N2)]; dom5=dom[(N1+3*N2+1):(N1+4*N2)] #Pushing the dominance deviation effects specific for each trait
    
    
    e=y-X%*%beta-Z%*%ad-Z%*%dom #Computing the vector of residuals
    e1=e[1:N1]; e2=e[(N1+1):(N1+N2)]; e3=e[(N1+N2+1):(N1+2*N2)]; e4=e[(N1+2*N2+1):(N1+3*N2)]; e5=e[(N1+3*N2+1):(N1+4*N2)] #Pushing the residual effects specific for each trait
    
    ##Obtaining the new variances components from the EM equations (EM - expection step)
    
    S=W%*%C%*%t(W)
    
    ve11.loop=(crossprod(e1,e1)+sum(diag(S[1:N1,1:N1])))/N1
    ve12.loop=(crossprod(e1,e2)+sum(diag(S[1:N1,(N1+1):(N1+N2)])))/N1
    ve13.loop=(crossprod(e1,e3)+sum(diag(S[1:N1,(N1+N2+1):(N1+2*N2)])))/N1
    ve14.loop=(crossprod(e1,e4)+sum(diag(S[1:N1,(N1+2*N2+1):(N1+3*N2)])))/N1
    ve15.loop=(crossprod(e1,e5)+sum(diag(S[1:N1,(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve21.loop=(crossprod(e2,e1)+sum(diag(S[(N1+1):(N1+N2),1:N1])))/N1
    ve22.loop=(crossprod(e2,e2)+sum(diag(S[(N1+1):(N1+N2),(N1+1):(N1+N2)])))/N1
    ve23.loop=(crossprod(e2,e3)+sum(diag(S[(N1+1):(N1+N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve24.loop=(crossprod(e2,e4)+sum(diag(S[(N1+1):(N1+N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve25.loop=(crossprod(e2,e5)+sum(diag(S[(N1+1):(N1+N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve31.loop=(crossprod(e3,e1)+sum(diag(S[(N1+N2+1):(N1+2*N2),1:N1])))/N1
    ve32.loop=(crossprod(e3,e2)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+1):(N1+N2)])))/N1
    ve33.loop=(crossprod(e3,e3)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve34.loop=(crossprod(e3,e4)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve35.loop=(crossprod(e3,e5)+sum(diag(S[(N1+N2+1):(N1+2*N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve41.loop=(crossprod(e4,e1)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),1:N1])))/N1
    ve42.loop=(crossprod(e4,e2)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+1):(N1+N2)])))/N1
    ve43.loop=(crossprod(e4,e3)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve44.loop=(crossprod(e4,e4)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve45.loop=(crossprod(e4,e5)+sum(diag(S[(N1+2*N2+1):(N1+3*N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    ve51.loop=(crossprod(e5,e1)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),1:N1])))/N1
    ve52.loop=(crossprod(e5,e2)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+1):(N1+N2)])))/N1
    ve53.loop=(crossprod(e5,e3)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+N2+1):(N1+2*N2)])))/N1
    ve54.loop=(crossprod(e5,e4)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+2*N2+1):(N1+3*N2)])))/N1
    ve55.loop=(crossprod(e5,e5)+sum(diag(S[(N1+3*N2+1):(N1+4*N2),(N1+3*N2+1):(N1+4*N2)])))/N1
    
    
    va11.loop=(crossprod(ad1,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1):(q+N1)])))/N1  
    va12.loop=(crossprod(ad1,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1):(q+N1+N2)])))/N1 
    va13.loop=(crossprod(ad1,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1+N2):(q+N1+2*N2)])))/N1 
    va14.loop=(crossprod(ad1,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1+2*N2):(q+N1+3*N2)])))/N1 
    va15.loop=(crossprod(ad1,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+1):(q+N1),(q+1+N1+3*N2):(q+N1+4*N2)])))/N1 
    
    va21.loop=(crossprod(ad2,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+1):(q+N1)])))/N1 
    va22.loop=(crossprod(ad2,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+1):(q+N1+N2)])))/N1  
    va23.loop=(crossprod(ad2,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va24.loop=(crossprod(ad2,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va25.loop=(crossprod(ad2,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+1):(q+N1+N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    va31.loop=(crossprod(ad3,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+1):(q+N1)])))/N1 
    va32.loop=(crossprod(ad3,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+1):(q+N1+N2)])))/N1  
    va33.loop=(crossprod(ad3,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va34.loop=(crossprod(ad3,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va35.loop=(crossprod(ad3,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+N2+1):(q+N1+2*N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    va41.loop=(crossprod(ad4,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+1):(q+N1)])))/N1 
    va42.loop=(crossprod(ad4,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+1):(q+N1+N2)])))/N1  
    va43.loop=(crossprod(ad4,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va44.loop=(crossprod(ad4,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va45.loop=(crossprod(ad4,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+2*N2+1):(q+N1+3*N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    va51.loop=(crossprod(ad5,Ga1)%*%ad1+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+1):(q+N1)])))/N1 
    va52.loop=(crossprod(ad5,Ga1)%*%ad2+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+1):(q+N1+N2)])))/N1  
    va53.loop=(crossprod(ad5,Ga1)%*%ad3+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+N2+1):(q+N1+2*N2)])))/N1  
    va54.loop=(crossprod(ad5,Ga1)%*%ad4+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+2*N2+1):(q+N1+3*N2)])))/N1  
    va55.loop=(crossprod(ad5,Ga1)%*%ad5+sum(diag(Ga1%*%C[(q+N1+3*N2+1):(q+N1+4*N2),(q+N1+3*N2+1):(q+N1+4*N2)])))/N1  
    
    
    vd11.loop=(crossprod(dom1,Gd1)%*%dom1+sum(diag(Gd1%*%C[(q+N1+4*N2+1):(q+N1+5*N2),(q+N1+4*N2+1):(q+N1+5*N2)])))/N1  
    vd12.loop=(crossprod(dom1,Gd1)%*%dom2+sum(diag(Gd1%*%C[(q+N1+4*N2+1):(q+N1+5*N2),(q+N1+5*N2+1):(q+N1+6*N2)])))/N1 
    vd13.loop=(crossprod(dom1,Gd1)%*%dom3+sum(diag(Gd1%*%C[(q+N1+4*N2+1):(q+N1+5*N2),(q+N1+6*N2+1):(q+N1+7*N2)])))/N1 
    vd14.loop=(crossprod(dom1,Gd1)%*%dom4+sum(diag(Gd1%*%C[(q+N1+4*N2+1):(q+N1+5*N2),(q+N1+7*N2+1):(q+N1+8*N2)])))/N1 
    vd15.loop=(crossprod(dom1,Gd1)%*%dom5+sum(diag(Gd1%*%C[(q+N1+4*N2+1):(q+N1+5*N2),(q+N1+8*N2+1):(q+N1+9*N2)])))/N1 
    
    vd21.loop=(crossprod(dom2,Gd1)%*%dom1+sum(diag(Gd1%*%C[(q+N1+5*N2+1):(q+N1+6*N2),(q+N1+4*N2+1):(q+N1+5*N2)])))/N1 
    vd22.loop=(crossprod(dom2,Gd1)%*%dom2+sum(diag(Gd1%*%C[(q+N1+5*N2+1):(q+N1+6*N2),(q+N1+5*N2+1):(q+N1+6*N2)])))/N1  
    vd23.loop=(crossprod(dom2,Gd1)%*%dom3+sum(diag(Gd1%*%C[(q+N1+5*N2+1):(q+N1+6*N2),(q+N1+6*N2+1):(q+N1+7*N2)])))/N1  
    vd24.loop=(crossprod(dom2,Gd1)%*%dom4+sum(diag(Gd1%*%C[(q+N1+5*N2+1):(q+N1+6*N2),(q+N1+7*N2+1):(q+N1+8*N2)])))/N1  
    vd25.loop=(crossprod(dom2,Gd1)%*%dom5+sum(diag(Gd1%*%C[(q+N1+5*N2+1):(q+N1+6*N2),(q+N1+8*N2+1):(q+N1+9*N2)])))/N1  
    
    vd31.loop=(crossprod(dom3,Gd1)%*%dom1+sum(diag(Gd1%*%C[(q+N1+6*N2+1):(q+N1+7*N2),(q+N1+4*N2+1):(q+N1+5*N2)])))/N1 
    vd32.loop=(crossprod(dom3,Gd1)%*%dom2+sum(diag(Gd1%*%C[(q+N1+6*N2+1):(q+N1+7*N2),(q+N1+5*N2+1):(q+N1+6*N2)])))/N1  
    vd33.loop=(crossprod(dom3,Gd1)%*%dom3+sum(diag(Gd1%*%C[(q+N1+6*N2+1):(q+N1+7*N2),(q+N1+6*N2+1):(q+N1+7*N2)])))/N1  
    vd34.loop=(crossprod(dom3,Gd1)%*%dom4+sum(diag(Gd1%*%C[(q+N1+6*N2+1):(q+N1+7*N2),(q+N1+7*N2+1):(q+N1+8*N2)])))/N1  
    vd35.loop=(crossprod(dom3,Gd1)%*%dom5+sum(diag(Gd1%*%C[(q+N1+6*N2+1):(q+N1+7*N2),(q+N1+8*N2+1):(q+N1+9*N2)])))/N1  
    
    vd41.loop=(crossprod(dom4,Gd1)%*%dom1+sum(diag(Gd1%*%C[(q+N1+7*N2+1):(q+N1+8*N2),(q+N1+4*N2+1):(q+N1+5*N2)])))/N1 
    vd42.loop=(crossprod(dom4,Gd1)%*%dom2+sum(diag(Gd1%*%C[(q+N1+7*N2+1):(q+N1+8*N2),(q+N1+5*N2+1):(q+N1+6*N2)])))/N1  
    vd43.loop=(crossprod(dom4,Gd1)%*%dom3+sum(diag(Gd1%*%C[(q+N1+7*N2+1):(q+N1+8*N2),(q+N1+6*N2+1):(q+N1+7*N2)])))/N1  
    vd44.loop=(crossprod(dom4,Gd1)%*%dom4+sum(diag(Gd1%*%C[(q+N1+7*N2+1):(q+N1+8*N2),(q+N1+7*N2+1):(q+N1+8*N2)])))/N1  
    vd45.loop=(crossprod(dom4,Gd1)%*%dom5+sum(diag(Gd1%*%C[(q+N1+7*N2+1):(q+N1+8*N2),(q+N1+8*N2+1):(q+N1+9*N2)])))/N1  
    
    vd51.loop=(crossprod(dom5,Gd1)%*%dom1+sum(diag(Gd1%*%C[(q+N1+8*N2+1):(q+N1+9*N2),(q+N1+4*N2+1):(q+N1+5*N2)])))/N1 
    vd52.loop=(crossprod(dom5,Gd1)%*%dom2+sum(diag(Gd1%*%C[(q+N1+8*N2+1):(q+N1+9*N2),(q+N1+5*N2+1):(q+N1+6*N2)])))/N1  
    vd53.loop=(crossprod(dom5,Gd1)%*%dom3+sum(diag(Gd1%*%C[(q+N1+8*N2+1):(q+N1+9*N2),(q+N1+6*N2+1):(q+N1+7*N2)])))/N1  
    vd54.loop=(crossprod(dom5,Gd1)%*%dom4+sum(diag(Gd1%*%C[(q+N1+8*N2+1):(q+N1+9*N2),(q+N1+7*N2+1):(q+N1+8*N2)])))/N1  
    vd55.loop=(crossprod(dom5,Gd1)%*%dom5+sum(diag(Gd1%*%C[(q+N1+8*N2+1):(q+N1+9*N2),(q+N1+8*N2+1):(q+N1+9*N2)])))/N1  
    
    
    iter=iter+1 #Updating the current state of the iteration
    
    ##Computing the difference between the new and the old parameters obtained in the iterative process
    dif=max(abs(c((ve11.loop-ve11),(ve12.loop-ve12),(ve13.loop-ve13),(ve14.loop-ve14),(ve15.loop-ve15),
                  (ve21.loop-ve21),(ve22.loop-ve22),(ve23.loop-ve23),(ve24.loop-ve24),(ve25.loop-ve25),
                  (ve31.loop-ve31),(ve32.loop-ve32),(ve33.loop-ve33),(ve34.loop-ve34),(ve35.loop-ve35),        
                  (ve41.loop-ve41),(ve42.loop-ve42),(ve43.loop-ve43),(ve44.loop-ve44),(ve45.loop-ve45),        
                  (ve51.loop-ve51),(ve52.loop-ve52),(ve53.loop-ve53),(ve54.loop-ve54),(ve55.loop-ve55),            
                  (va11.loop-va11),(va12.loop-va12),(va13.loop-va13),(va14.loop-va14),(va15.loop-va15),
                  (va21.loop-va21),(va22.loop-va22),(va23.loop-va23),(va24.loop-va24),(va25.loop-va25),
                  (va31.loop-va31),(va32.loop-va32),(va33.loop-va33),(va34.loop-va34),(va35.loop-va35),
                  (va41.loop-va41),(va42.loop-va42),(va43.loop-va43),(va44.loop-va44),(va45.loop-va45),
                  (va51.loop-va51),(va52.loop-va52),(va53.loop-va53),(va54.loop-va54),(va55.loop-va55),
                  (vd11.loop-vd11),(vd12.loop-vd12),(vd13.loop-vd13),(vd14.loop-vd14),(vd15.loop-vd15),
                  (vd21.loop-vd21),(vd22.loop-vd22),(vd23.loop-vd23),(vd24.loop-vd24),(vd25.loop-vd25),
                  (vd31.loop-vd31),(vd32.loop-vd32),(vd33.loop-vd33),(vd34.loop-vd34),(vd35.loop-vd35),
                  (vd41.loop-vd41),(vd42.loop-vd42),(vd43.loop-vd43),(vd44.loop-vd44),(vd45.loop-vd45),
                  (vd51.loop-vd51),(vd52.loop-vd52),(vd53.loop-vd53),(vd54.loop-vd54),(vd55.loop-vd55))))
    
    
    ##Updating variance components
    ve11=ve11.loop; ve12=ve12.loop; ve13=ve13.loop; ve14=ve14.loop; ve15=ve15.loop
    ve21=ve21.loop; ve22=ve22.loop; ve23=ve23.loop; ve24=ve24.loop; ve25=ve25.loop
    ve31=ve31.loop; ve32=ve32.loop; ve33=ve33.loop; ve34=ve34.loop; ve35=ve35.loop
    ve41=ve41.loop; ve42=ve42.loop; ve43=ve43.loop; ve44=ve44.loop; ve45=ve45.loop    
    ve51=ve51.loop; ve52=ve52.loop; ve53=ve53.loop; ve54=ve54.loop; ve55=ve55.loop
    va11=va11.loop; va12=va12.loop; va13=va13.loop; va14=va14.loop; va15=va15.loop
    va21=va21.loop; va22=va22.loop; va23=va23.loop; va24=va24.loop; va25=va25.loop
    va31=va31.loop; va32=va32.loop; va33=va33.loop; va34=va34.loop; va35=va35.loop
    va41=va41.loop; va42=va42.loop; va43=va43.loop; va44=va44.loop; va45=va45.loop
    va51=va51.loop; va52=va52.loop; va53=va53.loop; va54=va54.loop; va55=va55.loop
    vd11=vd11.loop; vd12=vd12.loop; vd13=vd13.loop; vd14=vd14.loop; vd15=vd15.loop
    vd21=vd21.loop; vd22=vd22.loop; vd23=vd23.loop; vd24=vd24.loop; vd25=vd25.loop
    vd31=vd31.loop; vd32=vd32.loop; vd33=vd33.loop; vd34=vd34.loop; vd35=vd35.loop
    vd41=vd41.loop; vd42=vd42.loop; vd43=vd43.loop; vd44=vd44.loop; vd45=vd45.loop
    vd51=vd51.loop; vd52=vd52.loop; vd53=vd53.loop; vd54=vd54.loop; vd55=vd55.loop
    
    ##Computing the model evaluation statistics 
    
    #narrow-sense additive heritability of all five traits
    ha1=cov(ad1,y1)/var(y1); ha2=cov(ad2,y2)/var(y2); ha3=cov(ad3,y3)/var(y3); ha4=cov(ad4,y4)/var(y4); ha5=cov(ad5,y5)/var(y5)
    
    
    #narrow-sense dominance heritability of all five traits
    hd1=cov(dom1,y1)/var(y1); hd2=cov(dom2,y2)/var(y2); hd3=cov(dom3,y3)/var(y3); hd4=cov(dom4,y4)/var(y4); hd5=cov(dom5,y5)/var(y5)
    
    g1=ad1+dom1; g2=ad2+dom2; g3=ad3+dom3; g4=ad4+dom4; g5=ad5+dom5
    
    #broad-sense heritability of all five traits
    hg1=cov(g1,y1)/var(y1); hg2=cov(g2,y2)/var(y2); hg3=cov(g3,y3)/var(y3); hg4=cov(g4,y4)/var(y4); hg5=cov(g5,y5)/var(y5)
    
    #Correlations between the estimated and parametric additive effects of all five traits
    cora1=cor(ad1,a.p1); cora2=cor(ad2,a.p2); cora3=cor(ad3,a.p3); cora4=cor(ad4,a.p4); cora5=cor(ad5,a.p5)
    
    #Correlations between the estimated and parametric dominance deviation effects of all five traits
    cord1=cor(dom1,d.p1); cord2=cor(dom2,d.p2); cord3=cor(dom3,d.p3); cord4=cor(dom4,d.p4); cord5=cor(dom5,d.p5)
    
    #Correlations between the estimated and parametric total genetic effects of all five traits
    corg1=cor(g1,g.p1); corg2=cor(g2,g.p2); corg3=cor(g3,g.p3); corg4=cor(g4,g.p4); corg5=cor(g5,g.p5)
    
    ##Organizing the results of the current iterative process
    results=rbind(c(ha1,ha2,ha3,ha4,ha5),
                  c(hd1,hd2,hd3,hd4,hd5),
                  c(hg1,hg2,hg3,hg4,hg5),
                  c(cora1,cora2,cora3,cora4,cora5),
                  c(cord1,cord2,cord3,cord4,cord5),
                  c(corg1,corg2,corg3,corg4,corg5))
    
    
    colnames(results) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(results) <- c('ha','hd','hg','cora','cord','corg')
    
    #Building variance-covariance matrix of the residual effects considering all traits in the current interative process
    vare=rbind(c(ve11,ve12,ve13,ve14,ve15),
               c(ve21,ve22,ve23,ve24,ve25),
               c(ve31,ve32,ve33,ve34,ve35), 
               c(ve41,ve42,ve43,ve44,ve45),        
               c(ve51,ve52,ve53,ve54,ve55))        
    
    colnames(vare) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(vare) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    if (i==1) {
      print("--------------------MV-GBLUP-AD-(h=0.3)--Jhonathan Pedroso-------------------")  
    }
    if (i==2) {
      print("--------------------MV-GBLUP-AD-(h=0.5)--Jhonathan Pedroso-------------------")  
    }
    if (i==3) {
      print("--------------------MV-GBLUP-AD-(h=0.7)--Jhonathan Pedroso-------------------")  
    }
    if (i==4) {
      print("--------------------MV-GBLUP-AD-(Mixed)--Jhonathan Pedroso-------------------")  
    }
    if (i==5) {
      print("---------------------MV-GBLUP-AD-(Hist)--Jhonathan Pedroso-------------------")  
    }
    
    print("-----------------------------------------------------------------------------")
    print("-----------------------------------------------------------------------------")
    
    ##Setting the function to count the iteration time
    time=proc.time()
    time.new=as.numeric(time[1])
    time_iter=(time.new-time.old)/60
    time.old=time.new
    
    ##Setting the iterative features to be printed in each iteration
    t=t(as.matrix(c(dif,iter,time_iter)))
    colnames(t) <- c("Dif:" , "Iter:", "Time/iter:")
    print(t)
    
    print("-----------------------------------------------------------------------------")
    
    print("Statistics:")
    print(results)
    
    print("-----------------------------------------------------------------------------")
    
    print("Residual Variance-Covariace Matrix:",)
    print(vare) 
    
    #Building the variance-covariance additive genetic matrix among all traits
    vara=rbind(c(va11,va12,va13,va14,va15),
               c(va21,va22,va23,va24,va25),
               c(va31,va32,va33,va34,va35),
               c(va41,va42,va43,va44,va45),
               c(va51,va52,va53,va54,va55))
    
    colnames(vara) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(vara) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("-----------------------------------------------------------------------------")
    
    print("Additive Variance-Covariace Matrix:",)
    print(vara) 
    
    print("-----------------------------------------------------------------------------")
    
    #Building the variance-covariance dominance genetic matrix among all traits
    vard=rbind(c(vd11,vd12,vd13,vd14,vd15),
               c(vd21,vd22,vd23,vd24,vd25),
               c(vd31,vd32,vd33,vd34,vd35),
               c(vd41,vd42,vd43,vd44,vd45),
               c(vd51,vd52,vd53,vd54,vd55))
    
    colnames(vard) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(vard) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("Dominance Variance-Covariace Matrix:",)
    print(vard) 
    
    print("-----------------------------------------------------------------------------")
    
    varg=vara+vard
    colnames(varg) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(varg) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    
    print("Genetic Variance-Covariace Matrix:",)
    print(varg) 
    
    print("-----------------------------------------------------------------------------")
    
    core.matrix=solve(diag(diag((vare))^0.5))%*%(vare)%*%solve(diag(diag((vare))^0.5)) #Obtaining the residual correlation matrix among all traits
    colnames(core.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(core.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    print("Residual Correlation Matrix:",)
    print(core.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    cora.matrix=solve(diag(diag((vara))^0.5))%*%(vara)%*%solve(diag(diag((vara))^0.5)) #Obtaining the additive effects correlation matrix among all traits
    colnames(cora.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(cora.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    print("Additive Correlation Matrix:",)
    print(cora.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    cord.matrix=solve(diag(diag((vard))^0.5))%*%(vard)%*%solve(diag(diag((vard))^0.5)) #Obtaining the dominance deviation effects correlation matrix among all traits
    colnames(cord.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(cord.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    print("Dominance Correlation Matrix:",)
    print(cord.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    corg.matrix=solve(diag(diag((varg))^0.5))%*%(varg)%*%solve(diag(diag((varg))^0.5)) #Obtaining the total genetic effects correlation matrix among all traits
    colnames(corg.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(corg.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    print("Genetic Correlation Matrix:",)
    print(corg.matrix) 
    
    print("-----------------------------------------------------------------------------")
    
    ad1=as.matrix(ad1); ad2=as.matrix(ad2); ad3=as.matrix(ad3); ad4=as.matrix(ad4); ad5=as.matrix(ad5)
    dom1=as.matrix(dom1); dom2=as.matrix(dom2); dom3=as.matrix(dom3); dom4=as.matrix(dom4); dom5=as.matrix(dom5)
    g1=as.matrix(g1); g2=as.matrix(g2); g3=as.matrix(g3); g4=as.matrix(g4); g5=as.matrix(g5)
    
    #Computing predicted residual error sum of squares (PRESS) of the additive, dominance desviation and total genetic effects, respectively. 
    PRESS.a1=sum(c((ad1-a.p1)^2)); PRESS.a2=sum(c((ad2-a.p2)^2)); PRESS.a3=sum(c((ad3-a.p3)^2)); PRESS.a4=sum(c((ad4-a.p4)^2)); PRESS.a5=sum(c((ad5-a.p5)^2))
    PRESS.d1=sum(c((dom1-d.p1)^2)); PRESS.d2=sum(c((dom2-d.p2)^2)); PRESS.d3=sum(c((dom3-d.p3)^2)); PRESS.d4=sum(c((dom4-d.p4)^2)); PRESS.d5=sum(c((dom5-d.p5)^2))  
    PRESS.g1=sum(c((g1-g.p1)^2)); PRESS.g2=sum(c((g2-g.p2)^2)); PRESS.g3=sum(c((g3-g.p3)^2)); PRESS.g4=sum(c((g4-g.p4)^2)); PRESS.g5=sum(c((g5-g.p5)^2))
    
    PRESS.matrix=rbind(c(PRESS.a1,PRESS.a2,PRESS.a3,PRESS.a4,PRESS.a5),
                       c(PRESS.d1,PRESS.d2,PRESS.d3,PRESS.d4,PRESS.d5),
                       c(PRESS.g1,PRESS.g2,PRESS.g3,PRESS.g4,PRESS.g5))
    
    colnames(PRESS.matrix) <-c("PlantHeight","EarHeight","EarLength","EarRowNumber","X20KernelWeight")
    rownames(PRESS.matrix) <-c("PRESS.a","PRESS.d","PRESS.g")
    
    print("PRESS values:",)
    print(PRESS.matrix) 
    
    print("-----------------------------------------------------------------------------") 
    
    results=list(results=list(results),vare=list(vare),vara=list(vara),vard=list(vard),varg=list(varg),core.matrix=list(core.matrix),cora.matrix=list(cora.matrix),cord.matrix=list(cord.matrix),corg.matrix=list(corg.matrix),PRESS.matrix=list(PRESS.matrix),dif=list(dif),iter=list(iter))  
    
    if (i==1) {
      save(results,file = "Results_MV-GBLUP-AD_h(0.3)_selec_markers_19-11-2014.RData")
    } 
    if (i==2) {
      save(results,file = "Results_MV-GBLUP-AD_h(0.5)_selec_markers_19-11-2014.RData")
    }
    if (i==3) {
      save(results,file = "Results_MV-GBLUP-AD_h(0.7)_selec_markers_19-11-2014.RData")
    }
    if (i==4) {
      save(results,file = "Results_MV-GBLUP-AD_h(Mixed)_selec_markers_19-11-2014.RData")
    }
    if (i==5) {
      save(results,file = "Results_MV-GBLUP-AD_h(Hist)_selec_markers_19-11-2014.RData")
    }
    
    if (iter==itermax) break
  }
  ######################################  
}

