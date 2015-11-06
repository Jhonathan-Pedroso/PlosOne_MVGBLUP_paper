############################################################################################################
######################################Building single-cross Hybrids#########################################

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

##Simulating the single-cross hybrids genotypes:

##All the phenotypic and genotypic data processed of the inbred lines coded above used to obtain the results of the manuscript can be load here:
setwd("/home/jhonathan/Documentos/MVGBLUP/Data/")
load("data_inbreds.RData")

gh.id=read.table("heterotic_group.txt",h=T)  #Inbred lines selected in the PCA analysis by visual graphical analysis of the inbred lines close to B73 and Mo17

id=as.matrix(data$`Inbred names`)

#Imputating molecular makers
require(rrBLUP)
Z.imp=as.numeric(Z)
Z.imp=matrix(Z.imp,nrow(Z),ncol(Z))
Z.imp=A.mat(Z.imp,min.MAF=F, max.missing=NULL,return.imputed=T)
Z.imp=matrix(round(Z.imp$imputed,2),nrow(Z),ncol(Z)) + 1  # Changing for the genotypic codification 2, 1 and 0
any(is.na(Z.imp)); Z=Z.imp; rm(Z.imp)

rownames(Z)<- data$`Inbred names`  #Labeling the rows of the marker matrix

Z.m=Z[as.factor(gh.id$Grupo1),]  # Separing he inbreds for the first hetoric group
Z.f=Z[as.factor(gh.id$Grupo2),]; rm(Z)  # Separing he inbreds for the second hetoric group

Z.h=matrix(0,(nrow(Z.m)*nrow(Z.f)),ncol(Z.f))  # Matrix of zeros that will be used to receive the built hybrids markers


#The construction of the hybrids in a 20x20 partial diallel crosses design
Z.f=Z.f[rep(1:20,20),]
Z.m=Z.m[rep(1:20,rep(20,20)),]

#Artificial crosses (from code lines 30 to 73) between the two heterotic groups. It is based on the expectation operation described in the manuscript
index0<- (Z.m==2 & Z.f==2)
Z.h[index0]<-2
index0<- (Z.m==2 & Z.f==1)
Z.h[index0]<-0.5*2+0.5*1
index0<- (Z.m==1 & Z.f==2)
Z.h[index0]<-0.5*2+0.5*1
index0<- (Z.m==0 & Z.f==2)
Z.h[index0]<-1
index0<- (Z.m==2 & Z.f==0)
Z.h[index0]<-1
index0<- (Z.m==0 & Z.f==0)
Z.h[index0]<-0
index0<- (Z.m>0 & Z.m<1 & Z.f>0 & Z.f<1)
Z.h[index0]<-(Z.m[index0]/1)*(0.5)*(Z.f[index0]/1)*0.5*2+((Z.m[index0]/1)*(0imputing gbs from whole genome sequecing data.5)*((Z.f[index0]/1)*(0.5)+(1-(Z.f[index0]/1))))*1+((Z.f[index0]/1)*(0.5)*((Z.m[index0]/1)*(0.5)+(1-(Z.m[index0]/1))))*1+(((Z.m[index0]/1)*(0.5)+(1-(Z.m[index0]/1)))*((Z.f[index0]/1)*(0.5)+(1-(Z.f[index0]/1))))*0
index0<- (Z.m>0 & Z.m<1 & Z.f>1 & Z.f<2)
Z.h[index0]<-(((Z.m[index0]/1)*(0.5))*((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*2)+(((Z.m[index0]/1)*(0.5))*(1-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5))*1)+(((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*((Z.m[index0]/1)*(0.5)+(1-(Z.m[index0]/1)))*1)+((Z.m[index0]/1)*(0.5)+(1-(Z.m[index0]/1)))*(1-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5))*0
index0<- (Z.m>1 & Z.m<2 & Z.f>1 & Z.f<2)
Z.h[index0]<-(((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*2)+(((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*(1-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5))*1)+(((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*(1-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5))*1)+((Z.m[index0]/1)*(0.5)+(1-(Z.m[index0]/1)))*((Z.f[index0]/1)*(0.5)+(1-(Z.f[index0]/1)))*0
index0<- (Z.m>1 & Z.m<2 & Z.f>0 & Z.f<1)
Z.h[index0]<-(((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*((Z.f[index0]/1)*(0.5))*2)+(((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*((Z.f[index0]/1)*(0.5)+(1-(Z.f[index0]/1)))*1)+((1-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5))*((Z.f[index0]/1)*(0.5))*1)+((1-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5))*((Z.f[index0]/1)*(0.5)+(1-(Z.f[index0]/1)))*0)
index0<- (Z.m>0 & Z.m<1 & Z.f==2)
Z.h[index0] <- ((Z.m[index0]/1)*(0.5))*2+((Z.m[index0]/1)*(0.5)+(1-(Z.m[index0]/1)))*1
index0<- (Z.m>1 & Z.m<2 & Z.f==2)
Z.h[index0] <- ((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*2+(1-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5))*1
index0<- (Z.m>0 & Z.m<1 & Z.f==1)
Z.h[index0]<-((Z.m[index0]/1)*(0.5))*0.5*2+((Z.m[index0]/1)*(0.5))*0.5*1+0.5*((Z.m[index0]/1)*(0.5)*1+(1-(Z.m[index0]/1)))*1
index0<- (Z.m>1 & Z.m<2 & Z.f==1)
Z.h[index0]<-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*0.5*2+((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*0.5*1+0.5*(1-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5))*1
index0<- (Z.m>0 & Z.m<1 & Z.f==0)
Z.h[index0]<-((Z.m[index0]/1)*(0.5))*1
index0<- (Z.m>1 & Z.m<2 & Z.f==0)
Z.h[index0]<-((Z.m[index0]/2)+(1-((Z.m[index0]/2)))*0.5)*1
index0<- (Z.m==2 & Z.f>0 & Z.f<1)
Z.h[index0]<- ((Z.f[index0]/1)*(0.5))*2+((Z.f[index0]/1)*(0.5)+(1-(Z.f[index0]/1)))*1
index0<- (Z.m==2 & Z.f>1 & Z.f<2)
Z.h[index0] <- ((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*2+(1-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5))*1
index0<- (Z.m==1 & Z.f>0 & Z.f<1)
Z.h[index0]<-((Z.f[index0]/1)*(0.5))*0.5*2+((Z.f[index0]/1)*(0.5))*0.5*1+0.5*((Z.f[index0]/1)*(0.5)*1+(1-(Z.f[index0]/1)))*1
index0<- (Z.m==1 & Z.f>1 & Z.f<2)
Z.h[index0]<-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*0.5*2+((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*0.5*1+0.5*(1-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5))*1
index0<- (Z.m==0 & Z.f>0 & Z.f<1)
Z.h[index0]<-((Z.f[index0]/1)*(0.5))*1
index0<- (Z.m==0 & Z.f>1 & Z.f<2)
Z.h[index0]<-((Z.f[index0]/2)+(1-((Z.f[index0]/2)))*0.5)*1

rm(gh.id,id,index0,Z.m,Z.f)

##Simulating the single-cross hybrids phenotypes in different heritabilities scenarios:

# The function below is used to deregressed the markers effects due to Bayes-B shrinkage specific variance procedure
correction_markers_g <- function (fmBB) {
  var.am.t=fmBB$ETA$Ad$varB
  var.dm.t=fmBB$ETA$Dom$varB
  ad=fmBB$ETA$Ad$b[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  var.a=fmBB$ETA$Ad$varB[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  dom=fmBB$ETA$Dom$b[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  var.d=fmBB$ETA$Dom$varB[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  hdj=0.7
  h.ca=matrix(0,100,1)
  for (i in 1:100) {
    var.a.dif.i <- var.a[-c(i)]
    sum.var.a.dif.i=sum(var.a.dif.i)+sum(var.d)
    var.a.t=sum(var.a)+sum(var.d)
    sum.var.t=sum(var.a)+sum(var.d)+fmBB$varE
    h.dif.i=sum.var.a.dif.i/sum.var.t
    h.t=var.a.t/sum.var.t
    h.ca[i,1]=hdj*(1-(h.dif.i/h.t))  
  }
  ad.new=ad/h.ca
  c=matrix(0,27000,1)
  c[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100]),] <-ad.new
  h.cd=matrix(0,100,1)
  for (i in 1:100) {
    var.d.dif.i <- var.d[-c(i)]
    sum.var.d.dif.i=sum(var.d.dif.i)+sum(var.a)
    sum.var.t=sum(var.a)+sum(var.d)+fmBB$varE
    var.d.t=sum(var.d)+sum(var.a)
    h.dif.i=sum.var.d.dif.i/sum.var.t
    h.t=var.d.t/sum.var.t
    h.cd[i,1]=hdj*(1-(h.dif.i/h.t))  
  }
  dom.new=dom/h.cd
  sum(h.cd)
  sum(rbind(h.ca,h.cd))
  d=matrix(0,27000,1)
  d[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100]),] <-dom.new
  ha=sum((var.a+var.d)/(var.a+var.d+fmBB$varE))
  d=d*ha; c=c*ha
  return(list(c=list(c),d=list(d)))
}   

#Building the phenotypic data using as reference the 0.3 heritability scenario
for (i in 1:5) {
  if (i==1) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_PH.RData") #Loading the results from the BayesB analysis
    y=fmBB$y; mu_PH=as.numeric(fmBB$mu)
  }
  if (i==2) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EH.RData")
    y=fmBB$y; mu_EH=as.numeric(fmBB$mu)
  }
  if (i==3) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EL.RData")
    y=fmBB$y; mu_EL=as.numeric(fmBB$mu)
  }
  if (i==4) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_ERN.RData")
    y=fmBB$y; mu_ERN=as.numeric(fmBB$mu)
  }
  if (i==5) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_KW.RData")
    y=fmBB$y; mu_KW=as.numeric(fmBB$mu)
  }
  
  mar=correction_markers_g(fmBB) # Deregress markers effects
  c=unlist(mar$c);d=unlist(mar$d)
  Z.h.a=Z.h[,(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]  #Separating the fraction of the genome with the 100 largest marker allele substitution effects
  Z.h.d=Z.h[,(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])] #Separating the fraction of the genome with the 100 largest marker dominance effects
  Wa=Wa.vit(Z.h.a)
  Wd=Wd.vit(Z.h.d)  
  c=c[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]  #Separating the 100 largest additive effects
  d=d[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])] #Separating the 100 largest dominant effects
  
  a <- c
  a <- Wa%*%a #Obtaining the parametric alelle substitution effects of the constructed single-cross hybrids 
  a.p=a
  Va=var(a) #Parametric additive genetic variance
  d <- Wd%*%d #Obtaining the parametric dominance deviations of the constructed single-cross hybrids
  d.p=d
  Vd=var(d) #Parametric dominance genetic variance
  Vg=Va+Vd #Parametric total genetic variance
  g=a+d
  h2 <- 0.3  #heritability
  Ve=(1-h2)/h2*Vg #Adjusting residual variance based on the heritability 
  y <- fmBB$mu + a + d + rnorm(nrow(Z.h),mean=0,sd=sqrt((1-h2)/h2*Vg)) #Simulating phenotypic effects
  comp=as.matrix(c(Ve,Va,Vd,Vg,a.p,d.p)) #saving parametric componentes 
  
  if (i==1) {
    y.h.PH.0.3=y; comp.PH.0.3=comp
  }
  if (i==2) {
    y.h.EH.0.3=y; comp.EH.0.3=comp
  }
  if (i==3) {
    y.h.EL.0.3=y; comp.EL.0.3=comp
  }
  if (i==4) {
    y.h.ERN.0.3=y; comp.ERN.0.3=comp
  }
  if (i==5) {
    y.h.KW.0.3=y; comp.KW.0.3=comp
  }  
}

#Building the phenotypic data using as reference the 0.5 heritability scenario
#The next three iterative processes below is the same as described above, but for the 0.5; 0,7 and historical heritability scenarios
for (i in 1:5) {
  if (i==1) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_PH.RData")
    y=fmBB$y; mu_PH=as.numeric(fmBB$mu)
  }
  if (i==2) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EH.RData")
    y=fmBB$y; mu_EH=as.numeric(fmBB$mu)
  }
  if (i==3) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EL.RData")
    y=fmBB$y; mu_EL=as.numeric(fmBB$mu)
  }
  if (i==4) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_ERN.RData")
    y=fmBB$y; mu_ERN=as.numeric(fmBB$mu)
  }
  if (i==5) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_KW.RData")
    y=fmBB$y; mu_KW=as.numeric(fmBB$mu)
  }
  
  mar=correction_markers_g(fmBB)
  c=unlist(mar$c);d=unlist(mar$d)
  Z.h.a=Z.h[,(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  Z.h.d=Z.h[,(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  Wa=Wa.vit(Z.h.a)
  Wd=Wd.vit(Z.h.d)  
  c=c[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  d=d[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  
  a <- c
  a <- Wa%*%a
  a.p=a
  Va=var(a)
  d <- Wd%*%d
  d.p=d
  Vd=var(d)
  Vg=Va+Vd
  g=a+d
  h2 <- 0.5  #heritability
  Ve=(1-h2)/h2*Vg
  y <- fmBB$mu + a + d + rnorm(nrow(Z.h),mean=0,sd=sqrt((1-h2)/h2*Vg))
  comp=as.matrix(c(Ve,Va,Vd,Vg,a.p,d.p))
  
  if (i==1) {
    y.h.PH.0.5=y; comp.PH.0.5=comp
  }
  if (i==2) {
    y.h.EH.0.5=y; comp.EH.0.5=comp
  }
  if (i==3) {
    y.h.EL.0.5=y; comp.EL.0.5=comp
  }
  if (i==4) {
    y.h.ERN.0.5=y; comp.ERN.0.5=comp
  }
  if (i==5) {
    y.h.KW.0.5=y; comp.KW.0.5=comp
  }  
}

#Building the phenotypic data using as reference the 0.7 heritability scenario
for (i in 1:5) {
  if (i==1) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_PH.RData")
    y=fmBB$y; mu_PH=as.numeric(fmBB$mu)
  }
  if (i==2) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EH.RData")
    y=fmBB$y; mu_EH=as.numeric(fmBB$mu)
  }
  if (i==3) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EL.RData")
    y=fmBB$y; mu_EL=as.numeric(fmBB$mu)
  }
  if (i==4) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_ERN.RData")
    y=fmBB$y; mu_ERN=as.numeric(fmBB$mu)
  }
  if (i==5) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_KW.RData")
    y=fmBB$y; mu_KW=as.numeric(fmBB$mu)
  }
  
  mar=correction_markers_g(fmBB)
  c=unlist(mar$c);d=unlist(mar$d)
  Z.h.a=Z.h[,(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  Z.h.d=Z.h[,(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  Wa=Wa.vit(Z.h.a)
  Wd=Wd.vit(Z.h.d)  
  c=c[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  d=d[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  
  a <- c
  a <- Wa%*%a
  a.p=a
  Va=var(a)
  d <- Wd%*%d
  d.p=d
  Vd=var(d)
  Vg=Va+Vd
  g=a+d
  h2 <- 0.7  #heritability
  Ve=(1-h2)/h2*Vg
  y <- fmBB$mu + a + d + rnorm(nrow(Z.h),mean=0,sd=sqrt((1-h2)/h2*Vg))
  comp=as.matrix(c(Ve,Va,Vd,Vg,a.p,d.p))
  
  if (i==1) {
    y.h.PH.0.7=y; comp.PH.0.7=comp
  }
  if (i==2) {
    y.h.EH.0.7=y; comp.EH.0.7=comp
  }
  if (i==3) {
    y.h.EL.0.7=y; comp.EL.0.7=comp
  }
  if (i==4) {
    y.h.ERN.0.7=y; comp.ERN.0.7=comp
  }
  if (i==5) {
    y.h.KW.0.7=y; comp.KW.0.7=comp
  }  
}

#Building the phenotypic data using as reference the historical heritability scenario
for (i in 1:5) {
  if (i==1) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_PH.RData")
    y=fmBB$y
    h2 <- 0.569; mu_PH=as.numeric(fmBB$mu)
  }
  if (i==2) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EH.RData")
    y=fmBB$y
    h2 <- 0.662; mu_EH=as.numeric(fmBB$mu)
  }
  if (i==3) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_EL.RData")
    y=fmBB$y
    h2 <- 0.381; mu_EL=as.numeric(fmBB$mu)
  }
  if (i==4) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_ERN.RData")
    y=fmBB$y
    h2 <- 0.57; mu_ERN=as.numeric(fmBB$mu)
  }
  if (i==5) {
    load("/home/jhonathan/Documentos/MVGBLUP/Data/fmBB_KW.RData")
    y=fmBB$y
    h2 <- 0.418; mu_KW=as.numeric(fmBB$mu)
  }
  
  mar=correction_markers_g(fmBB)
  c=unlist(mar$c);d=unlist(mar$d)
  Z.h.a=Z.h[,(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  Z.h.d=Z.h[,(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  Wa=Wa.vit(Z.h.a)
  Wd=Wd.vit(Z.h.d)  
  c=c[(order(fmBB$ETA$Ad$b,decreasing=T)[1:100])]
  d=d[(order(fmBB$ETA$Dom$b,decreasing=T)[1:100])]
  
  a <- c
  a <- Wa%*%a
  a.p=a
  Va=var(a)
  d <- Wd%*%d
  d.p=d
  Vd=var(d)
  Vg=Va+Vd
  g=a+d
  Ve=(1-h2)/h2*Vg
  y <- fmBB$mu + a + d + rnorm(nrow(Z.h),mean=0,sd=sqrt((1-h2)/h2*Vg))
  comp=as.matrix(c(Ve,Va,Vd,Vg,a.p,d.p))
  
  if (i==1) {
    y.h.PH.hist=y; comp.PH.hist=comp
  }
  if (i==2) {
    y.h.EH.hist=y; comp.EH.hist=comp
  }
  if (i==3) {
    y.h.EL.hist=y; comp.EL.hist=comp
  }
  if (i==4) {
    y.h.ERN.hist=y; comp.ERN.hist=comp
  }    
  if (i==5) {
    y.h.KW.hist=y; comp.KW.hist=comp
  }  
}

rm(Va,Vd,Ve,Vg,Wa,Wd,Z.h.a,Z.h.d,a,a.p,comp,d,d.p,g,y,c,fmBB,h2,i,mar,data,correction_markers_g,Wa.vit,Wd.vit)

save(list = ls(all = TRUE),file = "data_hybrids.RData") #Saving single-cross hybrids data

############################################################################################################
