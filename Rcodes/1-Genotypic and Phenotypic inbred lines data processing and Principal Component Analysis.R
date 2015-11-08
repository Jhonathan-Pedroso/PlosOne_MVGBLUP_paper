############################################################################################################

#All the genotypic data ZeaGBSv1.0 (Romay et al. 2003) are available at http://www.panzea.org/#!genotypes/cctl
#This data set was subdivided into 20 parts linearly along the chromessomes (2 parts each chromossome) (softwere TASSEL 4.0) in order to turn more 
#efficient the usage of the computer memory. 

#The final imputed inbred line molecular markers matrix and all phenotyic data (used to obtain the results of the manuscript) can be loaded 
#at the begining of the code with the tittle: Inbred lines analysis with the BayesB model

#All the code used to process the genotypic and phenotypic data of the inbred lines used to obtain the results found in the manuscript can be accessed bellow:

############################################################################################################
#####Systematic sampling of 27.000 markers uniformly distributed throughout the 10 maize linkage groups#####
                        
setwd("/home/jhonathan/Documentos/MVGBLUP/Data")

for (i in 1:20) {
  
  if (i==1) {
    sampl=sample(rank(runif(25)),1) #Drawing of the first marker in the range between the 1th and 25th marks
    M=read.table("ch1_1.txt",h=F)  #Pushing half of the molecular data of chromossome 1
    sequence=round((seq(from=sampl,to=ncol(M),length.out=1350)),0)  #Defining the markers position throughout the sequence of the first chromosome
    M.f=M[,sequence] #Pushing only the markers present in the range
  }
  if (i==2) {
    M=read.table("ch1_2.txt",h=F)
  }
  if (i==3) {
    M=read.table("ch2_1.txt",h=F)
  }
  if (i==4) {
    M=read.table("ch2_2.txt",h=F)
  }
  if (i==5) {
    M=read.table("ch3_1.txt",h=F)
  }
  if (i==6) {
    M=read.table("ch3_2.txt",h=F)
  }
  if (i==7) {
    M=read.table("ch4_1.txt",h=F)
  }
  if (i==8) {
    M=read.table("ch4_2.txt",h=F)
  }
  if (i==9) {
    M=read.table("ch5_1.txt",h=F)
  }
  if (i==10) {
    M=read.table("ch5_2.txt",h=F)
  }
  if (i==11) {
    M=read.table("ch6_1.txt",h=F)
  }
  if (i==12) {
    M=read.table("ch6_2.txt",h=F)
  }
  if (i==13) {
    M=read.table("ch7_1.txt",h=F)
  }
  if (i==14) {
    M=read.table("ch7_2.txt",h=F)
  }
  if (i==15) {
    M=read.table("ch8_1.txt",h=F)
  }
  if (i==16) {
    M=read.table("ch8_2.txt",h=F)
  }
  if (i==17) {
    M=read.table("ch9_1.txt",h=F)
  }
  if (i==18) {
    M=read.table("ch9_2.txt",h=F)
  }
  if (i==19) {
    M=read.table("ch10_1.txt",h=F)
  }
  if (i==20) {
    M=read.table("ch10_2.txt",h=F)
  }
  
  if (i!=1) {
    sequence=round((seq(from=1,to=ncol(M),length.out=1350)),0)  #Defining the markers position throughout the sequence of the other chromosomes
    M=M[,sequence]
    M.f=cbind(M.f,M) #Binding the columns of the markers matrix sampled until to obtain in the end of the iterative process the final matrix
  }
  
}

M=M.f; rm(M.f,sequence,sampl,i)  #Final markers matrix 

############################################################################################################

############################################################################################################
#####################Changing the markers codification from letters to numbers (-1,0,1).#####################

# Letters are encoding unphased diploid genotypes in the following manner: 
#A=A/A=-1; C=C/C=1; G=G/G=1; T=T/T=-1; K=G/T=0; M=A/C=0; R=A/G=0; S=C/G=1; W=A/T=-1; Y=C/T=0; N=missing;

#In the iterative process bellow will be changed the letters for the codification -1, 0 and 1 in the markers matrix:

M=as.matrix(M)
M=M[-c(1),]
index <- M=="A" |  M=="T" |  M=="W"  
M[index]<-(-1)
index <- M=="C" | M=="G" | M=="S"
M[index]<-(1)
index <- M=="K" | M=="M" | M=="R" | M=="Y"
M[index]<-(0)
index <- M=="N"
M[index]<-NA
Z=M; rm(index,M)
############################################################################################################

############################################################################################################
###########################Processing phenotypic and genotypic inbred lines data############################

#All the phenotypic data (traitMatrix_maize282NAM_v15-130212.txt) were download at http://www.panzea.org/#!phenotypes/c1m50. 
#This data is from the NAM and Maize 282 association collected between 2006 and 2009. 
#We filtered the data at TASSEL softwere TASSEL 4.0 related to the traits: Plant Height, Ear Insertion Height, Ear Length, Numbers of Rows per Ear and Kernel Weight.
#The data was colleted in assays conducted in 2006 during the summer harvest in the city of Aurora, New York.
#Trait descriptions and interpretations of the environment (location/year) codes are available at the last link above.

setwd("/home/jhonathan/Documentos/MVGBLUP/Data")
data=read.table("pheno_lines_filtered.txt",na.strings="",sep="\t") #Pushing the phenotypic data filtered
data=data[-c(1:3),]
colnames(data)<-c("Inbred names","PH","EH","EL","ERN","KW"); data=data[-c(1),]
ID1=as.matrix(data[,1]) #Identification of the lines from the phenotypic data filtered

ID2=read.table("markers_lines_id.txt",na.strings="",sep="\t") #Identification the lines from the genotypic data processed
colnames(ID2) <-c("ID", "Inbred names markers")

lines_select=merge(x=ID2,y=ID1,by.x="ID",by.y="V1") #Separating the inbred lines names common to both phenotypic and genotypic data

data=merge(x=data,y=lines_select,by.x="Inbred names", by.y="ID") #Separating the phenotypic inbred lines data common to both phenotypic and genotypic data
data=data[,c(1,7,2,3,4,5,6)]; #Reordering data columns
colnames(data) <- c("ID", "Inbred names", "PH", "EH", "EL", "ERN", "KW")

Z=data.frame(Z)
Z=cbind(ID2,Z)  
Z=Z[ID2[,1],] #Separating the genotypic inbred lines data common to both phenotypic and genotypic data

Z=Z[lines_select[,1],] 
rm(ID1,ID2,lines_select)
n2=ncol(Z)-2
Z=Z[,3:(n2+2)]  #Separating genotypic data from the IDs
Z=as.matrix(Z); rm(n2)

save(list = ls(all = TRUE),file = "data_inbreds.RData") #Saving the inbred lines data

############################################################################################################

############################################################################################################
##########################################PCA Analysis######################################################

##All the phenotypic and genotypic data processed above used to obtain the results of the manuscript can be download here:
setwd("/home/jhonathan/Documentos/MVGBLUP/Data")
load("data_inbreds.RData")

# Additive relationship matrix function
A.vit= function (K) {
  freq=colMeans(K)/2
  
  Fa=t(matrix(freq,ncol(K),nrow(K)))
  
  W=K-2*Fa
  cor=2*sum(freq*(1-freq))
  Av=W%*%t(W)/cor
  return (Av)
}

##Imputing inbred lines molecular marker data
require(rrBLUP)
Z.imp=as.numeric(Z)
Z.imp=matrix(Z.imp,nrow(Z),ncol(Z))
Z.imp=A.mat(Z.imp,min.MAF=F, max.missing=NULL,return.imputed=T)
Z.imp=matrix(round(Z.imp$imputed,2),nrow(Z),ncol(Z)) + 1  # Changing for the genotypic codification 2, 1 and 0
any(is.na(Z.imp))

A=A.vit(Z.imp) # Building additive relationship matrix function

dec=eigen(t(A))  # Performing spectral decompositon 
U=dec$vector; l=dec$values; v=dec$v  # Sepating eigenvectors and eigenvalues
pc1=U[,1]*(l[1]); pc2=U[,2]*(l[2]) #Obtaining the firt and second principal component
pc1=as.matrix(pc1); pc2= as.matrix(pc2)
rownames(pc1) <-data$`Inbred names`; rownames(pc2) <-data$`Inbred names` #Labeling the principal components

setwd("/home/jhonathan/Documentos/MVGBLUP/Data/")
gh.id=read.table("heterotic_group.txt",h=T)  #Inbred lines selected in the PCA analysis by visual graphical analysis of the inbred lines close to B73 and Mo17

##The content below is the code to obtain the plot showing the chosen and unchosen inbred lines as present in the manuscript 

g.m=as.factor(gh.id$Grupo1) 
g.f=as.factor(gh.id$Grupo2)

test=data$`Inbred names`[-c(g.m)];test=test[-c(g.f)]
test=as.matrix(test)
rownames(A) <- data$`Inbred names`
t1=A[test,]

pc1.t=pc1;pc1.t=pc1.t[-c(g.m)];pc1.t=pc1.t[-c(g.m)];
pc2.t=pc2;pc2.t=pc2.t[-c(g.m)];pc2.t=pc2.t[-c(g.m)];

setwd("/home/jhonathan/Documentos/MVGBLUP/Results")

val=cbind(as.matrix(data$`Inbred names`),1:length(data$`Inbred names`))
cod=cbind(as.numeric(val[,2]),as.matrix(c(pc1)),as.matrix(c(pc2)))
colnames(cod)<-c("lines","pc1","pc2")

cod1 <- cod[,2]<as.numeric(-2.5) 
b73g1<-as.matrix(val[,2][cod1])
b73g1 <- merge(cod,b73g1, by.x="lines", by.y="V1")
cod1 <- cod[,2]>as.numeric(-7.5) #b73
b73t<-as.matrix(val[,2][cod1])
b73g1 <- merge(b73g1,b73t, by.x="lines", by.y="V1")

p1 <- b73g1[,3]>10
p1 <- as.matrix(b73g1[,1][p1])
p1 <-merge(b73g1,p1, by.x="lines", by.y="V1")

g1 <-merge(val,p1, by.x="V2", by.y="lines")
g1 <-merge(g1,as.matrix(g.m), by.x="V1", by.y="V1"); g1=cbind(as.matrix(1:nrow(g1)),g1)
colnames(g1)=1:ncol(g1)
vec=as.matrix(c(1,3,4,6,7,8,10,45,46,47,50,51,52,53,55,56,58,59,62,63))
g1=merge(g1,vec, by.x="1", by.y="V1")

jpeg("Plot_PCA_lines.jpeg", width = 15, height = 10, units = 'in', res = 350)
plot(pc1,pc2,pch=16, type="p",col='gray81',xlab='pc1',ylab='pc2',font=2,font.lab=2,family="serif")
text(g1[,4], g1[,5], g1[,2], cex=0.6, pos=4,col="firebrick2",family="serif")

cod1 <- cod[,2]>as.numeric(-10) 
Mo17g1<-as.matrix(val[,2][cod1])
Mo17g1 <- merge(cod,Mo17g1, by.x="lines", by.y="V1")
cod1 <- cod[,2]<as.numeric(10) #Mo17
Mo17t<-as.matrix(val[,2][cod1])
Mo17g1 <- merge(Mo17g1,Mo17t, by.x="lines", by.y="V1")

p2 <- Mo17g1[,3]>-10
p2 <- as.matrix(Mo17g1[,1][p2])
p2 <-merge(Mo17g1,p2, by.x="lines", by.y="V1")

p3 <- p2[,3]<10
p3 <- as.matrix(p2[,1][p3])
p3 <-merge(p2,p3, by.x="lines", by.y="V1")
g2 <-merge(val,p3, by.x="V2", by.y="lines")
g2 <-merge(g2,as.matrix(g.f), by.x="V1", by.y="V1"); g2=cbind(as.matrix(1:nrow(g2)),g2)
colnames(g2)=1:ncol(g2)
vec=as.matrix(c(1,3,4,6,7,8,12,19,20,21,22,24,31,33,34,35,36,37,38,40))
g2=merge(g2,vec, by.x="1", by.y="V1")
g2=g2[-c(3),]
text(g2[,4], g2[,5], g2[,2], cex=0.6, pos=4,col="deepskyblue4",family="serif")
text(pc1[2534,],pc2[2534,],"Ames19008", cex=0.6, pos=4,col="deepskyblue4",family="serif")
dev.off()
rm(A,Mo17g1,Mo17t,U,b73g1,b73t,cod,g1,g2,gh.id,pc1,pc2,t1,test,val,vec,cod1,dec,g.f,g.m,l,pc1.t,pc2.t,v,p1,p2,p3,Z.imp)

############################################################################################################

