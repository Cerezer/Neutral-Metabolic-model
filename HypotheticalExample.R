######################################################
# Script used to simulate the dynamics of neutral and metabolic processes with hypothetical data.
# Estimate the likelihood of combined parameters.
# Created by FO Cerezer in August 2017

## packages required
library(reshape)
library(rgdal)
library(raster)
library(maptools)

## Creating a matrix that represents a hypothetical community
AbundCom<-matrix(c(10,20,30),3,2) 
colnames(AbundCom)<-1:ncol(AbundCom)

# Convert to a long table with 3 columns (i.e. number of the community, the species present in the community, and its abundance)
longCom<-melt(AbundCom) 

## Identifies in which community the individual is present
comInd<-rep(longCom$X1,longCom$value)

## Identifies the species of the individual 
comm1<-rep(longCom$X2,longCom$value) 

## Number of individuals per community
nind<-tabulate(comInd,nbins = nrow(AbundCom))

## Community indexer
coms<-sort(unique(comInd))

# Determine the position of each individual (binned within communities - grid cells)

Latitude<- 1:4
Longitude<- 1:4

## Plotting the hypothetical communities
plot(x<-rnorm(sum(nind),Longitude[comInd],0.1),y<-rnorm(sum(nind),Latitude[comInd],0.1),
pch=21,bg=comm1, ylab = "Latitude", xlab = "Longitude")

## Creates a geographic distance matrix between pairs of communities
Coord<- data.frame(Latitude,Longitude)

D<- as.matrix(dist(Coord))  

## Define model parameters

#z # Rate in which dispersal decays with geodistance
 
temp<-c(10,14,8)+273.15 # Hypothetical temperature for each community in Kelvin

vlogit # Basal speciation rate transformed into exponential

E<-0.65 # Activation energy of metabolic processes

K<- 8.617*10^-5 # Boltzmann's constant

b # The intensity of the temperature effect on the basal speciation rate

## Grid Aproximation function 
# This approach was used to investigate the probability of different values of the parameters (z, vlogit, and b) explain the observed data of the communities

## Creating vectors with range values for the parameters of the neutral-metabolic model

bvec<-seq(0,1,length.out = 2) # Parameter b ranges from 0 to 1. Only two values are used as an illustration

zvec<- seq(0,0.1,length.out = 2) #Parameter z includes values with minimal to maximum dispersal between communities. Only two values are used as an illustration

vlogitvec<- seq(-5,-2.44,length.out = 2) # Parameter v ranges from 0 (no speciation) to 1 (with speciation). Vlogit values represent the values of v transformed into exponential (see the NeutralTempModel function below). Only two values are used as an illustration

## Data frame with all combinations of parameters previously created
Pars<-expand.grid(bvec,zvec,vlogitvec)
Pars2<- data.frame(t(Pars))

## Simulating the dynamics of neutral and metabolic processes with the different combinations of values for the parameters z, vlogit, and b. 

# Packages required
library(snow)
library(Rmpi)

## Load functions
load(NeutralTempModel)
load(GridApproxFunct)


{
  clusterGridAp<-makeSOCKcluster(3)
  clusterExport(clusterGridAp,c("temp","comInd","coms","nind","comm1","D","GridApproxFunctWrap","NeutralTempModel","GridApproxFunct"))
  CluAppGridAp<-clusterApply(clusterGridAp,Pars2, GridApproxFunctWrap, temp = temp,comInd = comInd,coms = coms,nind = nind,comm1 = comm1,D = D,E = 0.65,K = 8.617*10^-5,ngen = 25,rep = 100)
  stopCluster(CluAppGridAp)
  
}

## Importing the files with the different parameter combinations rotated
## Averaging the 100 repetitions for each combination 
## Calculating the log-likelihood for each combination

## Finds the list of files containing the parameter combinations and saves it in the ParsOpen object   
ParsOpen<- list.dirs(".../")

## Saves the observed species richness in Sobs object 
Sobs<-tapply(comm1,comInd,function(x)length(unique(x))) 

## Saves the simulated species richness (average of 100 rep/parms) in the ParsLambda object
ParsLambda<-matrix(NA, nrow = length(ParsOpen), ncol = length(coms))

## Saves the loglikelihood values for each combination of parameters in a matrix
loglik<- {}

## For used to access the folder with all the repetitions of each combination of parameters

for (i in 1:length(ParsOpen)) {
  
  ParsOpenFiles<-list.files(ParsOpen[i],full.names = TRUE) 
  
  lambda<-matrix(NA,length(ParsOpenFiles),length(Sobs))
  
  for (k in 1:length(ParsOpenFiles)) {
    
    DataSpID<-readLines(ParsOpenFiles[k])
    
    lambda[k,]<-tapply(DataSpID,comInd,function(x)length(unique(x))) # Saves each repetition on a line of the lambda matrix to calculate the mean
    
    print(c(i,k))
    
  } 
  
  lambdaMean<-colMeans(lambda)
  ParsLambda[i,]<-lambdaMean
  
  # Estimating the likelihood of the predicted model with observed data
  loglik[i]<-sum(dpois(Sobs,lambdaMean,log = TRUE))
  
}
