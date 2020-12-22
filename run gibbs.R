# rm(list=ls())
library('mvtnorm')
set.seed(2)

#basic settings
nrep=5

#get functions
setwd('U:\\GIT_models\\occup_montalvo_SAM')
source('aux occup montalvo.R')
source('gibbs_occup_montalvo.R')

#get design matrix for occupancy
setwd('U:\\GIT_models\\occup_montalvo_SAM\\fake data')
xmat.occ=data.matrix(read.csv('fake data xmat occ.csv',as.is=T))
ind.reff=1:2
ind.gr=3:ncol(xmat.occ)
xmat.occ1=xmat.occ[,ind.reff]
xmat.occ2=xmat.occ[,ind.gr]
nloc=nrow(xmat.occ)

#get data
tmp=read.csv('fake data y.csv',as.is=T)
nspp=nrow(tmp)/(nrep*nloc); nspp
y=array(tmp$V1,dim=c(nloc,nspp,nrep))

#design matrix for detection
tmp=read.csv('fake data xmat det.csv',as.is=T)
nparam.det=nrow(tmp)/(nloc*nrep); nparam.det
xmat.det=array(tmp$V1,dim=c(nloc,nparam.det,nrep))

#basic settings
ngr=10
ngibbs=1000
nburn=ngibbs/2

#priors
tau2.a=0.1; tau2.b=0.1
gamma1=0.1

#run gibbs
mod1=gibbs_occup(y=y,xmat.occ1=xmat.occ1,xmat.occ2=xmat.occ2,
                 xmat.det=xmat.det,ngr=ngr,
                 tau2.a=tau2.a,tau2.b=tau2.b,gamma1=gamma1,
                 ngibbs=ngibbs,nburn=nburn)