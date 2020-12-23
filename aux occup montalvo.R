tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#------------------------------------------------
sample.z=function(xmat.occ,betas,nloc,nspp,nrep,y,xmat.det,gammas){
  media.occ=xmat.occ%*%betas
  prob.occ=pnorm(media.occ)
  prob.det=matrix(1,nloc,nspp)
  for (i in 1:nrep){
    media.det=xmat.det[,,i]%*%gammas
    tmp=1-pnorm(media.det)
    prob.det=prob.det*tmp
  }
  tmp=prob.det*prob.occ
  prob.fim=tmp/(tmp+(1-prob.occ))
  
  #sample z
  z=matrix(1,nloc,nspp)
  max.y=apply(y,c(1,2),max)
  cond=max.y==0
  z[cond]=rbinom(sum(cond),size=1,prob=prob.fim[cond])
  z
}
#----------------------------------------------
sample.zstar=function(z,xmat.occ,betas,nloc,nspp){
  media=xmat.occ%*%betas
  lo.mat=matrix(ifelse(z==1,0,-1000),nloc,nspp)
  hi.mat=matrix(ifelse(z==0,0, 1000),nloc,nspp)
  zstar=tnorm(n=nloc*nspp,lo=lo.mat,hi=hi.mat,mu=media,sig=1)
  matrix(zstar,nloc,nspp)
}
#--------------------------------------------
sample.ystar=function(nrep,xmat.det,gammas,y,nnloc,nspp){
  ystar=y
  for (j in 1:nrep){
    media=xmat.det[,,j]%*%gammas
    lo.mat=matrix(ifelse(y[,,j]==1,0,-1000),nloc,nspp)
    hi.mat=matrix(ifelse(y[,,j]==0,0, 1000),nloc,nspp)
    ystar[,,j]=tnorm(n=nloc*nspp,lo=lo.mat,hi=hi.mat,mu=media,sig=1)
  }
  ystar
}
#--------------------------------------------
sample.betas1=function(m.betas1,tau2.betas1,betas2,
                      xmat.occ1,xtx1,zss,nspp,nparam.occ1){
  prec1=diag(1/tau2.betas1)
  if (nparam.occ1==1) prec1=1/tau2.betas1
  var1=solve(xtx1+prec1)
  
  betas1=matrix(NA,nparam.occ1,nspp)
  for (i in 1:nspp){
    pmedia=t(xmat.occ1)%*%zss[,i]+prec1%*%m.betas1
    betas1[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
  }
  betas1
}
#----------------------------------------
sample.gammas=function(ystar,xmat.det,z,m.gamma,tau2.gamma,nparam.det,nspp,nrep){
  invTau=diag(1/tau2.gamma)
  gammas=matrix(NA,nparam.det,nspp)
  for (i in 1:nspp){
    cond=z[,i]==1
    for (j in 1:nrep){
      if (j==1) {
        xmat.det1=xmat.det[cond,,j]
        ystar1=ystar[cond,i,j]
      }
      if (j> 1) {
        xmat.det1=rbind(xmat.det1,xmat.det[cond,,j])
        ystar1=c(ystar1,ystar[cond,i,j])
      }
    }
    xtx=t(xmat.det1)%*%xmat.det1    
    var1=solve(xtx+invTau)
    pmedia1=t(xmat.det1)%*%ystar1+invTau%*%m.gamma
    gammas[,i]=rmvnorm(1,var1%*%pmedia1,var1)
  }
  gammas
}
#--------------------------------------------
sample.m.gamma=function(gammas,tau2.gamma,nspp,nparam.det){
  soma.gammas=rowSums(gammas)
  invTau=diag(1/tau2.gamma)
  prec=nspp*invTau + diag(rep(1/100,nparam.det))
  var1=solve(prec)
  pmedia=invTau%*%soma.gammas
  t(rmvnorm(1,mean=var1%*%pmedia,var1))
}
#--------------------------------------------
sample.tau2.gamma=function(gammas,m.gamma,nspp,tau2.a,tau2.b,nparam.det){
  a1=(nspp+2*tau2.a)/2
  m.gamma.mat=matrix(m.gamma,nparam.det,nspp)
  err2=rowSums((gammas-m.gamma.mat)^2)
  1/rgamma(nparam.det,a1,(err2/2)+tau2.b)
}
#----------------------------------------
sample.betas2=function(zss,w,nparam.occ2,ngr,xtx2,nloc,xmat.occ2){
  betas2=matrix(NA,nparam.occ2,ngr)
  for (i in 1:ngr){
    cond=w==i
    nw=sum(cond)
    prec1=(nw*xtx2)+diag(rep(1/100,nparam.occ2))
    var1=solve(prec1)
    if (nw==0) soma.zss=matrix(0,nloc,1)
    if (nw==1) soma.zss=zss[,cond]
    if (nw >1) soma.zss=matrix(rowSums(zss[,cond]),nloc,1)
    pmedia=t(xmat.occ2)%*%soma.zss
    betas2[,i]=rmvnorm(1,var1%*%pmedia,var1)
  }
  betas2
}
#-----------------------------------
sample.w=function(betas2,ltheta,w,ngr,nparam.occ2,xmat.occ2,zss,xtx2){
  
  #pre-compute some stuff
  media=xmat.occ2%*%betas2
  prec1=xtx2+diag(rep(1/100,nparam.occ2))
  var1=solve(prec1)
  
  tmp=-(nparam.occ2/2)*log(100)
  tmp1=-(1/2)*determinant(prec1,logarithm=T)$modulus[[1]]
  p1=tmp+tmp1
  for (i in 1:nspp){
    zss.mat=matrix(zss[,i],nloc,ngr)
    err2=colSums((zss.mat-media)^2)
    lprob=-(1/2)*err2+ltheta

    #do we need a new group?    
    max.w=max(w)
    if (max.w<ngr){
      lprob=lprob[1:max.w]
      ind=max.w+1
      
      #probability of new group
      mu=var1%*%t(xmat.occ2)%*%zss[,i]
      sum.zss2=sum(zss[,i]^2)
      p2=(-t(mu)%*%prec1%*%mu+sum.zss2)
      lprob=c(lprob,p1+(-1/2)*p2+ltheta[ind])
    }
    tmp=lprob-max(lprob)
    tmp=exp(tmp)
    prob=tmp/sum(tmp)
    ind=rmultinom(1,size=1,prob=prob)
    w[i]=which(ind==1)
  }
  w
}
#----------------------------------------------
sample.theta=function(gamma1,w,ngr){
  nk=rep(0,ngr)
  tmp=table(w)
  nk[as.numeric(names(tmp))]=tmp
  theta=v=rep(NA,ngr)
  aux=1
  for(i in 1:(ngr-1)){
    nk.maior=nk[(i+1):ngr]
    v[i]=rbeta(1,nk[i]+1,sum(nk.maior)+gamma1)
    theta[i]=v[i]*aux
    aux=aux*(1-v[i])
  }
  theta[ngr]=aux
  theta
}
#--------------------------------------------
get.llk=function(nloc,nspp,betas,xmat.occ,y,gammas,xmat.det){
  #get occupancy info
  media.occ=xmat.occ%*%betas
  Phi.occ=pnorm(media.occ)
  lPhi.occ=log(Phi.occ)
  One_Phi.occ=1-Phi.occ

  #get detection info
  lOne_Phi.det=lPhi.det=One_Phi.det=Phi.det=array(NA,dim=c(nloc,nspp,nrep))
  for (i in 1:nrep){
    media.det=xmat.det[,,i]%*%gammas
    Phi.det[,,i]=pnorm(media.det)
    lPhi.det[,,i]=log(Phi.det[,,i])
    One_Phi.det[,,i]=1-Phi.det[,,i]
    lOne_Phi.det[,,i]=log(One_Phi.det[,,i])
  }
  
  #get llk
  fim=rep(NA,nspp)
  for (i in 1:nspp){
    lprob=rep(NA,nloc)
    y1=y[,i,]
    tmp=rowSums(y1)
    
    #all observations are equal to zero
    cond=tmp==0
    # lp1=rowSums(log(One_Phi.det[cond,i,]))+log(Phi.occ[cond,i])
    p1=apply(One_Phi.det[cond,i,],1,prod)*Phi.occ[cond,i]
    p2=One_Phi.occ[cond,i]
    lprob[cond]=log(p1+p2)
    
    #at least one observations is equal to 1
    lp1=y1[!cond,]*lPhi.det[!cond,i,]+(1-y1[!cond,])*lOne_Phi.det[!cond,i,]
    lp2=rowSums(lp1)+lPhi.occ[!cond,i]
    lprob[!cond]=lp2
    fim[i]=sum(lprob)
  }
  sum(fim)
}
#--------------------------------
sample.m.betas1=function(nspp,betas1,tau2.betas1,nparam.occ1){
  prec=(nspp/tau2.betas1)+(1/100)
  var1=1/prec
  if (nparam.occ1 >1) pmedia=rowSums(betas1)/tau2.betas1
  if (nparam.occ1==1) pmedia=sum(betas1)/tau2.betas1
  rnorm(nparam.occ1,mean=var1*pmedia,sd=sqrt(var1))
}
#----------------------------------------
sample.tau2.betas1=function(nspp,betas1,m.betas1,tau2.a,
                            tau2.b,nparam.occ1){
  a1=(nspp/2)+tau2.a
  m.betas1.mat=matrix(m.betas1,nparam.occ1,nspp)
  if (nparam.occ1 >1) err2=rowSums((betas1-m.betas1.mat)^2)
  if (nparam.occ1==1) err2=sum((betas1-m.betas1.mat)^2)
  b1=tau2.b+(err2/2)
  1/rgamma(nparam.occ1,a1,b1)
}
