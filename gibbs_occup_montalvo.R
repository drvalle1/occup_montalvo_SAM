gibbs_occup=function(y,xmat.occ1,xmat.occ2, xmat.det,ngr,tau2.a,tau2.b,
                     gamma1,ngibbs,nburn){

  #pre-compute some stuff
  xtx1=t(xmat.occ1)%*%xmat.occ1
  xtx2=t(xmat.occ2)%*%xmat.occ2
  xmat.occ=cbind(xmat.occ1,xmat.occ2)
  
  #initial values
  nspp=ncol(y[,,1])
  nloc=nrow(xmat.occ1)
  nparam.occ1=ncol(xmat.occ1)
  nparam.occ2=ncol(xmat.occ2)
  nparam.det=ncol(xmat.det[,,1])
  ystar=y
  ystar=ifelse(y==1,1,-1)
  m.betas1=rep(0,nparam.occ1)
  tau2.betas1=rep(1,nparam.occ1)
  w=sample(1:ngr,size=nspp,replace=T)
  betas1=matrix(0,nparam.occ1,nspp)
  betas2=matrix(0,nparam.occ2,ngr)
  theta=rep(1/ngr,ngr)
  z=apply(y,c(1,2),max)
  zstar=matrix(ifelse(z==1,1,-1),nloc,nspp)
  m.gamma=rep(0,nparam.det)
  tau2.gamma=rep(1,nparam.det)
  gammas=matrix(0,nparam.det,nspp)

  #MCMC settings
  store.betas2=matrix(NA,ngibbs,nparam.occ2*ngr)
  store.betas1=matrix(NA,ngibbs,nparam.occ1*nspp)
  store.m.betas1=matrix(NA,ngibbs,nparam.occ1)
  store.tau2.betas1=matrix(NA,ngibbs,nparam.occ1)
  store.gammas=matrix(NA,ngibbs,nparam.det*nspp)  
  store.m.gamma=matrix(NA,ngibbs,nparam.det)
  store.tau2.gamma=matrix(NA,ngibbs,nparam.det)
  store.w=matrix(NA,ngibbs,nspp)
  store.theta=matrix(NA,ngibbs,ngr)
  store.llk=matrix(NA,ngibbs,1)
  
  options(warn=2)

  #start gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    print(table(w))
    
    #update latent variables
    betas=rbind(betas1,betas2[,w])
    z=sample.z(xmat.occ=xmat.occ,betas=betas,nloc=nloc,
               nspp=nspp,nrep=nrep,y=y,xmat.det=xmat.det,gammas=gammas)
    # z=z.true
    zstar=sample.zstar(z=z,xmat.occ=xmat.occ,betas=betas,nloc=nloc,nspp=nspp)
    # zstar=zstar.true
    
    ystar=sample.ystar(nrep=nrep,xmat.det=xmat.det,gammas=gammas,
                       y=y,nspp=nspp)
    # ystar=ystar.true

    #sample betas1 and associated prior parameters
    zss=zstar-(xmat.occ2%*%betas2)[,w]
    betas1=sample.betas1(m.betas1=m.betas1,tau2.betas1=tau2.betas1,
                         betas2=betas2,xmat.occ1=xmat.occ1,xtx1=xtx1,
                         zss=zss,nspp=nspp,nparam.occ1=nparam.occ1)
    
    m.betas1=sample.m.betas1(nspp=nspp,betas1=betas1,tau2.betas1=tau2.betas1,
                             nparam.occ1=nparam.occ1)
    # m.betas1=m.betas1.true
    tau2.betas1=sample.tau2.betas1(nspp=nspp,betas1=betas1,
                                   nparam.occ1=nparam.occ1,
                                   m.betas1=m.betas1,tau2.a=tau2.a,tau2.b=tau2.b)
    
    #sample betas2
    zss=zstar-xmat.occ1%*%betas1
    betas2=sample.betas2(zss=zss,w=w,nparam.occ2=nparam.occ2,ngr=ngr,
                         xtx2=xtx2,nloc=nloc,xmat.occ2=xmat.occ2)

    #update gammas and associated prior parameters
    gammas=sample.gammas(ystar=ystar,xmat.det=xmat.det,z=z,m.gamma=m.gamma,tau2.gamma=tau2.gamma,
                         nparam.det=nparam.det,nspp=nspp,nrep=nrep)
    m.gamma=sample.m.gamma(gammas=gammas,tau2.gamma=tau2.gamma,nspp=nspp,nparam.det=nparam.det)
    tau2.gamma=sample.tau2.gamma(gammas=gammas,m=m.gamma,nspp=nspp,tau2.a=tau2.a,
                                 tau2.b=tau2.b,nparam.det=nparam.det)
    
    #sample other parameters
    w=sample.w(betas2=betas2,ltheta=log(theta),w=w,
               ngr=ngr,nparam.occ2=nparam.occ2,xmat.occ2=xmat.occ2,
               zss=zss,xtx2=xtx2)
    # w=w.true
    theta=sample.theta(gamma1=gamma1,w=w,ngr=ngr)
    # theta=rep(1/ngr,ngr)

    betas=rbind(betas1,betas2[,w])
    llk=get.llk(nloc=nloc,nspp=nspp,betas=betas,
                xmat.occ=xmat.occ,xmat.det=xmat.det,y=y,gammas=gammas)  
      
    #re-order w from time to time
    if (i<nburn & i%%50==0){
      k=table(w)
      k1=rep(0,ngr)
      k1[as.numeric(names(k))]=k
      ind=order(k1,decreasing=T)
      
      #if re-ordering needs to be done
      if (sum(ind!=1:ngr)>0) {
        theta=theta[ind]
        betas2=betas2[,ind]
        if (nparam.occ2==1) betas2=matrix(betas2,1,ngr)
        wnew=w
        for (j in 1:ngr){
          cond=w==ind[j]
          wnew[cond]=j
        }
        w=wnew
      }
      
      #if re-ordering does not need to be done
      ind=sample(1:nspp,size=1)
      uni.gr=sample(unique(w),size=1)
      w[ind]=uni.gr
    }
  
    #store results
    store.betas2[i,]=betas2
    store.betas1[i,]=betas1
    store.m.betas1[i,]=m.betas1
    store.tau2.betas1[i,]=tau2.betas1
    store.gammas[i,]=gammas    
    store.m.gamma[i,]=m.gamma
    store.tau2.gamma[i,]=tau2.gamma
    store.w[i,]=w
    store.theta[i,]=theta
    store.llk[i]=llk
  }
  seq1=nburn:ngibbs
  list(betas1=store.betas1[seq1,],m.betas1=store.m.betas1[seq1,],tau2.betas1=store.tau2.betas1[seq1,],
       betas2=store.betas2[seq1,],
       gammas=store.gammas[seq1,],m.gamma=store.m.gamma[seq1,],tau2.gamma=store.tau2.gamma[seq1,],
       w=store.w[seq1,],theta=store.theta[seq1,],
       llk=store.llk)
}