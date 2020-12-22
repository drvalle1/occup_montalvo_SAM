compare1=function(estim,true){
  rango=range(c(true,estim))
  plot(true,estim,ylim=rango,xlim=rango)
  lines(rango,rango,col='red',lwd=2)  
}

plot(mod1$llk,type='l')
seq1=500:ngibbs
plot(mod1$llk[seq1],type='l')

#how many groups
plot(mod1$theta[ngibbs,],type='h')

#compare group membership
aux=data.frame(w.true=w.true,w.estim=mod1$w[ngibbs,])
k=table(aux); k
seq1=c(2,1,3)
k[,seq1]

#compare betas
betas1.estim=matrix(mod1$betas1[ngibbs,],nparam.occ1,nspp)
compare1(estim=betas1.estim,true=betas1.true)

#compare betas2
ngr=10
betas2.estim=matrix(mod1$betas2[ngibbs,],nparam.occ2,ngr)
compare1(betas2.estim[,seq1],betas2.true)

#compare m.betas1
m.betas1.estim=mod1$m.betas1[ngibbs,]
compare1(estim=m.betas1.estim,true=m.betas1.true)

#compare tau2.betas1
tau2.betas1.estim=mod1$tau2.betas1[ngibbs,]
compare1(estim=tau2.betas1.estim,true=tau2.betas1.true)

#compare gammas
gammas.estim=matrix(mod1$gammas[ngibbs,],nparam.det,nspp)
compare1(estim=gammas.estim,true=gammas.true)

#compare m.gammas
m.gammas.estim=mod1$m.gamma[ngibbs,]
compare1(estim=m.gammas.estim,true=m.gamma.true)

#compare tau2.gammas
tau2.gammas.estim=mod1$tau2.gamma[ngibbs,]
compare1(estim=tau2.gammas.estim,true=tau2.gamma.true)