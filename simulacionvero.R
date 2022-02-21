library(spatstat)
library(cubature)
N=500
sig=0.1
W=square(1)
Warea=area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W)

Z<-attr(X,"Lambda")
#covarianza=varianza^2*exp(-d*escala)
plot(Z)
points(X,pch=16)
X
#usar runif
rma=0.2
#Funcion close pairs
cove<-function(d,par){
  return(par[1]^2*exp(-d/exp(par[2])))
}
integrand<-function(arg,par){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par))*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}
wcuad<-function(d,rma){
  return((d+rma)*(d-rma)*-1/rma)
}
P<-closepairs(X,0.15,twice = FALSE)
a=sum(log(cove(d=P$d,par=c(0.1,0.15))))
print(wcuad(d=dU$d,rma=rma))
dU$d

lgv<- function(par,rmax=0.2,cova,w0f,X){
  P<-closepairs(X,rmax,twice = FALSE)
  #w=w0f(d=P$d,rma=rmax)
  #Sumatoria de los primeros terminos
  a=cova(d=P$d,par=par)
  #si no entra la simulacion Montecarlo
  #x0<-runifpoint(mc,win=square(1))
  #y0<-runifpoint(mc,win=square(1))
  #dU<-crosspairs(x0,y0,rmax,what="ijd")
  int=vegas(integrand, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par)
  ver=(a-log(int$int))
  ver=-sum(ver)
  return(ver)
}

rv<-optim(par=c(sqrt(0.7),log(0.05)),X=Y,rmax=rma,cova=cove,w0f=wcuad,fn=lgv)
rv2<-c((rv$par[1])^2,exp(rv$par[2]))
rv2
Y=X
npoints=c()

Y=X
Y
rsig=c()
rphi=c()
print(lgv(X=Y,rmax=rma,mc=nmc,cova=cove,w0f=wcuad,dU=dU$d,par=c(0.1,0.15)))
nmc=1000000
x0<-runifpoint(nmc,win=square(1))
y0<-runifpoint(nmc,win=square(1))
dU<-crosspairs(x0,y0,rma,what="ijd")
for (i in 1:20){
x0<-runifpoint(nmc,win=square(1))
y0<-runifpoint(nmc,win=square(1))
dU<-crosspairs(x0,y0,rma,what="ijd")
rv<-optim(par=c(sqrt(0.5),log(0.5)),X=Y,rmax=rma,mc=nmc,cova=cove,w0f=wcuad,dU=dU$d,fn=lgv)
rsig=c(rsig,rv$par[1])
rphi=c(rphi,rv$par[2])
} 
si20<-data.frame(sig=rsig,phi=rphi)
si20$var=si20$sig^2
si20
Y
npoints=c(npoints,883)
ressim=data.frame(sig,phi)

si1
int=vegas(integrand, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
      relTol=1e-3, absTol=1e-12,
      flags=list(verbose=0, final=1),par=c(1,0.15))
integrand<-function(arg,par){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par))*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}
integrand(arg=c(0.05,0.002,0.003,0.004),par=c(1,0.15))
