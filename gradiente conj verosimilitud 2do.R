library(spatstat)
library(cubature)
#simulacion LGCP
N=500
sig=0.1
W=square(1)
Warea=area(W)
mu0=log(N/Warea)-sig/2
phi=0.15
X <- rLGCP(model="exp", mu=mu0, var=sig, scale=phi, win = W)

Z<-attr(X,"Lambda")
plot(Z)
points(X,pch=16)
X

#funciones necesarias


#Funcion de covarianza
cove<-function(d,par,par2){
  a=0
  for (i in 1:length(par)){
    a=a+(par[i])^2*exp(-d/par2[i]) 
    # falta la funcion k
  }
  return(a)
}
#integral1 
integrand1<-function(arg,par,par2){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}

#integral2
integrand2<-function(arg,par,par2,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*2*par[j]*exp(-d/par2[j])*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}

#gradiente
grad<-function(par,par2,j){
  int=vegas(integrand1, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par,par2=par2)
  int2=vegas(integrand2, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
             relTol=1e-3, absTol=1e-6,
             flags=list(verbose=0, final=1),par=par,par2=par2,j=j)
  P<-closepairs(X,0.2,twice = FALSE)
  a=2*par[j]*exp(-P$d/par2[j])
  ver=a-int$int/int2$int
  ver=sum(ver)
  return(ver)
}

#varianzas (parametros a estimar)
par=c(sqrt(0.5),sqrt(0.5),sqrt(0.5),sqrt(0.5),sqrt(0.5))

#phi
par2=c(0.05,0.1,0.15,0.2,0.25)
alp=0.00001
lam=0
for (k in 1:200){
  para=par
  print((par)^2)
  for (i in 1:(length(par))){
      x0=para[i]+alp*grad(par,par2,i)
      par[i]=sign(x0)*pmax(abs(x0)-lam*alp,0)}
  }

