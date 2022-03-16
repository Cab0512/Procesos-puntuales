library(spatstat)
library(cubature)
inten<- function(x,y) (25*x^2+16*x+4*y+9*y^2)
X <-rpoispp(100,lambda=inten)
plot(X)
integrand<-function(arg,j){
  x<-arg[1]
  y<-arg[2]
  d=c(0,0,0,0,0)
  d[j]=1
  w=d[1]*x^2+d[2]*x^1+d[3]+d[4]*y+d[5]*y^2
  return(w)
}
grad<-function(X,b,j){
  int=vegas(integrand, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
            relTol=1e-3, absTol=1e-12,
            flags=list(verbose=0, final=1) ,j=j)
return(2*X$n/b[j]-2*b[j]*int$int)
}
int=vegas(integrand, lowerLimit = rep(0, 2), upperLimit = rep(1, 2),
          relTol=1e-3, absTol=1e-12,
          flags=list(verbose=0, final=1) ,j=1)
int$int
par=c(sqrt(1),sqrt(1),sqrt(1),sqrt(1),sqrt(1))
lam=0
alp=0.001
grad(X,par,5)
for (k in 1:100){
  para=par
  print(par^2)
  for (i in 1:(length(par))){
    x0=para[i]-alp*grad(X,par,i)
    par[i]=sign(x0)*pmax(abs(x0)-lam*alp,0)
  }
}
X$n

