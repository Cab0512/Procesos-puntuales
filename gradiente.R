cove<-function(d,par,par2){
  a=0
  for (i in 1:length(par)){
  a=a+par[i]*cos(par2[i]*d) 
  # falta la funcion k
  }
  return(a)
}

integrand1<-function(arg,par,par2,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}

integrand2<-function(arg,par,par2,j){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par,par2))*cos(par2[j]*d)*(1*(d<=0.2) + 0* (d>=0.2))
  return(w)
}

grad<-function(par,par2,j,rmax){
  int=vegas(integrand1, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par,par2=par2,j=j)
  int2=vegas(integrand2, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
             relTol=1e-3, absTol=1e-6,
             flags=list(verbose=0, final=1),par=par,par2=par2,j=j)
  P<-closepairs(X,rmax,twice = FALSE)
  a=cos(par2[j]*P$d)
  ver=a-int$int/int2$int
  ver=-sum(ver)
  return(ver)
}
par=c(1,2,3,4,5,6,7)
par2=c(2,3,4,5,6,7)
grad(par,2,0.2)

#para cualquier f

inte<-function(arg,par,cova,fgrad){
  x1<-arg[1]
  x2<-arg[2]
  x3<-arg[3]
  x4<-arg[4]
  d=sqrt(((x1-x3)^2)+((x2-x4)^2))
  w=exp(cove(d,par))*fgrad(d,par)
  return(w)
}

#no se agrega el patron puntual para ahorrar tiempo
est<-function(par,f,cova){
  int=vegas(inte, lowerLimit = rep(0, 4), upperLimit = rep(1, 4),
            relTol=1e-3, absTol=1e-6,
            flags=list(verbose=0, final=1),par=par,cova=cova,f=f)
  a=f(P$d,par)
  a=sum(a)-int
  return(a)        
}

