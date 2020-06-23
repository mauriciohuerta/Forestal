### Ejecutar RequiredPackages en primera instancia una unica vez de ser necesario.
### Este código genera la función plsGAMLSS.

library(rbsmodels)
library(gamlss)
library(bootBCa)

#########################################################
#########################################################
##                                                     ##
##--------------------- plsGAMLSS ---------------------##
##--------- Partial least squares for GAMLSS ----------##
##                                                     ##
#########################################################
#########################################################
##-----------------------------------------------------##
## Inputs                                              ##
##-----------------------------------------------------##
##                                                     ##
## y       response variable                           ##
## x       covariates                                  ##
## nc      number of components                        ##
## family  a gamlss.family object, which is used to    ##  
#          define the distribution and the link        ##
#          functions of the various parameters         ##
## ...     extra arguments for gamlss function         ##
##                                                     ##
##-----------------------------------------------------##
## Outputs                                             ##
##-----------------------------------------------------##
##                                                     ##
## A              unnormalized weigth matrix           ##
## W              weigth matrix                        ##
## R              corrected weigth matrix              ##
## P              a matrix of X loadings               ##
## Q              Y loading vector                     ##
## T              factor matrix (or PLS components)    ##
## FullModels     a list with the nc models            ##
## FinalModel     a latest model with nc components    ##
## betas.center   regression coefficients for          ##
##                standarized covariates               ##
## betas          regression coefficients for original ##
##                covariates                           ##
## fitted.values  fitted values of the model           ##
## scaled.X       matrix of standarized covariates     ##
## dataY          response variable                    ##
## X.residual     a matrix of residual covariates      ##
## nc             number of components                 ## 
## family         family use in gamlss. See            ##
##                family.gamlss() for details          ##
##                                                     ##
#########################################################
#########################################################

plsGAMLSS<- function(y, x, nc, family, ...){
  cat(paste("\nComputing",as.character(nc),"PLS components\n"))
  if (family$family[1]=="RBS"){
    cat(paste("\nFitted distribution: ","Reparameterized Birnbaum-Saunders","\n")) 
  }
  else{
    cat(paste("\nFitted distribution: ",family$family[2],"\n"))}
  cat(paste("Mu Link: ",family$mu.link,"\n"))
  cat("\nStatus: Computing\n\n")
  X<-as.matrix(x)
  dataY<-as.matrix(y)
  dataX<-scale(X)
  n<- nrow(dataX)
  k<-ncol(dataX)
  AA<-matrix(rep(0, k*nc), nrow=k)
  WW<-matrix(rep(0, k*nc), nrow=k)
  TT<-matrix(rep(0, n*nc),  nrow=n)
  PP<-matrix(rep(0, k*nc),  nrow=k)
  X.residual<-dataX
  models<-list()
  for (i in 1:nc){
    cat(paste("________Component",as.character(i),"________\n"))
    if (i == 1){
      for (j in 1:k){
        AA[j,1]<-gamlss(dataY~dataX[,j], family=family, trace=FALSE, ...)$mu.coeff[2]
      }
      WW[,1]<-AA[,1]/sqrt(sum(AA[,1]^2))
      TT[,1]<-dataX%*%WW[,1]
      models[[1]]<-gamlss(dataY~TT[,1], family=family, trace=FALSE, ...)
      for (j in 1:k){
        PP[j,1]<-lm(X.residual[,j]~0+TT[,1])$coeff[1]
        X.residual[,j]<-lm(X.residual[,j]~TT[,1])$resid
      }
    }
    if (i > 1){
      for (j in 1:k){
        AA[j,i]<-gamlss(dataY~dataX[,j]+TT[,1:i], family=family, trace=FALSE, ...)$mu.coeff[2]
      }
      WW[,i]<-AA[,i]/sqrt(sum(AA[,i]^2))
      TT[,i]<-X.residual%*%WW[,i]
      models[[i]]<-gamlss(dataY~TT[,1:i], family=family, trace=FALSE, ...)
      for (j in 1:k){
        PP[j,i]<-lm(X.residual[,j]~0+TT[,i])$coeff[1]
        X.residual[,j]<-lm(X.residual[,j]~TT[,i])$resid
      }
    }
  }
  FM<-gamlss(dataY~TT, family=family, trace=FALSE, ...)
  QQ<-FM$mu.coeff
  RR<-WW%*%solve(t(PP)%*%WW)
  betas.center<-RR%*%as.matrix(QQ[-1])
  intercept<-QQ[1]-sum(betas.center*apply(X, 2, mean)/apply(X, 2, sd))
  betas<-rbind(intercept,betas.center/apply(X, 2, sd))
  FV<-fitted(FM)
  cat("\nStatus: OK\n")
  output<-list(A=AA,W=WW,R=RR,T=TT,P=PP,FinalModel=FM,Q=as.matrix(QQ),betas.center=betas.center,betas=betas,FullModels=models,
               fitted.values=FV, scaled.X=dataX, dataY=dataY,X.residual=X.residual,nc=nc,family=family, dataX=X)
}

#########################################################
#########################################################
##                                                     ##
##------------------- plsGAMLSS.add -------------------##
##-------- Computes additional components for ---------##
##--------- Partial least squares for GAMLSS ----------##
##                                                     ##
#########################################################
#########################################################
##-----------------------------------------------------##
## Inputs                                              ##
##-----------------------------------------------------##
##                                                     ##
## obj     an plsGAMLSS object                         ##
## nc.add  number of additional PLS components to add  ##
## nc      number of components                        ##
## ...     extra arguments for gamlss function         ##
##                                                     ##
##-----------------------------------------------------##
## Outputs                                             ##
##-----------------------------------------------------##
##                                                     ##
## A              unnormalized weigth matrix           ##
## W              weigth matrix                        ##
## R              corrected weigth matrix              ##
## P              a matrix of X loadings               ##
## Q              Y loading vector                     ##
## T              factor matrix (or PLS components)    ##
## FullModels     a list with the nc models            ##
## FinalModel     a latest model with nc components    ##
## betas.center   regression coefficients for          ##
##                standarized covariates               ##
## betas          regression coefficients for original ##
##                covariates                           ##
## fitted.values  fitted values of the model           ##
## scaled.X       matrix of standarized covariates     ##
## dataY          response variable                    ##
## X.residual     a matrix of residual covariates      ##
## nc             number of components                 ## 
## family         family use in gamlss. See            ##
##                family.gamlss() for details          ##
##                                                     ##
#########################################################
#########################################################

plsGAMLSS.add<-function(obj, nc.add,...){
  X<-obj$dataX
  dataX<-obj$scaled.X
  dataY<-obj$dataY
  n<- nrow(dataX)
  k<-ncol(dataX)
  AA<-cbind(obj$A,matrix(rep(0, k*nc.add), nrow=k))
  WW<-cbind(obj$W,matrix(rep(0, k*nc.add), nrow=k))
  TT<-cbind(obj$T,matrix(rep(0, n*nc.add),  nrow=n))
  PP<-cbind(obj$P,matrix(rep(0, k*nc.add),  nrow=k))
  X.residual<-obj$X.residual
  models<-obj$FullModels
  nc<-obj$nc+nc.add
  family<-obj$family
  cat(paste("\nComputing",as.character(nc.add),
            "additional PLS components: from nc=",
            as.character(obj$nc),"to",as.character(nc)," \n"))
  if (family$family[1]=="RBS"){
    cat(paste("\nFitted distribution: ","Reparameterized Birnbaum-Saunders","\n")) 
  }
  else{
    cat(paste("\nFitted distribution: ",family$family[2],"\n"))}
  cat(paste("Mu Link: ",family$mu.link,"\n"))
  cat("\nStatus: Computing\n\n")
  for (i in (obj$nc+1):nc){
    cat(paste("________Component",as.character(i),"________\n"))
    for (j in 1:k){
      AA[j,i]<-gamlss(dataY~dataX[,j]+TT[,1:i], family=family, trace=FALSE, ...)$mu.coeff[2]
    }
    WW[,i]<-AA[,i]/sqrt(sum(AA[,i]^2))
    TT[,i]<-X.residual%*%WW[,i]
    models[[i]]<-gamlss(dataY~TT[,1:i], family=family, trace=FALSE, ...)
    for (j in 1:k){
      PP[j,i]<-lm(X.residual[,j]~0+TT[,i])$coeff[1]
      X.residual[,j]<-lm(X.residual[,j]~TT[,i])$resid
    }
  }
  FM<-gamlss(dataY~TT, family=family, trace=FALSE, ...)
  QQ<-FM$mu.coeff
  RR<-WW%*%solve(t(PP)%*%WW)
  betas.center<-RR%*%as.matrix(QQ[-1])
  intercept<-QQ[1]-sum(betas.center*apply(X, 2, mean)/apply(X, 2, sd))
  betas<-rbind(intercept,betas.center/apply(X, 2, sd))
  FV<-fitted(FM)
  cat("\nStatus: OK\n")
  output<-list(A=AA,W=WW,R=RR,T=TT,P=PP,FinalModel=FM,Q=as.matrix(QQ),betas.center=betas.center,betas=betas,FullModels=models,
               fitted.values=FV, dataX=X, scaled.X=dataX, dataY=dataY,X.residual=X.residual,nc=nc,family=obj$family)
}

bootX<-function(x,data){
  p<-as.numeric(lm(data[x,1]~data[x,-1])$coeff[ncol(data)])
  return(p)
}
bootY<-function(x,data,obj){
  q<-as.numeric(gamlss(data[x,1]~data[x,-1],family=obj$family,trace=FALSE)$mu.coeff[ncol(data)])
  return(q)
}

max.nc<-function(obj,B,alpha){
  n<-nrow(obj$scaled.X)
  FM<-rep(NA,ncol(obj$scaled.X))
  c<-0
  count=ncol(obj$scaled.X)
  while((count>0)&(c<(obj$nc))){
    c<-c+1
    print(c)
    i.X<-rep(0,ncol(obj$scaled.X))
    for(j in 1:ncol(obj$scaled.X)){
      boot.lu<-BCa(x=1:n,M=B,theta=bootX,alpha=c(alpha/2,(1-(alpha/2))),data=cbind(obj$scaled.X[,j],obj$T[,1:c]),delta=NA)[4:5]
      i.X[j]<-!((boot.lu[1]<0)&(0<boot.lu[2]))
    }
    FM<-cbind(FM,i.X)
    count<-sum(i.X)
  }
  if((c==obj$nc)&(sum(i.X)>0)){
    cat(paste("\nThe maximum number of components is",as.character(c)," \n"))
  }
  else{
    cat(paste("\nThe maximum number of components is",as.character(c-1)," \n"))
  }
  output<-list(c.max=c,FullMatrix=FM[,-1])
}

optim.nc<-function(cmax,alpha,B,obj){
  c=0
  lower=9999
  CI<-as.numeric()
  while((lower>0)&(c<cmax)){
    c<-c+1
    print(c)
    lower<-BCa(x=1:nrow(obj$scaled.X),M=B,theta=bootY,alpha=alpha,data=cbind(obj$dataY,obj$T[,1:c]),delta=NA,obj=obj)[4]
    CI[c]<-lower
  }
  if((c==cmax)&(lower>0)){
    cat(paste("\nThe maximum number of components is",as.character(c)," \n"))
  }
  else{
    cat(paste("\nThe maximum number of components is",as.character(c-1)," \n"))
  }
  return(CI)
}

betas.pls<-function(obj,nc){
  FM<-gamlss(obj$dataY~obj$T[,nc], family=obj$family, trace=FALSE)
  QQ<-FM$mu.coeff
  RR<-obj$W[,nc]%*%solve(t(obj$P[,nc])%*%obj$W[,nc])
  betas.center<-RR%*%as.matrix(QQ[-1])
  intercept<-QQ[1]-sum(betas.center*apply(obj$dataX, 2, mean)/apply(obj$dataX, 2, sd))
  betas<-rbind(intercept,betas.center/apply(obj$dataX, 2, sd))
  FV<-fitted(FM)
  output<-list(betas=betas,fitted=FV)
}
