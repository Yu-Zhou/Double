library(rgenoud)
library(quantreg)
library(xtable)
library(faraway)
library(Matrix)
library(lpSolve)
Mean_Function <- function(x1,x2,a)
{
  y = exp( 0.3*(1+x1+x2)+0.3*a*(x1+x2-1)) + sigma*(2+0.8*a*((x1+x2)>1) )*rnorm(n,0,1)
  y
}

Mean_Est<-function(eta,x,y,a,prob)
{
  #mean estimator
  g<-as.numeric(I(x%*%eta>0))
  c<-a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  val <- mean(c*y*wts)
  #print(c(eta,val))
  val
}

R_Mean_Est <- function(eta,x,y,a,prob)
{   # optimization objective function for robust method
  g<-as.numeric(I(x%*%eta>0))
  c<-a*g + (1-a)*(1-g)
  wts <- g*1/prob + (1-g)*(1/(1-prob))
  val <- mean(c*y*wts)+ mean((1-c*wts)*(Mreg(x,y,a,eta)))
  val
}

Mreg <- function(x,y,a,eta){
  g<-as.numeric(I(x%*%eta>0))
  x1 <-x[,2]
  x2 <-x[,3]
  reg.lm <- lm(y~x1*x2 + I(x1^2) + I(x2^2) +a*(x1+x2) )
  newdata <- data.frame(x1,x2,a=g)
  predict(reg.lm,newdata)
}

Mreg2 <- function(x,y,a,eta){
  g<-as.numeric(I(x%*%eta>0))
  x1 <-x[,2]
  x2 <-x[,3]
  
  reg.lm <- lm(y~a*(x1*x2) )
  newdata <- data.frame(x1,x2,a=g)
  predict(reg.lm,newdata)
}

Quant_Est <- function(eta,x,y,a,prob,tau){
  #quantile estimator
  g <- as.numeric( I(x %*% eta > 0))
  c <- a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  wts <- c*wts
  model <- rq(y ~ 1, weights=wts, tau=tau)
  coefficients(model)[1]
}

### for R_Quant_Est###
Betahat<- function(tau, sol){ # The solution b(tau_i) 
  #prevails from tau_i to tau_i+1
  Taus <-sol[1,]
  ind <- Taus<tau
  which.max(Taus[ind])
}
Check<- function(tau,u){  u*(tau-as.numeric(u<0))  }

Test <- function(arg1,tau,qrwts.y,qrwts.a,y,y.a){
  return(sum(qrwts.y * Check(tau,(y-arg1)))
         +sum(qrwts.a * Check(tau,(y.a - arg1)) ))
}

R_Quant_Est <- function(eta,x,y,a,prob,tau, y.a.0, y.a.1){
  #quantile estimator
  g <- as.numeric( I(x %*% eta > 0))
  c <- a*g+(1-a)*(1-g)
  wts <- g*1/prob+(1-g)*(1/(1-prob))
  qrwts.y <- c*wts
  qrwts.a <- 1- qrwts.y
  
  y.a <- rep(0,n)#####  MARK
  for (i in 1:n){
     y.a[i] <- g[i]*y.a.1[i] + (1-g[i])*y.a.0[i,]
  }
    
  xleft = min(y); #need modify
  xright = max(y);
  phi = 0.34
  while (xright - xleft > 10^(-5)){
    x1=phi*xright+(1-phi)*xleft
    x2=phi*xleft+(1-phi)*xright
    if(Test(x1,tau,qrwts.y,qrwts.a,y,y.a) >  Test(x2,tau,qrwts.y,qrwts.a,y,y.a) )  
      xleft = x1
    else
      xright = x2
  }
  (xleft + xright) / 2
   #optimize(test,tau=tau,qrwts.y=qrwts.y,qrwts.a=qrwts.a,y=y,y.a=y.a,lower=0,upper=4,)
}

Qestimate<-function(x,y,a,prob,tau,p_level,hard_limit=FALSE)
{
  data <- data.frame(Y=y,X1=x[,2],X2=x[,3])
  data0 <- data[a==0,]
  data1 <- data[a==1,]
  rqp0 <- rq(Y~1+X1+X2,data=data0,tau=2)$sol  # quantreg process
  rqp1 <- rq(Y~1+X1+X2,data=data1,tau=2)$sol  # quantreg process
  
  nvars<-length(sq) 
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  
  est<-genoud(fn=Quant_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,tau=tau,print.level=p_level,max=TRUE,pop.size=pop.size,wait.generations=it.num,
              gradient.check=FALSE, BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=sq,
              hard.generation.limit=hard_limit,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  perf_est_mean <- Mean_Est(eta,x,y,a,prob)
  
  perf_est_20   <- Quant_Est(eta,x,y,a,prob, .2)
  perf_est_50   <- Quant_Est(eta,x,y,a,prob, .50)
  perf  <- c(perf_est_mean,   perf_est_20, perf_est_50)
  
  summary<-c(eta,hatQ,perf)
  names(summary)<-c(sum_names,"mean","q20","q50")
  
  return(summary)
}

R_Qestimate<-function(x,y,a,prob,tau,p_level,hard_limit=FALSE)
{
  #
  data <- data.frame(Y=y,X1=x[,2],X2=x[,3])
  data0 <- data[a==0,]
  data1 <- data[a==1,]
  rqp0 <- rq(Y~1+X1+X2,data=data0,tau=2)$sol  # quantreg process
  rqp1 <- rq(Y~1+X1+X2,data=data1,tau=2)$sol  # quantreg process
  u0 <- runif(n);u1 <- runif(n)
  betahat_index_0 <- unlist(lapply(u0,Betahat,sol= process_0))
  beta0 <- process_0[4:(3+dim(x)[2]) , betahat_index_0 ]
  betahat_index_1 <- unlist(lapply(u1,Betahat,sol = process_1))
  beta1 <- process_1[4:(3+dim(x)[2]) , betahat_index_1 ]
  y.a.1 <- rep(0,n);   y.a.0 <- rep(0,n)#####  MARK
  for (i in 1:n){
    y.a.1[i] <- x[i,]%*%beta1[,i]
    y.a.0[i] <- x[i,]%*%beta0[,i]
  }
  #
  nvars<-length(sq) 
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  
  est<-genoud(fn=R_Quant_Est,nvars=nvars,
              x=x,y=y,a=a,prob=prob,tau=tau,y.a.0=y.a.0,y.a.1=y.a.1,
              max=TRUE,print.level=p_level,pop.size=pop.size,
              wait.generations=it.num,gradient.check=FALSE, BFGS=FALSE,
              P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
              Domains=Domains,starting.values=sq,hard.generation.limit=hard_limit,
              solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  
  #estimate mean , .2 quantile and .5 quantile
  perf_est_mean <- Mean_Est(eta,x,y,a,prob)
  perf_est_20   <- R_Quant_Est(eta,x,y,a,prob, .2, y.a.0 , y.a.1)
  perf_est_50   <- R_Quant_Est(eta,x,y,a,prob, .5, y.a.0 , y.a.1)
  perf  <- c(perf_est_mean,  perf_est_20, perf_est_50)
  
  summary<-c(eta,hatQ,perf)
  names(summary)<-c(sum_names,"mean","q20","q50")
  
  return(summary)
}

Mestimate<-function(x,y,a,prob,p_level,hard_limit=FALSE)
{
  data <- data.frame(Y=y,X1=x[,2],X2=x[,3])
  data0 <- data[a==0,] 
  data1 <- data[a==1,]
  #rqp0 <- rq(Y~1+X1+X2,data=data0,tau=2,)$sol  # quantreg process
  #rqp1 <- rq(Y~1+X1+X2,data=data1,tau=2)$sol  # quantreg process
  
  nvars<-length(sq)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  
  est<-genoud(fn=Mean_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,print.level=p_level,max=TRUE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=sq,hard.generation.limit=hard_limit
              ,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  #########  for evaluating, calculate robust quantile
  rqp0 <- rq(Y~1+X1+X2,data=data0,tau=2)$sol  # quantreg process
  rqp1 <- rq(Y~1+X1+X2,data=data1,tau=2)$sol  # quantreg process
  u0 <- runif(n);u1 <- runif(n)
  betahat_index_0 <- unlist(lapply(u0,Betahat,sol= process_0))
  beta0 <- process_0[4:(3+dim(x)[2]) , betahat_index_0 ]
  betahat_index_1 <- unlist(lapply(u1,Betahat,sol = process_1))
  beta1 <- process_1[4:(3+dim(x)[2]) , betahat_index_1 ]
  y.a.1 <- rep(0,n);   y.a.0 <- rep(0,n)#####  MARK
  for (i in 1:n){
    y.a.1[i] <- x[i,]%*%beta1[,i]
    y.a.0[i] <- x[i,]%*%beta0[,i]
  }
  
  #########  estimated  eta ####################
  eta<-est$par
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  
  perf_est_mean <- Mean_Est(eta,x,y,a,prob)
  #R_Quant_Est(eta,x,y,a,prob,tau, y.a.0, y.a.1)  syntax
  perf_est_20   <- R_Quant_Est(eta,x,y,a,prob, .2, y.a.0 , y.a.1)
  perf_est_50   <- R_Quant_Est(eta,x,y,a,prob, .5, y.a.0 , y.a.1)
  perf  <- c(perf_est_mean,  perf_est_20, perf_est_50)
  
  summary<-c(eta,hatQ,perf)
  names(summary)<-c(sum_names,"mean","q20","q50")
  
  return(summary)
}

R_Mestimate <- function(x,y,a,prob,p_level,hard_limit=FALSE)   # robust estimator
{ 
  data <- data.frame(Y=y,X1=x[,2],X2=x[,3])
  data0 <- data[a==0,]
  data1 <- data[a==1,]
  
  nvars<-length(sq)
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  
  est<-genoud(fn=R_Mean_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,print.level=p_level,max=TRUE,pop.size=pop.size,wait.generations=it.num,gradient.check=FALSE,
              BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=sq,hard.generation.limit=hard_limit
              ,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
  #########  for evaluating, calculate robust quantile
  rqp0 <- rq(Y~1+X1+X2,data=data0,tau=2)$sol  # quantreg process
  rqp1 <- rq(Y~1+X1+X2,data=data1,tau=2)$sol  # quantreg process
  u0 <- runif(n);u1 <- runif(n)
  betahat_index_0 <- unlist(lapply(u0,Betahat,sol= process_0))
  beta0 <- process_0[4:(3+dim(x)[2]) , betahat_index_0 ]
  betahat_index_1 <- unlist(lapply(u1,Betahat,sol = process_1))
  beta1 <- process_1[4:(3+dim(x)[2]) , betahat_index_1 ]
  y.a.1 <- rep(0,n);   y.a.0 <- rep(0,n)#####  MARK
  for (i in 1:n){
    y.a.1[i] <- x[i,]%*%beta1[,i]
    y.a.0[i] <- x[i,]%*%beta0[,i]
  }
  
  #########  estimated  eta ####################
  eta<-est$par
  eta<-eta/sqrt(sum(eta^2)) 
  hatQ<-est$value
  
  perf_est_mean <- Mean_est(eta,x,y,a,prob)
  perf_est_20   <- R_Quant_Est(eta,x,y,a,prob, .2, y.a.0 , y.a.1)
  perf_est_50   <- R_Quant_Est(eta,x,y,a,prob, .5, y.a.0 , y.a.1)
  perf  <- c(perf_est_mean,    perf_est_20, perf_est_50)
  
  summary<-c(eta,hatQ,perf)
  names(summary)<-c(sum_names,"mean","q20","q50")
  
  return(summary)
}

Qestimate_Pop<-function(x,y,a,prob,tau,p_level,hard_limit=FALSE)
{
  data <- data.frame(Y=y,X1=x[,2],X2=x[,3])
  data0 <- data[a==0,]
  sdata0 <- data0[sample(nrow(data0),10^4),]
  data1 <- data[a==1,]
  sdata1 <- data1[sample(nrow(data1),10^4),]
  rqp0 <- rq(Y~1+X1+X2,data = sdata0,tau=2)$sol  # quantreg process
  rqp1 <- rq(Y~1+X1+X2,data = sdata1,tau=2)$sol  # quantreg process
  
  nvars<-length(sq) 
  Domains<-cbind(rep(-1,nvars),rep(1,nvars))
  
  est<-genoud(fn=Quant_Est,nvars=nvars,x=x,y=y,a=a,prob=prob,tau=tau,print.level=p_level,max=TRUE,pop.size=pop.size,wait.generations=it.num,
              gradient.check=FALSE, BFGS=FALSE,P1=50, P2=50, P3=10, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,Domains=Domains,starting.values=sq,
              hard.generation.limit=hard_limit,solution.tolerance=0.0001,optim.method="Nelder-Mead")
  
  #########  estimated  eta ####################
  eta<-est$par
  eta<-eta/sqrt(sum(eta^2))
  hatQ<-est$value
  #eta,x,y,a,prob,tau
  perf_est_mean <- Mean_est(eta,x,y,a,prob)
  perf_est_20   <- Quant_Est(eta,x,y,a,prob, .2)
  perf_est_50   <- Quant_Est(eta,x,y,a,prob, .5)
  perf  <- c(perf_est_mean,  perf_est_10,  perf_est_20, perf_est_50)
  
  summary<-c(eta,hatQ,perf)
  names(summary)<-c(sum_names,"mean","q20","q50")
  
  return(summary)
}
