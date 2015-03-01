## double robust estimator two covariates
#
source("basic2.R")

sigma = 0.3 #global parameter
n = 10^5 #global parameter
sum_names <- c("intercept","X1","X2","hatQ")
p_level <- 1
eta_range <- 1:3
sq<-rep(0,3)
pop.size<-1500
it.num<-8
nvars<-length(sq)
Domains<-cbind(rep(-1,nvars),rep(1,nvars)) 
x1 <- runif(n)
x2 <- runif(n)
x<-cbind(1,x1,x2) 

tp <- ilogit( -0.5 + 0.8*x1^2+ x2^2 )
# incorrect model： tp2 <- ilogit(-1 + x1+ x2 )
a<-rbinom(n,1,tp)                                                        

y <- mean_function(x1,x2,a)


############################## Propensity score model    ########################
logit.true <- glm(a~I(x1^2)+I(x2^2),family=binomial, epsilon=1e-14)
ph.true<- as.vector(logit.true$fit)


summary_20_IPWE_PST <-   qestimate_pop(x,y,a,ph.true,.2,p_level)
#eta_25_IPWE_PST <-summary_25_IPWE_PST[eta_range]
#perf_25_IPWE_PST     <- summary_25_IPWE_PST[5:7]
save.image("popu_0.20q.RData")






