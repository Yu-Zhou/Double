## double robust estimator two covariates
#
source("basic2.R")

sigma = 0.3 #global parameter
n = 500 #global parameter
sum_names <- c("intercept","X1","X2","hatQ")
p_level <- 1
eta_range <- 1:3
sq<-rep(0,3)
pop.size<-600
it.num<-7
nvars<-length(sq)
Domains<-cbind(rep(-1,nvars),rep(1,nvars)) 

sim_function <- function(sim_num=2,p_level=0,hard_limit=FALSE){
  sim_data_IPWE_PST<- sim_data_IPWE_PSM <- sim_data_DR_PST <- sim_data_DR_PSM <- NULL
  #######################  data generation #######################
  for(i in 1:sim_num){
    ptm0 <- proc.time()
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
    
    logit.m <- glm(a~ x1 ,family=binomial, epsilon=1e-14)
    ph.m<- as.vector(logit.m$fit)
    
    ###################################### 1. estimate eta for IPWE ， true PS   ########
    mean_summary_IPWE_PST <- mestimate(x,y,a,ph.true,p_level)
    #mean_eta_IPWE_PST   <- mean_summary_IPWE_PST[eta_range]
    #perf_mean_IPWE_PST  <- mean_summary_IPWE_PST[5:7]
    
    #summary_10_IPWE_PST <-   qestimate(x,y,a,ph.true,.1,p_level)
    #eta_10_IPWE_PST     <- summary_10_IPWE_PST[eta_range]    
    #perf_10_IPWE_PST     <- summary_10_IPWE_PST[5:7]
    
    summary_20_IPWE_PST <-   qestimate(x,y,a,ph.true,.2,p_level)
    #eta_25_IPWE_PST <-summary_25_IPWE_PST[eta_range]
    #perf_25_IPWE_PST     <- summary_25_IPWE_PST[5:7]
    summary_50_IPWE_PST <-   qestimate(x,y,a,ph.true,.5,p_level)
    
    ###################################### 2. estimate eta for IPWE ， misspecified PS   ########
    mean_summary_IPWE_PSM <- mestimate(x,y,a,ph.m,p_level)
    #mean_eta_IPWE_PSM    <- mean_summary_IPWE_PSM[eta_range]
    #perf_mean_IPWE_PSM   <- mean_summary_IPWE_PSM[5:7]
    
    #summary_10_IPWE_PSM <-   qestimate(x,y,a,ph.m,.1,p_level)
    #eta_10_IPWE_PSM     <- summary_10_IPWE_PST[eta_range]    
    #perf_10_IPWE_PSM    <- summary_10_IPWE_PST[5:7]
    
    summary_20_IPWE_PSM <-   qestimate(x,y,a,ph.m,.20,p_level)
    #eta_25_IPWE_PSM <-summary_25_IPWE_PST[eta_range]
    #perf_25_IPWE_PSM    <- summary_25_IPWE_PST[5:7]
    summary_50_IPWE_PSM <-   qestimate(x,y,a,ph.m,.5,p_level)
    
    ##################################### 3. DR estimate optimal eta for quantile criterion, true PS ################
    #summary_10_DR_PST <-   r_qestimate(x,y,a,ph.true, .1,p_level)
    #eta_10_DR_PST     <- summary_10_DR_PST[eta_range] 
    #perf_10_DR_PST    <- summary_10_DR_PST[5:7]
    
    summary_20_DR_PST <-   r_qestimate(x,y,a,ph.true,.2,p_level)
    #eta_25_DR_PST     <- summary_25_DR_PST[eta_range] 
    #perf_25_DR_PST    <- summary_25_DR_PST[5:7]
    summary_50_DR_PST <-   r_qestimate(x,y,a,ph.true, .5,p_level)
    
    ##################################### 4. DR estimate optimal eta for quantile criterion, misspecified PS ################
    #summary_10_DR_PSM <-   r_qestimate(x,y,a,ph.m,.1,p_level)
    #eta_10_DR_PSM     <- summary_10_DR_PSM[eta_range] 
    #perf_10_DR_PSM    <- summary_10_DR_PSM[5:7]
    
    summary_20_DR_PSM <-   r_qestimate(x,y,a,ph.m,.2,p_level)
    #eta_25_DR_PSM     <- summary_25_DR_PSM[eta_range] 
    #perf_25_DR_PSM    <- summary_25_DR_PSM[5:7]
    
    summary_50_DR_PSM <-   r_qestimate(x,y,a,ph.m,.50,p_level)
    
    sim_data_IPWE_PST <- rbind(sim_data_IPWE_PST, 
                               c(mean_summary_IPWE_PST ,   1),
                  #             c(summary_10_IPWE_PST   ,  10), 
                               c(summary_20_IPWE_PST   ,  20),
                               c(summary_50_IPWE_PST   ,  50)
    )
    
    sim_data_IPWE_PSM <- rbind(sim_data_IPWE_PSM,
                               c(mean_summary_IPWE_PSM  ,1),
                               #c(summary_10_IPWE_PSM   ,10),
                               c(summary_20_IPWE_PSM   ,20),
                               c(summary_50_IPWE_PSM   ,50)
    )
    
    sim_data_DR_PST <- rbind(sim_data_DR_PST,
                             #c(summary_10_DR_PST ,10),
                             c(summary_20_DR_PST ,20),
                             c(summary_50_DR_PST ,50)
    )
    
    sim_data_DR_PSM <- rbind(sim_data_DR_PSM,
                             #c(summary_10_DR_PSM ,10),
                             c(summary_20_DR_PSM ,20),
                             c(summary_50_DR_PSM ,50)
    )
    ptm1=proc.time() - ptm0
    jnk=as.numeric(ptm1[3])
    cat('\n','It took ', jnk, "seconds, Iteration:", i,'\n')
    print(paste0("Current working dir: ", i))
  }
  #comments start here  
  list(sim_data_IPWE_PST = sim_data_IPWE_PST ,
       sim_data_IPWE_PSM = sim_data_IPWE_PSM ,
       sim_data_DR_PST   = sim_data_DR_PST   ,
       sim_data_DR_PSM   = sim_data_DR_PSM)
}
# sim_function(1)
cores <- 8
p_level <- 0
results <- mclapply(rep(15,cores), sim_function, p_level,hard_limit=TRUE, mc.cores = cores)

save.image("Sigma_0.3_500_ver2.RData")
