sim_data_IPWE_PST = NULL
sim_data_IPWE_PSM = NULL
sim_data_DR_PST   = NULL
sim_data_DR_PSM   = NULL
for(i in 1:cores){
  sim_data_IPWE_PST <- rbind(sim_data_IPWE_PST, results[[i]]$sim_data_IPWE_PST )
  sim_data_IPWE_PSM <- rbind(sim_data_IPWE_PSM, results[[i]]$sim_data_IPWE_PSM )
  sim_data_DR_PST  <- rbind(sim_data_DR_PST , results[[i]]$sim_data_DR_PST  )
  sim_data_DR_PSM  <- rbind(sim_data_DR_PSM , results[[i]]$sim_data_DR_PSM  )
}

sim_data_IPWE_PST <- data.frame(sim_data_IPWE_PST)
colnames(sim_data_IPWE_PST)[8] <- "method"

sim_data_IPWE_PSM <- data.frame(sim_data_IPWE_PSM)
colnames(sim_data_IPWE_PSM)[8] <- "method"

sim_data_DR_PST <- data.frame(sim_data_DR_PST)
colnames(sim_data_DR_PST)[8] <- "method"

sim_data_DR_PSM <- data.frame(sim_data_DR_PSM)
colnames(sim_data_DR_PSM)[8] <- "method"

raw_analysis <- function(sim_data){
  eta_int <- with(sim_data, tapply(intercept, method, mean))
  eta_x1 <- with(sim_data, tapply(X1, method, mean))
  eta_x2 <- with(sim_data, tapply(X2, method, mean))
  perf_mean <- with(sim_data, tapply(mean, method, mean))
  #perf_q10 <- with(sim_data, tapply(q10, method, mean))
  perf_q20 <- with(sim_data, tapply(q20, method, mean))
  perf_q50 <- with(sim_data, tapply(q50, method, mean))
  
  eta_int_sd <- with(sim_data, tapply(intercept, method, sd))
  eta_x1_sd <- with(sim_data, tapply(X1, method, sd))
  eta_x2_sd <- with(sim_data, tapply(X2, method, sd))
  perf_mean_sd <- with(sim_data, tapply(mean, method, sd))
  #perf_q10_sd <- with(sim_data, tapply(q10, method, sd))
  perf_q20_sd <- with(sim_data, tapply(q20, method, sd))
  perf_q50_sd <- with(sim_data, tapply(q50, method, sd))
  performance_hatQ <- with(sim_data, tapply(hatQ, method, mean))
  performance_hatQ_sd <- with(sim_data, tapply(hatQ, method, sd))
  
  table <- data.frame(methods=c("mean","Q.20","Q.50"),
                               eta_int=paste(round(eta_int,3), "(",round(eta_int_sd,3),")",sep=""),
                               eta_x1 =paste(round(eta_x1,3),  "(",round(eta_x1_sd,3) ,")",sep=""),
                               eta_x2 =paste(round(eta_x2,3),  "(",round(eta_x2_sd,3) ,")",sep=""),
                               hatQ   =paste(round(performance_hatQ,3), "(",round(performance_hatQ_sd,3) ,")",sep=""),
                               perf_mean = paste(round(perf_mean,3),  "(",round(perf_mean_sd,3) ,")",sep=""),
                               #perf_q10 =  paste(round(perf_q10,2),  "(",round(perf_q10_sd,2) ,")",sep=""),
                               perf_q20 =  paste(round(perf_q20,3),  "(",round(perf_q20_sd,3) ,")",sep=""),
                               perf_q50 =  paste(round(perf_q50,3),  "(",round(perf_q50_sd,3) ,")",sep="") )
  return(table)
}

raw_analysis2 <- function(sim_data){
  eta_int <- with(sim_data, tapply(intercept, method, mean))
  eta_x1 <- with(sim_data, tapply(X1, method, mean))
  eta_x2 <- with(sim_data, tapply(X2, method, mean))
  perf_mean <- with(sim_data, tapply(mean, method, mean))
  #perf_q10 <- with(sim_data, tapply(q10, method, mean))
  perf_q20 <- with(sim_data, tapply(q20, method, mean))
  perf_q50 <- with(sim_data, tapply(q50, method, mean))
  
  eta_int_sd <- with(sim_data, tapply(intercept, method, sd))
  eta_x1_sd <- with(sim_data, tapply(X1, method, sd))
  eta_x2_sd <- with(sim_data, tapply(X2, method, sd))
  perf_mean_sd <- with(sim_data, tapply(mean, method, sd))
  #perf_q10_sd <- with(sim_data, tapply(q10, method, sd))
  perf_q20_sd <- with(sim_data, tapply(q20, method, sd))
  perf_q50_sd <- with(sim_data, tapply(q50, method, sd))
  performance_hatQ <- with(sim_data, tapply(hatQ, method, mean))
  performance_hatQ_sd <- with(sim_data, tapply(hatQ, method, sd))
  
  table <- data.frame(methods=c("Q.20","Q.50"),
                      eta_int=paste(round(eta_int,3), "(",round(eta_int_sd,3),")",sep=""),
                      eta_x1 =paste(round(eta_x1,3),  "(",round(eta_x1_sd,3) ,")",sep=""),
                      eta_x2 =paste(round(eta_x2,3),  "(",round(eta_x2_sd,3) ,")",sep=""),
                      hatQ   =paste(round(performance_hatQ,3), "(",round(performance_hatQ_sd,3) ,")",sep=""),
                      perf_mean = paste(round(perf_mean,3),  "(",round(perf_mean_sd,3) ,")",sep=""),
                      #perf_q10 =  paste(round(perf_q10,2),  "(",round(perf_q10_sd,2) ,")",sep=""),
                      perf_q20 =  paste(round(perf_q20,3),  "(",round(perf_q20_sd,3) ,")",sep=""),
                      perf_q50 =  paste(round(perf_q50,3),  "(",round(perf_q50_sd,3) ,")",sep=""))
  return(table)
}



table_IPWE_PST <-raw_analysis(sim_data_IPWE_PST)
table_IPWE_PSM <-raw_analysis(sim_data_IPWE_PSM)
table_DR_PST <-raw_analysis2(sim_data_DR_PST)
table_DR_PSM <-raw_analysis2(sim_data_DR_PSM)
table_IPWE_PST 
table_IPWE_PSM 
table_DR_PST 
table_DR_PSM


Table[ order(Table$methods), ]
library(xtable)

save(table_IPWE_PST,table_IPWE_PSM,
     table_DR_PST,  table_DR_PSM,file="1000.RData")

