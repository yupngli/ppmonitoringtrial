###################
## R配置
###################

# Benchmark
## library("microbenchmark")

# compiler
library(compiler)

# data handling
library(dplyr)

# data output
library(openxlsx)

###################
## 载入环境变量
###################
setwd("C:/Users/yli01/OneDrive - lianbio/working/Bayes2StagePP/MthdByAstrazeneca/Git_PP_monitor/ppmonitoringtrial")
source("environment_values.R")

###################
## prior先验设定
###################

# 先验设定的参考文献: http://dx.doi.org/10.21037/tcr.2019.05.17 

CalPrior <- function(prop,sd){
  a <- ((1-prop)/sd^2-1/prop)*prop^2
  b <- a*(1/prop-1)
  return(c(a,b))
}

##################################################
## Predictive Prob 贝叶斯预测概率使用的函数
##################################################

PredPower <- function(N, N_I, ns, ni) {
  #	N	:	整个试验的计划样本量
  #	N_I	:	当前中期分析（包括PP）时计划的样本量
  #	ns	:	整个试验需要的responders数量
  #	ni	:	中期分析观测到的responders数量
  
  #ISS:	中期分析时responders的数量向量,0 to N_I
  ISS <- 0:N_I		
  #mid:	中期分析时response rate向量
  mid <- ISS / N_I
  mid_e <- length(ISS)
  # Exception when ISS=0
  mid[1] <- ifelse(mid[1] == 0, 0.5 / N_I, ISS / N_I)
  # Exception when ISS=N_I
  mid[mid_e] <- (ISS[mid_e] - 0.5) / N_I
  
  # Binomial Likelihood: (N_I,ISS)*p^(ISS)*(1-p)^(N_I-ISS)
  LIK <- dbinom(ISS, N_I, ni / N_I, log = FALSE)
  
  # Prob of Success,也就是PP
  POS <-
    1 - pbinom((ns - ni) - 1,
               N - N_I ,
               mid,
               lower.tail = TRUE,
               log.p = FALSE)
  
  #	probability of success at final given interim over all interim results
  assu <- LIK * POS
  
  return(sum(assu))
}

#######################
## Compile R function
#######################

cmp_CalPrior <- cmpfun(CalPrior)
cmp_PredPower <- cmpfun(PredPower)

##################################################
## Simulation 试验模拟主体
##################################################

#================================ 试验设计 ================================
## phase A -> phase B -> phase C -> phase D
## phase A: Na个患者
## phase B: 从Na+1开始，执行PP监测，到第Na+Nb人
## phase C: 前面a,b两阶段已有Na+Nb人，若最后一个人PP没有问题，入第Na+Nb+1人做Interim
## phase D: Interim的决策决定了是否入剩余的所有人,即总样本Na+Nb+1+Nd
#================================ 试验设计 ================================

#================================ 输入参数 ================================
## Nt: 整个试验计划的总样本量
## N_i: 中期分析计划的样本量
## priora: Beta先验参数a
## priorb: Beta先验参数b
## tv: Target value
## lrv: Lower reference value
## dc_lrv: -
## ar_tv: -
## Tar: 在Phase c时，最终达成非停决定的response rate
## strt: 从该观测后开始监测PP
## stop: 从该观测后停止监测PP
## step: 每几个看一次PP？
## PPmon: 继续前进的最低PP要求
## repsim: 模拟次数
#================================ 输入参数 ================================

#================================ 决策 ================================
# 双决策(efficacy,futility)
#..... Go: Pr[p>=lrv|x] >= dc_lrv
#..... Stop: Pr[p>=tv|x] <= ar_tv

# 单决策(futility stop)
#..... Stop: Pr[p>=tv|x] <= ar_tv,等价于Pr[p<tv|x] > 1-ar_tv
#............此处可令(1-ar_tv)为futility stopping的概率边界
#================================ 决策 ================================

SimTrial <- function(Nt,
                   N_i,
                   priora,
                   priorb,
                   true_1,
                   tv,
                   lrv,
                   dc_lrv,
                   ar_tv,
                   Tar,
                   strt,
                   stop,
                   step,
                   PPmon,
                   repsim,
                   seed) { 
  
  set.seed(seed)
  
  # 创建模拟矩阵
  source_matrix <- matrix(data=NA,nrow=repsim,ncol=Nt)
  
  # 创建决策矩阵:记录决策概率
  ## 构建指示向量，标记出在哪个人的时候开始做监测（PP以及IA,FA）
  monitor_index_1 <- seq(from=strt,to=stop,by=step)
  if (max(monitor_index_1)<stop){
    monitor_index <- c(monitor_index_1,stop,N_i,Nt)
  } else {
    monitor_index <- c(monitor_index_1,N_i,Nt)
  }
  
  ## 计算监测过程中，需要看几次（每次PP其实都是一次期中分析）
  numeber_of_pp <- length(monitor_index)
  ## 构建决策概率矩阵：行数为模拟次数，列数为 PP次数+2 (一次是IA，一次是FA)
  decision_matrix <- matrix(data=NA,nrow=repsim,ncol=numeber_of_pp)
  # 创建决策矩阵:记录决策,1为succ,0为fail
  decision_matrix_decis <- matrix(data=NA,nrow=repsim,ncol=numeber_of_pp)
  
  #	Response rate的先验分布
  Bpar1 <- c(priora,priorb)
  
  # 创建模拟数据集
  ## 创建向量记录beta分布产生的response rate
  response_rate <- vector(mode = "numeric",length = repsim)
  for (row in 1:repsim){
    theta_1 <- rbeta(1, Bpar1[1] + 1, Bpar1[2] + 1)
    response_rate[row] <- theta_1
    # 为全人群Nt产生模拟的response
    all.response <- rbinom(Nt, 1, theta_1)
    source_matrix[row,] <- all.response
  }
  
  #########################################################################
  # phase A & phase B
  #########################################################################
  
  # 中期分析时,记录截至第一次PP失败时，有多少样本量
  ## 初始化向量
  pp.pass.N <- vector(mode = "numeric",length = repsim)
  
  for (row in 1:repsim){
    # 截取模拟数据集中的前N_i个数据作为Phase C中期分析的依据
    all.response.IA <- source_matrix[row,1:N_i]
    # 数据容器 - 存储PP监测的结果，因为包含了IA和FA，这里要减2
    pp.monitor <- decision_matrix[row,1:(numeber_of_pp-2)]
    # 数据容器 - 存储PP监测时，患者数量
    pp.monitor.N <- monitor_index[1:(numeber_of_pp-2)]
    # 数据容器 - 存储PP监测时，步数
    for (g in 1:(length(pp.monitor))) {
      pp.monitor[g] <-
        cmp_PredPower(Nt, pp.monitor.N[g], round(Nt * Tar, 0), max(sum(all.response.IA[1:pp.monitor.N[g]])))
    }
    # PP计算结果填入决策概率矩阵
    decision_matrix[row,1:(length(pp.monitor))] <- pp.monitor
    # 将决策结果填入决策矩阵,TRUE为成功前进,FALSE为失败停止
    pp.monitor.result <- pp.monitor>=PPmon
    decision_matrix_decis[row,1:(length(pp.monitor))] <- pp.monitor.result
    # 如果当前情况第16人成功，第20人失败，则此时实际参与监测的样本量为20
    if (min(pp.monitor.result)==TRUE){
      pp.pass.N[row] <- pp.monitor.N[numeber_of_pp-2]
    } else if (max(pp.monitor.result)==FALSE){
      pp.pass.N[row] <- pp.monitor.N[1]
    } else {
      pp.pass.N[row] <- pp.monitor.N[min(which(pp.monitor.result==FALSE))]
    }
  }
  
  #########################################################################
  
  actual_samplesize_vec <- vector(mode = "numeric",length = repsim)
  
  for (row in 1:repsim){
  
  # calculate the REAL number of corresponders among all patients
  all.response.sum <- sum(source_matrix[row,])
  
  # calculate the REAL number of NON-corresponders among all patients
  all.nonresponse.sum <- Nt - all.response.sum
  
  # calculate the REAL number of corresponders among patients in Interim Analysis
  all.response.ia.sum <- sum(source_matrix[row,1:N_i])
  
  # calculate the REAL number of NON-corresponders among patients in Interim Analysis
  all.nonresponse.ia.sum <- N_i - all.response.ia.sum
  
  # Final analysis
  ### 1:GO, 2:STOP
  all.response_go <- pbeta(lrv, 1 + all.response.sum, 1 + all.nonresponse.sum, lower.tail = FALSE)
  all.response_stop <- pbeta(tv, 1 + all.response.sum, 1 + all.nonresponse.sum, lower.tail = FALSE)
  FA_GO_FLAG <- ifelse(all.response_go >= dc_lrv, 1,0)
  FA_STOP_FLAG <- ifelse(all.response_stop <= ar_tv, 2,0)
  ### GO和STOP同时成立，则为STOP(值得商榷，最好不要让这种事情发生，所以要谨慎对待决策空间),还需考虑是否有两个决策都不成立的可能性
  all.response_rag <- ifelse(FA_GO_FLAG*FA_STOP_FLAG==0,max(FA_GO_FLAG,FA_STOP_FLAG),2)
  
  # For interim analysis
  iall.response_go <- pbeta(lrv, 1 + all.response.ia.sum, 1 + all.nonresponse.ia.sum, lower.tail = FALSE)
  iall.response_stop <- pbeta(tv, 1 + all.response.ia.sum, 1 + all.nonresponse.ia.sum, lower.tail = FALSE)
  IA_GO_FLAG <- ifelse(iall.response_go >= dc_lrv, 1,0)
  IA_STOP_FLAG <- ifelse(iall.response_stop <= ar_tv, 2,0)
  iall.response_rag <- ifelse(IA_GO_FLAG*IA_STOP_FLAG==0,max(IA_GO_FLAG,IA_STOP_FLAG),2)
  
  decision_matrix_decis[row,length(pp.monitor)+1] <- ifelse(iall.response_rag==1,TRUE,FALSE)
  decision_matrix_decis[row,length(pp.monitor)+2] <- ifelse(all.response_rag==1,TRUE,FALSE)
  
  # 考虑进中期分析和FA后，记录决策
  current_sample_size <- pp.pass.N[row]
  
  if (current_sample_size==(N_i-1)){
    if (iall.response_rag==1){
      actual_samp_size <- Nt
    } else {
      actual_samp_size <- N_i
    }
  } else {
    actual_samp_size <- current_sample_size
  }
  
  actual_samplesize_vec[row] <- actual_samp_size
  
  }
  
  #########################################################################
  # 结果打包输出
  
  result_package <- list(trueP = true_1,
                         orr = response_rate,
                         source = as.data.frame(source_matrix),
                         ppmatrix = as.data.frame(decision_matrix),
                         decisionmatrix = as.data.frame(decision_matrix_decis),
                         actualsamplesize = actual_samplesize_vec)
  
  return(result_package)
}


# prior
prior_vec <- cmp_CalPrior(prop=0.52,sd=0.31)
priora <- prior_vec[1]
priorb <- prior_vec[2]

sim1 <- SimTrial(Nt,
                 N_i,
                 priora,
                 priorb,
                 true_1,
                 tv,
                 lrv,
                 dc_lrv,
                 ar_tv,
                 Tar,
                 strt,
                 stop,
                 step,
                 PPmon,
                 repsim,
                 seed)

resprate <- sim1$orr
# Kernel Density Plot
resprate_d <- density(resprate) # returns the density data
plot(resprate_d) 