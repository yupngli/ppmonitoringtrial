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
  #	N_I	:	中期分析时计划的样本量
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

cmp_CalPrior = cmpfun(CalPrior)
cmp_PredPower = cmpfun(PredPower)

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
## PPmon: 继续前进的最低PP要求
## decision: 决策模式,Dual or Single(futility only)
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
                   PPmon,
                   decision,
                   repsim,
                   seed=1234) { 
  
  set.seed(seed)
  
  # 创建模拟矩阵
  source_matrix <- matrix(data=NA,nrow=repsim,ncol=Nt)
  
  # 创建决策矩阵:记录决策概率
  decision_matrix <- matrix(data=NA,nrow=repsim,ncol=Nt)
  
  # 创建决策矩阵:记录决策
  decision_matrix_decis <- matrix(data=NA,nrow=repsim,ncol=Nt)
  
  #	Response rate的先验分布
  Bpar1 <- c(priora,priorb)
  
  for (row in 1:repsim){
    theta_1 <- rbeta(1, Bpar1[1] + 1, Bpar1[2] + 1)
    # 为全人群Nt产生模拟的response
    all.response <- rbinom(Nt, 1, theta_1)
    source_matrix[row,] <- all.response
  }
  
  #########################################################################
  # phase A & phase B
  #########################################################################
  
  for (row in 1:repsim){
    # 截取前N_i个数据作为Phase C中期分析的依据
    all.response.IA <- source_matrix[row,1:N_i]
    # 数据容器 - 存储PP监测的结果
    pp.monitor <- vector(mode = "numeric", length = stop - strt + 1)
    # 数据容器 - 存储PP监测时，患者数量
    pp.monitor.N <- seq(from=strt,to=stop,by=1)
    # 数据容器 - 存储PP监测时，步数
    h <- seq(from=1,to=(stop - strt + 1),by=1)
    for (g in 1:(stop - strt + 1)) {
      pp.monitor[g] <-
        cmp_PredPower(Nt, strt + (g - 1), round(Nt * Tar, 0), max(sum(all.response.IA[1:(strt + (g - 1))])))
    }
    # PP计算结果填入决策矩阵
    decision_matrix[row,strt:stop] <- pp.monitor
  }
  
  # 截取前N_i个数据作为Phase C中期分析的依据 - retire
  all.response.IA <- all.response[1:N_i]
  
  # 数据容器 - 存储PP监测的结果
  pp.monitor <- vector(mode = "numeric", length = stop - strt + 1)
  
  # 数据容器 - 存储PP监测时，患者数量
  pp.monitor.N <- seq(from=strt,to=stop,by=1)
  
  # 数据容器 - 存储PP监测时，步数
  h <- seq(from=1,to=(stop - strt + 1),by=1)
  
  for (g in 1:(stop - strt + 1)) {
    pp.monitor[g] <-
      cmp_PredPower(Nt, strt + (g - 1), round(Nt * Tar, 0), max(sum(all.response.IA[1:(strt + (g - 1))])))
  }
  
  # 组装pp.monitor与pp.monitor.N成dataframe: pp.df
  pp.df <- data.frame(pp.monitor, pp.monitor.N)
  
  # pp.df中筛选出小于pp标准的，停止当前招募
  pp.stop.df <- pp.df %>% filter(pp.monitor<PPmon)
  
  # 中期分析时,累积有多少人通过了pp监测
  pp.pass.N.temp <- ifelse(length(pp.stop.df$pp.monitor.N) == 0, 0, min(pp.stop.df$pp.monitor.N))
  if (min(pp.df$pp.monitor) >= PPmon){
    pp.pass.N <- N_i 
  } else {
    pp.pass.N <- pp.pass.N.temp
  }
  
  #########################################################################
  
  # calculate the REAL number of corresponders among all patients
  all.response.sum <- sum(all.response)
  
  # calculate the REAL number of NON-corresponders among all patients
  all.nonresponse.sum <- Nt - all.response.sum
  
  # calculate the REAL number of corresponders among patients in Interim Analysis
  all.response.ia.sum <- sum(all.response.IA)
  
  # calculate the REAL number of NON-corresponders among patients in Interim Analysis
  all.nonresponse.ia.sum <- N_i - all.response.ia.sum
  
  #########################################################################
  
  # 	Decide
  
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
  
  #########################################################################
  # 	Assemble the results   
  
  ## Interim Analysis
  #### actual_samp_size: 记录早停时候的样本量
  #### 情况1:IA之前,pp监测中停止,样本量实际为pp.pass.N
  #### 情况2:IA之时停止,样本量实际为N_i
  #### 情况3:不早停，样本量为Nt
  # actual_samp_size <-
  #   ifelse(iall.response_rag == 2, ifelse(pp.pass.N < N_i, pp.pass.N, N_i), Nt)
  actual_samp_size <- pp.pass.N
  if (actual_samp_size == N_i){
    if (iall.response_rag != 2){
      actual_samp_size <- Nt
    }
  }
  
  pre_int1 <- ifelse(pp.pass.N < N_i, "STOP", "GO")
  p <- true_1
  ne <- N_p
  
  # change flag to char
  FA_FL <- recode(all.response_rag, `2` = "STOP", `1` = "GO")
  IA_FL <- recode(iall.response_rag, `2` = "STOP", `1` = "GO")
  
  #	Package results
  tr <-
    data.frame(
      ne,
      p,
      all.response.sum, # the REAL number of corresponders among all patients
      all.nonresponse.sum, # the REAL number of NON-corresponders among all patients
      all.response_go, # boundary for GO at FA
      all.response_stop, # boundary for STOP at FA
      FA_FL, # Decision at FA
      all.response.ia.sum, # the REAL number of corresponders among patients in IA
      all.nonresponse.ia.sum, # the REAL number of NON-corresponders among patients in IA
      iall.response_go, # boundary for GO at IA
      iall.response_stop, # boundary for STOP at IA
      IA_FL, # Decision at IA
      pp.pass.N, # the number of patients before NO Continue decision is made
      pre_int1, # indicate if trial stops before IA
      actual_samp_size # the sample size when trial stops
    )
  
  return(tr)
}