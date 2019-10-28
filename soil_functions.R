

#adapted from ciTools to work with my model matrix
my_ci_survreg<-function(fit,mat,alpha=0.05){
  # 
  crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
  #
  #variance matrix
  cov_mat <- vcov(fit)
  
  #estimates
  beta <- coef(fit)
  scale <- fit$scale
  #
  nPred <- dim(mat)[1]
  seYhat <- rep(NA, nPred)
  for (i in 1:nPred) {
    d_g_beta <- c(exp(mat[i, ] %*% beta)) * mat[i, ] * gamma(1 + scale)
    d_g_delta <- exp(mat[i, ] %*% beta) * digamma(1 + scale) * gamma(1 + scale) * scale
    d_g_vec <- c(d_g_beta, d_g_delta)
    seYhat[i] <- sqrt(t(d_g_vec) %*% cov_mat %*% d_g_vec)
  }
  # distr <- fit$dist
  pred=exp(mat %*% beta) * gamma(1 + scale)
  # pred <- ciTools:::calc_surv_mean(mat = mat, distr = distr, beta = beta, 
  #                        scale = scale)
  w <- exp(crit_val * seYhat/pred)
  mat=as.data.frame(mat)
  mat$pred <- as.vector(pred)
  mat$lwr <- as.vector(pred/w)
  mat$upr <- as.vector(pred * w)
  return(mat)
}



library(survival)
library(tidyverse)
library(frailtypack)
library(tidyr)

#read in data created by data wrangling
soil_all=read.csv("probes/all_probes.csv")
soil=soil_all[soil_all$counted,-1]

#take out "unknown" values
soil=soil[soil$Impact!="Unknown",] #took this out
#specify random effect
soil$re=with(soil,paste0(Swamp,Probe_no))#random effect
# soil$re=with(soil,paste0(Probe_no))#random effect

#create cvariable for ratio of swamp and catchment areas
soil$area_ratio=soil$swamp_area/soil$catchment_area
soil$lsoil_moist_day_before=log(soil$soil_moist_day_before)
soil$ln.days.rain=log(soil$n.days.rain)
soil$ltotal.rain.volume=log(soil$total.rain.volume)
soil$laverage.days.no.rain=log(soil$average.days.no.rain)
soil_narrow=dplyr::select(soil,days.above.50pcSM,da50C,days.above.25pcSM,da25C,
                   days.above.75pcSM,da75C,Impact,veg_type,re,Probe_depth,
                   lsoil_moist_day_before,ltotal.rain.volume,laverage.days.no.rain,
                   scarpdist,temp,precip,area_ratio,Swamp,Probe_no,date)

#scale these variables to stabalize models
soil_narrow[,11:17]=apply(soil_narrow[,11:17],2,scale)

#make sure these are factors
soil_narrow$Impact=factor(soil_narrow$Impact)
soil_narrow$veg_type=factor(soil_narrow$veg_type)

soil_narrow$veg_type=relevel(soil_narrow$veg_type,ref = "bt")

#the code for this function I'm using in R isn't quite working
#so here I hack it a bit
A=model.matrix(~veg_type*Impact,data=soil_narrow)[,-1]
veg_type_=A[,1:3]
Impact_=A[,4]
veg_type_Impact_=A[,5:7]
U_subdat=subset(soil_narrow, soil_narrow$Impact=="Undermined")
R_subdat=subset(soil_narrow, soil_narrow$Impact=="Reference")

#prediction data
nd=cbind(soil_narrow[,11:17]%>%summarise_all(mean),
         expand.grid(Probe_depth=c(20),
               veg_type=levels(soil_narrow$veg_type),
               Impact=levels(soil_narrow$Impact)))


