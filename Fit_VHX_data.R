## Use the AS data to estimate the prior on the time of merozoite release ##
rm(list=ls())
library(Rcpp)
library(gplots)
library(dplyr)
library(boot)
library(doParallel)
setwd('D:/Dropbox/MORU/MICmalaria/Bunching/data/VHX study Maesot/')
setwd("~/Dropbox/MORU/MICmalaria/Bunching/code")
#Rcpp::sourceCpp('D:/Dropbox/MORU/MICmalaria/Bunching/code/compute_relapseTimes.cpp')

source('~/Dropbox/MORU/MICmalaria/Bunching/code/functions.R', echo=TRUE)
source('~/Dropbox/MORU/MICmalaria/Bunching/code/compute_relapseTimesR.R')
source('~/Dropbox/MORU/MICmalaria/Bunching/code/compute_ParasiteProfile.R')

pkdata = read.csv('~/Dropbox/MORU/MICmalaria/Bunching/data/VHX study Maesot/pk_data_cleaned.csv')

hist(pkdata$timetorelapse, breaks = seq(0,400,by=1), xlim=c(0,100),
     main='Relapse times in CQ treated patients', xlab = 'days')

# Only select patients with relapse before 60 days
pkdata = filter(pkdata, timetorelapse<60)
N = length(pkdata$timetorelapse)

# write_csv(x = pkdata[,c("log10pk","timetorelapse","theta_p","alpha1","alpha2","log10_beta","patency_pcount")], 
#           path = 'D:/Dropbox/MORU/MICmalaria/Bunching/data/DataForModel.csv')

ts = 0:150
rhos = t(apply(pkdata, 1, function(x) 10^biexp_f(ts, log10_beta = x['log10_beta'], 
                                                 p = x['theta_p'], alpha1 = x['alpha1'],
                                                 alpha2 = x['alpha2'])))

data_XY = list(relapse_time = pkdata$timetorelapse, drug_concentrations=rhos,
               patency_pcount = pkdata$patency_pcount)

data_XY$patency_pcount[is.na(data_XY$patency_pcount)] = 10^mean(log10(data_XY$patency_pcount),na.rm=T)

#####################################################################################################################################
#####################################################################################################################################
K = 10^5
# Run the approximation algorithm
res=abc_simulation(K=K, data_XY = data_XY, random_gen = 'unif')
res_weib=abc_simulation(K=K, data_XY = data_XY, random_gen = 'weibull')
load('ABCsims_unifLatent.RData')
load('ABCsims_weibullLatent.RData')


#####################################################################################################################################
#####################################################################################################################################
load('ABCsims_unifLatent.RData')
load('ABCsims_weibullLatent.RData')

res = res_weib

K = length(res$scores)

hist(res$scores)
####
pdf('~/Dropbox/MORU/MICmalaria/Bunching/code/ABC_algorithm_results.pdf')
par(las=1, bty='n', mfrow=c(1,2))
### Log k ###
ind = res$scores > -2600
sum(ind)/N
# dens_log_k = kde2d(res$thetas$log_k[ind], res$scores[ind]/100)
# contour(dens_log_k, nlevels = 15, labels = '',  ylab ='L1 score', xlab='log k')
hist(res$thetas$log_k[ind], main='', xlab='log k', yaxt='n', ylab='', col='grey', freq = F)
lines(seq(0,5,by=.1), dnorm(seq(0,5,by=.1), mean=1.9, sd=0.75), lwd=3, col='blue')
abline(v=median(res$thetas$log_k[ind]), col='red', lwd=3, lty=2)
# sm_av_k =wapply(x = res$thetas$log_k[ind], y = res$scores[ind]/100, fun = mean,
#                 width = .05, n = 500)
# lines(sm_av_k$x, sm_av_k$y, lwd=3,col='red')
log_k_hat = median(res$thetas$log_k[ind])#sm_av_k$x[which.max(sm_av_k$y)]

# bsresk = boot(data = cbind(res$thetas$log_k,res$scores/100),
#               statistic = find_smooth_max, R = 200, 
#               parallel = 'multicore', ncpus = 4)
# abline(v=quantile(bsresk$t, probs = c(.025,.975)), col='red', lty=2,lwd=3)
# legend('topleft', col='red',lty=c(1,2),lwd=3, bty='n',
#        legend = c('Estimated expected score', '95% B.I. of maximum score'))

### Log EC 50 #### 

# densEC50 = kde2d(res$thetas$log_EC_50[ind], res$scores[ind]/100)
# contour(densEC50, nlevels = 15, labels = '', ylab ='L1 score', xlab='log EC50')
hist(res$thetas$log_EC_50[ind], main='', xlab='log EC50', yaxt='n', ylab='', col='grey',freq = F)
abline(v=median(res$thetas$log_EC_50[ind]),col='red',lwd=3,lty=2)
lines(seq(0,10,by=.1), dnorm(seq(0,10,by=.1), mean=5.5, sd=0.75), lwd=3, col='blue')
# sm_avEC50=wapply(x = res$thetas$log_EC_50, y = res$scores/100, fun = mean,
#                  width = .05, n = 500)
# lines(sm_avEC50$x, sm_avEC50$y, lwd=3,col='red')
logEC50_hat = median(res$thetas$log_EC_50[ind])#sm_avEC50$x[which.max(sm_avEC50$y)]

# bsresEC50 = boot(data = cbind(res$thetas$log_EC_50,res$scores/100),
#                  statistic = find_smooth_max, R = 200, 
#                  parallel = 'multicore', ncpus = 4)
# abline(v=quantile(bsresEC50$t, probs = c(.025,.975)), col='red', lty=2,lwd=3)

dev.off()

# EC50s = seq(quantile(bsresEC50$t, probs = c(.025,.975))[1],
#             quantile(bsresEC50$t, probs = c(.025,.975))[2], length.out = 100)
# slopes = seq(quantile(bsresk$t, probs = c(.025,.975))[1],
#              quantile(bsresk$t, probs = c(.025,.975))[2], length.out = 100)
# mics = array(dim=c(100,100))
# for(i in 1:100){
#   for(j in 1:100){
#     mics[i,j] = mic_inpute(logEC50 = EC50s[i], logk = slopes[j])
#   }
# }
# range(mics)
mic_hat = mic_inpute(logEC50 = logEC50_hat, logk = log_k_hat)
exp(logEC50_hat)
exp(log_k_hat)


#### Plot the concentration-response curve ####
pdf('~/Dropbox/MORU/MICmalaria/Bunching/code/CQ_MIC_results.pdf')
par(las=1,bty='n', mfrow=c(1,2))
quantiles = log10(apply(rhos, 2, quantile, probs=c(.05,.5,.95)))
plot(0:80, quantiles[2,1:81], type='l',lwd=3,yaxt='n',
     ylab = 'CQ (ug/L)', xlab='days post treatment',
     main='Fitted CQ concentrations', ylim=c(log(2),3))
axis(2, at = c(1,2,3), labels = c(10,100,1000))
lines(0:80, quantiles[1,1:81], type='l',lwd=2,lty=2)
lines(0:80, quantiles[3,1:81], type='l',lwd=2,lty=2)
abline(h=log10(mic_hat),col='red', lty=2,lwd=3)
cross1 =  (0:80)[which.min(abs(log10(mic_hat) - quantiles[1,1:81]))]
cross2 =  (0:80)[which.min(abs(log10(mic_hat) - quantiles[3,1:81]))]

abline(v=c(cross1, cross2),lty=2)


cq_rho = 10^seq(1, 3, by=.01)
yys = dose_response(cq_rho, logk=log_k_hat, logEC50 = logEC50_hat)
plot(log10(yys), log10(cq_rho), yaxt='n', type='l',  ylim=c(log(2),3),
     xlab ='Parasite Reduction Ratio', ylab='CQ concentration (ug/L)',lwd=3,
     main='Response - Concentration', xaxt='n')
axis(2, at=c( 1,2,3), labels = c(10,100,1000))
axis(1, at= -1:3, labels = c(0.1, 1, 10, 100, 1000))
lines(c(-2,0), log10(c(mic_hat,mic_hat)),col='red',lwd=3, lty=2); 
lines(c(0,0),log10(c(1,mic_hat)),col='red',lwd=3, lty=2)
text(x = .3, y = 1, labels = 'MIC',lwd=3,col='black', cex=.9)

IC50_hat =  cq_rho[which.min(abs(0.2 - yys))]
lines(c(-2,log10(0.2)), log10(c(IC50_hat,IC50_hat)),col='red',lwd=3, lty=2); 
lines(log10(c(0.2,0.2)),log10(c(1,IC50_hat)),col='red',lwd=3, lty=2)
text(x = -.42, y = 1, labels = expression('IC'[50]),lwd=3,col='black',cex=.9)


EC50_hat =  cq_rho[which.min(abs(500 - yys))]
lines(c(-2,log10(500)), log10(c(EC50_hat,EC50_hat)),col='red',lwd=3, lty=2); 
lines(log10(c(500,500)),log10(c(1,EC50_hat)),col='red',lwd=3, lty=2)
text(x = 2.3, y = 1, labels = expression('EC'[50]),lwd=3,col='black',cex=.8)


dev.off()


latent_vars = round(rtrunc(N, spec='weibull', b=(data_XY$relapse_time-6),
                           shape = 2.12, scale = 19.86))
res_mat = compute_ParasiteProfiles(drug_concentrations = rhos,
                                   log_EC_50 = logEC50_hat, 
                                   log_k = log_k_hat,
                                   merozoite_release = latent_vars,
                                   patency_pcount = data_XY$patency_pcount)
res_mat = log10(exp(res_mat))
par(mfrow=c(1,2))
plot(NA,NA, xlim=c(0,60), ylim = c(0,10), xlab='days post treatment',
     ylab = 'log10 Parasite load')
for(indv in 1:100){
  ind = !is.na(res_mat[indv,])
  lines(which(ind), res_mat[indv, ind])
}

times = compute_relapse_TimesR(drug_concentrations = rhos,
                               log_EC_50 = logEC50_hat, 
                               log_k = log_k_hat,
                               merozoite_release = latent_vars,
                               patency_pcount = data_XY$patency_pcount)
plot(density(times[!is.na(times)]),main='',xlab='days',yaxt='n',
     lwd=2,col='blue', ylab='')
lines(density(data_XY$relapse_time),lwd=2,col='red')
legend('topright', col=c('red','blue'), lwd=2, legend=c('Observed','Simulated'),bty='n')


##### What happens if the MIC doubles?? ######
par(mfrow=c(1,1))
NNN = 20000
ind_bb = sample(1:nrow(rhos), NNN, replace = T)
Large_trial = rhos[ind_bb,]
latent_vars = rweibull(NNN, shape = 2.12, scale = 19.8)
parasites = 10^rnorm(NNN, mean = mean(log10(data_XY$patency_pcount)), sd = sd(log10(data_XY$patency_pcount)))


logEC50s = seq(logEC50_hat/2, 2*logEC50_hat, by = 0.001)
mics = mic_inpute(logEC50 = logEC50s, logk = log_k_hat)
logEC50_doubleMIC = logEC50s[which.min(abs(mics - 2*mic_hat))]
logEC50_halfMIC = logEC50s[which.min(abs(mics - mic_hat/2))]


current_bunching = compute_relapse_TimesR(drug_concentrations = Large_trial,
                                          log_EC_50 = logEC50_hat, 
                                          log_k = log_k_hat,
                                          merozoite_release = latent_vars, 
                                          patency_pcount = parasites)
doubleMIC_bunching = compute_relapse_TimesR(drug_concentrations = Large_trial,
                                            log_EC_50 = logEC50_doubleMIC, 
                                            log_k = log_k_hat,
                                            merozoite_release = latent_vars, 
                                            patency_pcount = parasites)
halfMIC_bunching = compute_relapse_TimesR(drug_concentrations = Large_trial,
                                          log_EC_50 = logEC50_halfMIC, 
                                          log_k = log_k_hat,
                                          merozoite_release = latent_vars, 
                                          patency_pcount = parasites)

pdf('DoublingMICchloroquine.pdf')
par(bty='n',las=1)
hist(data_XY$relapse_time, freq=F, breaks = seq(0, 60, by=2),
     main='Distribution of patent relapse times',
     xlab='days post chloroquine treatment',ylab='', yaxt='n',
     density = 10, angle=60)

lines(density(current_bunching[!is.na(current_bunching)]),lwd=4)

round(median(current_bunching,na.rm=T) - median(doubleMIC_bunching,na.rm=T))

lines(density(doubleMIC_bunching[!is.na(doubleMIC_bunching)]), col='red',
      lwd=4, lty=2)

lines(density(halfMIC_bunching[!is.na(halfMIC_bunching)]), col='blue',
      lwd=4, lty=2)
# relapse_GORDON_data = c(rep(4, 1),rep(5, 1),rep(6, 18),rep(7, 36),rep(8, 18),rep(9, 23))
# hist(7*relapse_GORDON_data - 3.5, xlim=c(1, 12), freq=F, border = 'blue', add=T,
#      breaks = seq(21,84,by=7), density = 10, angle=-30, lwd=3)
# 

legend('topleft', lwd=3, col=c('black','red','blue'), lty=c(1,2,2),
       legend = c('Predictions for estimated MIC', 
                  'Predictions for 2 fold increase in MIC',
                  'Predictions for 2 fold decrease in MIC'),bty='n')

dev.off()


relapse_GORDON_data = c(rep(4, 1),rep(5, 1),rep(6, 18),rep(7, 36),rep(8, 18),rep(9, 23))
hist(7*relapse_GORDON_data - 3.5, xlim=c(1, 12), freq=F, col='blue', add=T)
hist(data_XY$relapse_time/7, add=T, border='red',freq = F)
