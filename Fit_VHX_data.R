# Load necessary R packages, these need to be installed beforehand with the function install.packages
library(Rcpp)
library(gplots)
library(dplyr)
library(doParallel)
library(compiler)

# Source R code for fitting 

source('functions.R')
source('compute_relapseTimesR.R')
source('compute_ParasiteProfile.R')

pkdata = readxl::read_excel('VHX_PKPD_data.xlsx', sheet='VHX_PKPD_data')
hist(pkdata$Time_until_Relapse, breaks = seq(0,60,by=1), xlim=c(0,60),
     main='Relapse times in CQ treated patients', xlab = 'days')

N = length(pkdata$Time_until_Relapse)

ts = 0:60
rhos = t(apply(pkdata, 1, function(x) 10^biexp_f(ts, log10_beta = x['log10_beta'], 
                                                 p = x['theta_p'], alpha1 = x['alpha1'],
                                                 alpha2 = x['alpha2'])))

data_XY = list(relapse_time = pkdata$Time_until_Relapse, drug_concentrations=rhos,
               patency_pcount = pkdata$Parasite_Biomass)

data_XY$patency_pcount[is.na(data_XY$patency_pcount)] = 10^mean(log10(data_XY$patency_pcount),na.rm=T)

#####################################################################################################################################
#####################################################################################################################################


# Run the approximation algorithm. This can be resource intensive (depending on K)
K = 10^4
res = abc_simulation(K=K, data_XY = data_XY, random_gen = 'unif')

save(res, file = 'ABCsims_output.RData')


#####################################################################################################################################
#####################################################################################################################################
# Look at model output

load('ABCsims_output.RData')
K = length(res$scores)

hist(res$scores, main='distance scores', xlab='Sum of L1 distance')
####

################################################################################################################
################################################################################################################

par(las=1, bty='n', mfrow=c(1,2))

# Choose nearest parameter set
ind = res$scores > -3000
sum(ind)/N

###################################################### Log k #### 


hist(res$thetas$log_k[ind], main='', xlab='log k', yaxt='n', ylab='', col='grey', freq = F, xlim=c(0,7))
lines(seq(0,10,by=.1), dnorm(seq(0,10,by=.1), mean=1.5, sd=1.5), lwd=3, col='blue')
abline(v=median(res$thetas$log_k[ind]), col='red', lwd=3, lty=2)
log_k_hat = median(res$thetas$log_k[ind])

###################################################### Log EC 50 #### 
hist(res$thetas$log_EC_50[ind], main='', xlab='log EC50', yaxt='n', ylab='', col='grey',freq = F, xlim=c(3,8))
abline(v=median(res$thetas$log_EC_50[ind]),col='red',lwd=3,lty=2)
lines(seq(0,10,by=.1), dnorm(seq(0,10,by=.1), mean=5.5, sd=1.5), lwd=3, col='blue')
logEC50_hat = median(res$thetas$log_EC_50[ind])#sm_avEC50$x[which.max(sm_avEC50$y)]

mic_hat = mic_inpute(logEC50 = logEC50_hat, logk = log_k_hat)
exp(logEC50_hat)
exp(log_k_hat)

# Posterior credible intervals
quantile(mic_inpute(logEC50 = res$thetas$log_EC_50[ind], logk = res$thetas$log_k[ind]), probs=c(0.025,.975))
quantile(exp(res$thetas$log_EC_50[ind]), probs=c(0.1,.9))
quantile(exp(res$thetas$log_k[ind]), probs=c(0.1,.9))




############################################################################################################
#### Plot the concentration-response curve ####
############################################################################################################
par(las=1,bty='n', mfrow=c(1,2))
quantiles = log10(apply(rhos, 2, quantile, probs=c(.05,.5,.95)))
plot(0:60, quantiles[2,], type='l',lwd=3,yaxt='n',
     ylab = 'CQ concentration (ug/L)', xlab='Days post treatment',
     main='Fitted CQ concentrations', ylim=c(log(2),3))
axis(2, at = c(1,2,3), labels = c(10,100,1000))
lines(0:60, quantiles[1,], type='l',lwd=2,lty=2)
lines(0:60, quantiles[3,], type='l',lwd=2,lty=2)
abline(h=log10(mic_hat),col='red', lty=2,lwd=3)
cross1 =  (0:60)[which.min(abs(log10(mic_hat) - quantiles[1,]))]
cross2 =  (0:60)[which.min(abs(log10(mic_hat) - quantiles[3,]))]

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


##################################################################################################################################################################
##### What happens if the MIC doubles?? ##################################################################################################################
##################################################################################################################################################################
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




par(bty='n',las=1)
hist(data_XY$relapse_time, freq=F, breaks = seq(0, 60, by=2),
     main='Distribution of patent relapse times',
     xlab='days post chloroquine treatment',ylab='', yaxt='n',
     density = 10, angle=60)

lines(density(current_bunching[!is.na(current_bunching)]),lwd=4)

round(median(current_bunching,na.rm=T) - median(doubleMIC_bunching,na.rm=T))

lines(density(doubleMIC_bunching[!is.na(doubleMIC_bunching)]), col='red',
      lwd=4, lty=2)


legend('topleft', lwd=3, col=c('black','red'), lty=c(1,2,2), cex=1.1,
       legend = c('Predictions for estimated MIC', 
                  'Predictions for 2 fold increase in MIC'),bty='n')


