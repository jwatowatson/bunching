library(MASS)
library(truncdist)

biexp_f = function(t, log10_beta, p, alpha1, alpha2){
  beta = 10^log10_beta
  log10(p*beta*exp(alpha1 * t) + (1-p)*beta*exp(alpha2 * t))
}

biexp_f_cov = function(t, log10_beta, p, alpha1, alpha2, theta_age, age){
  beta = theta_age * age + 10^log10_beta
  log10(p*beta*exp(alpha1 * t) + (1-p)*beta*exp(alpha2 * t))
}

surface_area = function(w,h){
  #Du Bois formula
  w^0.425 * h^0.725 * .007184
}


compute_likelihoodR = function(theta, latent_theta, data_all){
  relapse_hat = array(dim=c(10,N))
  for(i in 1:ncol(latent_theta)){
    relapse_hat[i,] = compute_relapse_TimesR(drug_concentrations = data_all$drug_concentrations, 
                                             log_EC_50 = theta['log_EC_50'], 
                                             log_k = theta['log_k'],
                                             merozoite_release = latent_theta[,i], 
                                             patency_pcount = data_all$patency_pcount)
  }
  
  relapse_hat = apply(relapse_hat, 2, mean, na.rm=T)
  relapse_hat[is.nan(relapse_hat)] = 0
  L1_liks = -abs(relapse_hat - data_all$relapse_time)
  
  return(sum(L1_liks))
  
}

sim_prior = function(N){
  a1 = rnorm(N, mean = 5.5, sd = 1.5)
  a2 = rnorm(N, mean = 1.5, sd = 1.5)
  return(list(log_EC_50=a1, log_k=a2))
}



abc_simulation = function(K, data_XY, random_gen = c('unif', 'weibull'), Ncores = 4){
  
  registerDoParallel(cores=Ncores)
  
  if(!(length(data_XY$relapse_time) == length(data_XY$patency_pcount)) | !(length(data_XY$patency_pcount) == nrow(data_XY$drug_concentrations)) ){
    stop('problemo: not the right length inputs')
  }
  
  thetas_PD = sim_prior(K)
  N = length(data_XY$patency_pcount)
  scores <- foreach(i=1:K, .combine = 'c') %dopar% {
    if(random_gen == 'unif'){
      latent_vars = matrix(round(runif(n = 10*N, min=0, max = (data_XY$relapse_time-6))), nrow=N)
    } else if (random_gen == 'weibull'){
      latent_vars1 = matrix(round(rtrunc(10*N, spec='weibull', b=(data_XY$relapse_time-6),
                                         shape = 2.12, scale = 19.86)), nrow=N)
    }
    compute_likelihoodR(theta = c(log_k=thetas_PD$log_k[i],
                                  log_EC_50=thetas_PD$log_EC_50[i]),
                        latent_theta = latent_vars, 
                        data_all = data_XY)
  }
  
  return(list(thetas=thetas_PD, scores = scores))
}

dose_response = function(rho, logk, logEC50, max_effect=10^3, min_effect = 10^(-1)){
  y = min_effect + (max_effect-min_effect) / (1 + exp(-exp(logk)*( log(rho) - logEC50)))
  return(y)
}


mic_inpute = function(logEC50, logk, max_effect=10^3, min_effect = 10^(-1)){
  #exp( logEC50 - log(max_effect-1)/exp(logk) )
  exp(logEC50 - (log(max_effect-1) - log(1-min_effect))/exp(logk))
}

find_smooth_max = function(data, ind){
  
  sm_av=wapply(x = data[ind,1], y = data[ind,2], fun = mean,
               width = .1, n = 500)
  
  max_value = sm_av$x[which.max(sm_av$y)]
  return(max_value)
}

bootstrap_summary = function(xs, ys, BB = 200){
  N = length(xs)
  boot_res = array(dim = c(BB, N))
  for(i in 1:BB){
    ind = sample(1:N, N, replace=T)
    boot_res[i,] = smooth_summary(xs[ind], ys[ind])
  }
}
