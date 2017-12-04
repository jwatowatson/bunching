
# This is the workhorse function for the ABC method
# It forward simulates the fates of N infections given the latent merozoite release time 
# and the observed total parasite biomass of patent infection


compute_relapse_TimesR = function(drug_concentrations,  
                                  log_EC_50, log_k,
                                  merozoite_release,
                                  patency_pcount,
                                  Tmax = 60){
  N = length(merozoite_release)
  N_merozoites = 10000
  log_threshold = log(sqrt(0.1));
  max_effect = 10^3
  
  recrudescence_times = array(dim=N)
  k = exp(log_k)
  
  # iterate through all subjects
  for(i in 1:N){
    
    log_parasitaemia = log(N_merozoites);
    days = merozoite_release[i];
    log_patency_threshold = log( patency_pcount[i] );
    
    while(log_parasitaemia < log_patency_threshold & log_parasitaemia >= 0 & days < Tmax){
      
       drug_level = drug_concentrations[i,days]; 
      
      
      dose_response = sqrt( max_effect / (1 + exp(-k*( log(drug_level) - log_EC_50))) );
      
      log_parasitaemia = log_parasitaemia - max(log_threshold, log(dose_response));
      days = days + 1;
      
    }

    if(log_parasitaemia <= 0){
      recrudescence_times[i] = NA;
    } else{
      recrudescence_times[i] = days;
    }
  }
  return(recrudescence_times)
}

# compile function
cmpfun(compute_relapse_TimesR)