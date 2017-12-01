compute_relapse_TimesR = function(drug_concentrations,  
                                  log_EC_50, log_k,
                                   merozoite_release,
                                   patency_pcount){
  N = length(merozoite_release)
  N_merozoites = 10000
  log_threshold = log(sqrt(0.1));
  Tmax = 150
  max_effect = 10^3
  
  recrudescence_times = array(dim=N)
  #parasite_profiles = array(dim = c(N, 150))
  k = exp(log_k)
  
  for(i in 1:N){
    
    log_parasitaemia = log(N_merozoites);
    days = merozoite_release[i];
    log_patency_threshold = log( patency_pcount[i] );
    
    while(log_parasitaemia < log_patency_threshold & log_parasitaemia >= 0 & days < Tmax){
      
      #parasite_profiles[i, days] = log_parasitaemia
      
      drug_level = drug_concentrations[i,days]; 
      
      
      dose_response = sqrt( max_effect / (1 + exp(-k*( log(drug_level) - log_EC_50))) );
      
      log_parasitaemia = log_parasitaemia - max(log_threshold, log(dose_response));
      days = days + 1;
      
    }
    #parasite_profiles[i, days] = log_parasitaemia
    
    if(log_parasitaemia <= 0){
      recrudescence_times[i] = NA;
    } else{
      recrudescence_times[i] = days;
    }
  }
  #return(list(times =recrudescence_times, prof=parasite_profiles));
  return(recrudescence_times)
}