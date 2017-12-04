# This function computes the parasite biomass profiles for a given set of PD parameters and PK profiles
# drug_concentrations : N x Tmax matrix; N is number of subjects, Tmax is maximum number of days to forward simulate
# log_EC_50 : log value of the EC50 parameter of PD model
# log_k : log value of the slope parameter k in PD model
# merozoite_release : days at which merozoites erupt from liver (start of blood stage infection)
# patency_pcount : parasite biomass at which the infection becomes patent

# This function returns parasite profiles : N x Tmax matrix
compute_ParasiteProfiles = function(drug_concentrations,  
                                    log_EC_50, log_k,
                                    merozoite_release,
                                    patency_pcount){
  N = length(merozoite_release)
  N_merozoites = 10000
  log_threshold = log(sqrt(0.1));
  Tmax = 150
  max_effect = 10^3
  
  parasite_profiles = array(dim = c(N, 150))
  k = exp(log_k)
  
  for(i in 1:N){
    
    log_parasitaemia = log(N_merozoites);
    days = merozoite_release[i];
    log_patency_threshold = log( patency_pcount[i] );
    
    while(log_parasitaemia < log_patency_threshold & log_parasitaemia >= 0 & days < Tmax){
      
      parasite_profiles[i, days] = log_parasitaemia
      
      drug_level = drug_concentrations[i,days]; 
      
      
      dose_response = sqrt( max_effect / (1 + exp(-k*( log(drug_level) - log_EC_50))) );
      
      log_parasitaemia = log_parasitaemia - max(log_threshold, log(dose_response));
      days = days + 1;
      
    }
    parasite_profiles[i, days] = log_parasitaemia
    
    
  }
  return(parasite_profiles);
  
}