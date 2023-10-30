#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA_INTEGER(disp_model); // 1 = half-normal, 2 = exponential, 3 = half-cauchy
  // DATA_INTEGER(like_model) // 1 = poisson, 2 = negative binomial. Not currently used.
  // DATA_INTEGER(ncreeks);  // ncreeks is the number of release locations (creeks). Not currently used.
  DATA_INTEGER(nrel);  // nrel is the number of tagged fish released
  DATA_INTEGER(nsites); // the number of trapping sites at each location
  DATA_INTEGER(nperiods); // the number of recap events at each location
  DATA_INTEGER(ntraps); // the number of traps at each site
  DATA_MATRIX(countmat); //(1,nsites*nperiods,1,ntraps) // the number of recap events at each location
  DATA_VECTOR(distances); //(1,nsites) // vector of distances from release site to center of recap site
  DATA_VECTOR(times); //(1,nperiods) // vector of times since release
  DATA_SCALAR(kappa);


  // PARAMETER(logit_survival);
  // PARAMETER(log_detectability);
  PARAMETER(sig_disp_alpha);
  PARAMETER(sig_disp_beta);
  PARAMETER(survival);
  PARAMETER(detectability);
  // PARAMETER_VECTOR(sig_disp);
  // PARAMETER(sig_disp);

  // Type survival = exp(logit_survival)/(1+exp(logit_survival));
  // Type detectability = exp(logit_detectability)/(1+exp(logit_detectability));
  // vector<Type> sig_disp(nperiods);
  // for(int period=0; period<nperiods; period++) {
  //   sig_disp(period) = exp(log_sig_disp(period));
  // }

  // Type site_width = full_site_width / Type(2.0);

  matrix<Type> predcount(nperiods, nsites * ntraps + 1);
  matrix<Type> obscount(nperiods, nsites * ntraps + 1);
  // vector<Type> dist_factor(nsites * nperiods * ntraps);
  vector<Type> hazard(nsites);
  matrix<Type> outmat(nsites * nperiods * ntraps, 4);
  vector<Type> pcapture(nperiods);
  vector<Type> sig_disp(nperiods);
  matrix<Type> multinom_probs(nperiods, nsites*ntraps + 1);
  vector<Type> total_caught(nperiods);

  Type f = 0; // objective function value

  int counter = 0;
  for (int i=0; i<nperiods; i++)
  {
    sig_disp(i) = sig_disp_alpha * times(i)/(1 + sig_disp_beta * times(i)); 
    for(int site=0; site<nsites; site++) {
      hazard(site) = detectability * exp(-pow(distances(site)/sig_disp(i), kappa));
    }
    pcapture(i) = 1 - exp(-hazard.sum()); 
    total_caught(i) = 0;

    for (int j=0; j<nsites; j++) {
      for (int k=0; k<ntraps; k++) {
        // prob of capture is evenly divided among traps @ each site
        multinom_probs(i, k + j*ntraps) = (hazard(j)/ntraps)/hazard.sum();  
        obscount(i, k + j * ntraps) = countmat(j+(i*nsites), k);
        total_caught(i) += obscount(i, k + j * ntraps);
      }
    }
    multinom_probs(nsites*ntraps) = (1-pcapture(i)) + pow(survival, times(i));
    obscount(i, nsites*ntraps) = nrel - total_caught(i);

     // predcount(counter) = nrel * pow(survival, times(i)) * pcapture(i) * (hazard(j)/hazard.sum());
     // obscount(counter) = countmat(j+(i*nsites), k);
     f -= dmultinom(obscount.row(i), multinom_probs.row(i), true);

     outmat(counter, 0) = obscount(counter);
     outmat(counter, 1) = predcount(counter);
     outmat(counter, 2) = distances(j); // times(i);
     outmat(counter, 3) = times(i); //distances(j);

     counter++;
    }
   }   
  }
 
  ADREPORT(survival);
  ADREPORT(detectability);
  ADREPORT(sig_disp);
  REPORT(outmat);
  return f;
}
