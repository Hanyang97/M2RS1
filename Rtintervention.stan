data {
  int<lower=0> N;                                     // total no. of days 
  int<lower=0> incidents[N];                          // daily incidents
  vector<lower=0>[N] infec_pot;                       // infection potential in vector form
  int<lower=0> intervention_times;                    // total no. of interventions
  matrix[N,intervention_times] intervention;          // intervention matrix containing 0 and 1
  real<lower=0> alphasd;                            // standard deviation for normal distribution for alpha_k 
  real<lower=0> betasd;                             // standard deviation for normal distribution for beta_l
  real alpha0mean;                                  // mean for alpha_0
  real alpha0sd;                                   // s.d. for alpha_0
}

parameters {
  vector[N] alpha;                         // vecotr contain alpha value
  vector[intervention_times] beta;      // vector contain beta value
}

transformed parameters {
  vector[N] Rt;                             // vector Rt
  real<lower=0> It[N];                      // real array approximated It
  Rt = exp(alpha + intervention * beta);                          // calculate Rt (+ intervention * beta)
  It = to_array_1d(Rt .* infec_pot);        // calculate average It and transform to array
}

model {
  alpha[1] ~ normal(alpha0mean,alpha0sd);                        // first alpha value follows Normal
  for (k in 2:N){          
    alpha[k] ~ normal(alpha[k-1], alphasd);     // alpha_t | alpha_t-1 ~ N(alpha_{t-1}, alphasd^2)
  }
  for (l in 1:intervention_times){
    beta[l] ~ normal(0, betasd);                    // beta ~ N(0, betasd^2)
  }
  for (n in 1:N){              
    incidents[n] ~ poisson(It[n]);
  }
}

generated quantities{
  int incidents_rep[N];                            // for posterior predictive checks
  for (i in 1:N){
    incidents_rep[i] = poisson_rng(It[i]);
  }
}

