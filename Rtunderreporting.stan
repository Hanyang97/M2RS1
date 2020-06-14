data {
    int<lower = 1> N;               // total number of days
    int<lower=0> I[N];                    // daily new cases
    int<lower=0> T[N];                    // number of tests daily
    int<lower = 0> num_pop;          // total population
    real<lower = 0> infection_potential[N];   // infection potential
}

parameters {
    real<lower =0> Rt[N];                        // R value daily
}

transformed parameters{
    real<lower=0> theta[N];
    for (m in 1:N){
        theta[m] = (Rt[m])*(infection_potential[m])/(num_pop);
    }
}
model {
    for (n in 1:N){
        I[n] ~ binomial(T[n], theta[n]);
    }
}

generated quantities{
    int incidents_rep[N];                            // for posterior predictive checks
    for (i in 1:N){
        incidents_rep[i] = binomial_rng(T[i], theta[i]);
  }
}
