// detection_efficiency_script
// this script applies to all tributaries, such that only one of these models must be run 
// and the outputs can be used as inputs for each of the DPS-level primary models.
// This stan model must be run prior to the primary stan model, in order to get posteriors for 
// detection efficiencies in tributaries that are then passed as data to the primary stan model.

data {
  // Fish detections at upstream and river mouth arrays
  // number of observations
  int N;
  // response (detected or not detected)
  int y[N];
  
  // total number of eras
  int J;
  
  // total number of tributaries
  int K;
  
  // The design matrix, X; populated with discharge and intercept terms for eras
  matrix[N, J+K] X;

}


parameters {
  // Initialize era (categorical site configuration) terms - alphas
  vector[J] alpha;
  
  // Initialize discharge slope terms (betas)
  vector[K] beta;
  

}

// this is where we estimate eta, our linear predictor
transformed parameters {
  // Make one vector that contains all parameters
  vector[J + K] params;
  
  // concatenate alphas and betas into params vector
  params[1:J] = alpha;
  params[J+1:J+K] = beta;
  
  // calculate linear predictor
  vector[N] eta;
  eta = X * params;
  
}


model {
  // priors
  alpha ~ normal(0, 5); // this is the prior on the intercepts, aka the inverse logit of the 
  //detection probability at zero discharge (a biologically nonsensical parameter)
  // note that a value of 5 gives you 99%; a value of 0 gives you 0.5, and a value of -5 gives you
  // 0.6%. So I think normal(0,5) is reasonable.
  beta ~ normal(0, 2); // this is the prior on the slopes; we would expect it to be around 0
  
  
  
  y ~ bernoulli_logit(eta);
}
