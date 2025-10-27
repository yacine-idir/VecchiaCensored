functions {
  vector get_B_F(matrix mini_cov) {
    int l = rows(mini_cov);
    vector[l] result;
    if (l == 1) {
      result[1] = mini_cov[1, 1];
    } else {
      // More efficient way to compute regression coefficients
      row_vector[l-1] Bsi = mdivide_right_spd(mini_cov[1, 2:l], mini_cov[2:l, 2:l]);
      // Compute conditional variance more efficiently using quad_form
      real Fsi = mini_cov[1, 1] - Bsi * mini_cov[2:l, 1];
      // Store results
      result[1:(l-1)] = to_vector(Bsi);
      result[l] = Fsi;
    }
    return result;
  }
  
  real vecchia_lpdf(vector Z, vector sigma, vector phi, real sigma_e, 
                    array[,,] real A, matrix dist, vector beta, 
                    array[] int nb_NN, array[,] int NN, matrix X, 
                    array[] int uncensored_idx){
    int n = size(uncensored_idx);
    int p = size(sigma);
    real my_target = 0;
    real expectation;
    
    // Pre-compute X * beta once
    vector[n] Xbeta = X * beta;
    
    for (i in 1:n) {
      int number_of_neigh = nb_NN[i];
      
      // Extract neighbor indices
      array[number_of_neigh] int neigh;
      for (j in 1:number_of_neigh) {
        neigh[j] = NN[i, j];
      }
      
      // Compute the covariance submatrix directly (only what we need)
      matrix[number_of_neigh, number_of_neigh] minicov = rep_matrix(0.0, number_of_neigh, number_of_neigh);
      
      for (j in 1:number_of_neigh) {
        for (k in j:number_of_neigh) {  // Only upper triangle
          real cov_val = 0.0;
          real dist_jk = dist[neigh[j], neigh[k]];
          
          // Sum over all covariance components
               for (l in 1:p) {
            cov_val += sigma[l] * (exp(-dist_jk/phi[l])) * A[l, neigh[j], neigh[k]];
          }
          

          
          minicov[j,k] = cov_val;
         minicov[k,j] = cov_val;  // Symmetry
        }
      }
      minicov = add_diag(minicov, sigma_e);
      vector[number_of_neigh] bf = get_B_F(minicov);
      
      
      if (number_of_neigh == 1) {
      expectation = Xbeta[i];
      } else {
      expectation = Xbeta[i] + dot_product(bf[1:(number_of_neigh-1)],Z[neigh[2:number_of_neigh]] - Xbeta[neigh[2:number_of_neigh]]);
      }
      
      
      
      real cond_var = bf[number_of_neigh];

// Check for numerical issues
if (is_nan(expectation) || is_inf(expectation)) {
  return negative_infinity();  // Reject this parameter proposal
}

if (cond_var <= 0 || is_nan(cond_var)) {
  return negative_infinity();  // Reject this parameter proposal
}
      
      
      if (uncensored_idx[i] == 0) {
      my_target += normal_lpdf(Z[i] | expectation, sqrt(bf[number_of_neigh]));
      } else {
      my_target += normal_lcdf(Z[i] | expectation, sqrt(bf[number_of_neigh]));
      }
      
    }
    return my_target;
  }
}

data {
  int<lower=0> N;
  vector[N] Z;
  int M;
  matrix[N, N] dist;
  int<lower=0> p; 
  int<lower=0> pX; 
  matrix[N, pX] X; 
  array[N] int uncensored_idx;
  array[N, M] int NN;
  array[N] int nb_NN;
  vector[pX] means_prior_means;
  vector[pX] var_prior_means;
  real mean_prior_sigma;
  real var_prior_sigma;
  real mean_prior_phi;
  real var_prior_phi;
  real mean_prior_tau;
  real var_prior_tau;
  array[p, N, N] real A;
}

parameters {
  vector[pX] alpha;
  vector<lower=0.01>[p] sigma;  // Slightly higher lower bound for numerical stability
  vector<lower=0.01>[p] phi;  // Bounded range for better mixing
  real<lower=0.01, upper=100> sigma_e;  // Bounded noise for better mixing
}

model {
  // More informative priors for better mixing
  alpha ~ normal(means_prior_means, sqrt(var_prior_means));
  sigma ~lognormal(log(mean_prior_sigma),sqrt(var_prior_sigma));
  phi ~lognormal(log(mean_prior_phi),sqrt(var_prior_phi));
  sigma_e ~ lognormal(log(mean_prior_tau),sqrt(var_prior_tau));
  
  // Likelihood - no full covariance matrix needed!
  Z ~ vecchia(sigma, phi, sigma_e, A, dist, alpha, nb_NN, NN, X, uncensored_idx);
}

