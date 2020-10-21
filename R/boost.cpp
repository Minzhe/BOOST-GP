#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

static arma::mat K_builder(arma::mat dist, double l);
static double rnorm_trunc(double mu, double sigma, double lower, double upper);
static arma::mat create_ones(int n, int p);
static arma::mat create_eye(int n);
static int num(IntegerVector x, int c);
static arma::uvec which(IntegerVector x, int c);
static arma::uvec which_stop(IntegerVector x);
static arma::uvec as_uvec(int j);
  
// [[Rcpp::export]]
Rcpp::List boost(arma::mat Y, arma::mat dist, IntegerMatrix nei, NumericVector s, int iter, int burn, double init_b_sigma, double init_h, double update_prop) {
  // Read data information
  int n = Y.n_rows;
  int p = Y.n_cols;
  // Calculate the minimum and maximum distances so that to set up prior setting
  arma::mat dist_temp = dist;
  dist_temp.diag() = max(dist);
  double t_min = dist_temp.min();
  double t_max = dist_temp.max();
  
  // Set hyperparameters
  double a_pi = 1;
  double b_pi = 1;
  double a_phi = 0.001; 
  double b_phi = 0.001; 
  double omega = 0.05;
  // double a_omega = 0.1;
  // double b_omega = 2.0 - a_omega;
  double h = init_h;
  double a_sigma = 3.0;
  double b_sigma = init_b_sigma; //change 20200717
  double a_l = t_min/2;
  double b_l = 2*t_max;
  
  // Set algorithm settings
  // int iter = 1000; //change 20200717
  bool nb = true;
  bool zi = true;
  // int burn = 500;
  int E = max(NumericVector::create(1, p*update_prop)); 
  int F = max(NumericVector::create(1, n*update_prop)); 
  double tau_log_phi = 1.0;
  double tau_log_lambda = 0.1; //change 20200717
  double tau_log_l = 0.1;
  // Do not display the warning messages
  std::ostream nullstream(0);
  set_cerr_stream(nullstream);
  
  // Set temporary variables
  int it, i, j, e, f, count = 0, gamma_sum_temp = 0, count_temp, H_temp, n_temp;
  double hastings, phi_temp, loglambda_temp, loglambda_null_temp, gamma_temp, l_temp, max_temp, sum_temp;
  double log_K_det_temp, log_K_det_sign_temp, K_inv_sum_temp, G1_temp, G2_temp;
  IntegerVector n_eff(p);
  IntegerVector H_sum_temp(n, 0);
  NumericVector prob_temp(2);
  NumericVector K_inv_sum(p), log_K_det(p), G1(p), G2(p), loglambda_null(p);
  NumericVector K_inv_sum_full(p), log_K_det_full(p), G1_full(p), G2_full(p), loglambda_null_full(p);
  arma::vec logLambda_temp, logLambda_temp_0;
  arma::mat K, temp, K_temp;
  double accept_phi = 0, accept_lambda = 0, accept_gamma = 0, accept_l = 0;
  double try_phi = 0, try_lambda = 0, try_gamma = 0, try_l = 0;
  
  IntegerMatrix H(n, p);
  NumericVector pi(n, 0.5);
  NumericVector phi(p);
  arma::mat logLambda(n, p);
  IntegerVector gamma(p);
  NumericVector l(p);
  
  NumericMatrix H_ppi(n, p);
  IntegerVector H_sum(iter);
  NumericVector gamma_ppi(p);
  IntegerVector gamma_sum(iter);
  NumericMatrix logBF(iter, p);
  NumericMatrix pi_store(iter, n);
  NumericMatrix phi_store(iter, p);
  NumericVector omega_store(iter);
  IntegerMatrix gamma_store(iter, p);
  NumericMatrix l_store(iter, p);
  
  // Initialization
  double phi_start = 10, loglambda_start = log(mean(mean(Y))), l_start = t_min/2;
  for (j = 0; j < p; j++) {
    for (i = 0; i < iter; i++) {
      logBF(i, j) = 0;
    }
    for (i = 0; i < n; i++) {
      logLambda(i, j) = loglambda_start;
      H_ppi(i, j) = 0;
      if (Y(i, j) != 0) {
        H(i, j) = 0;
      } else {
        H(i, j) = 0;
        H_sum_temp(i) = H_sum_temp(i) + H(i, j);
      }
    }
    gamma_ppi(j) = 0;
    n_eff(j) = num(H(_, j), 0);
    if (nb) {
      phi(j) = phi_start;
    } else {
      phi(j) = -1;
    }
    l(j) = l_start;
    //gamma(j) = rbinom(1, 1, 0.5)(0);
    gamma(j) = 0;
    if (gamma(j) == 1) {
      gamma_sum_temp++;
      K = K_builder(dist, l(j));
      K = K.submat(which(H(_, j), 0), which(H(_, j), 0));
      log_det(log_K_det_temp, log_K_det_sign_temp, K);
      log_K_det(j) = log_K_det_temp*log_K_det_sign_temp;
      solve(temp, K, create_ones(n_eff(j), 1), solve_opts::likely_sympd + solve_opts::no_approx);
      K_inv_sum(j) = (create_ones(1, n_eff(j))*temp).eval()(0, 0);
      solve(temp, K, logLambda.submat(which(H(_, j), 0), as_uvec(j)), solve_opts::likely_sympd + solve_opts::no_approx);
      G1(j) = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*temp).eval()(0, 0);
      G2(j) = (temp.t()*create_ones(n_eff(j), n_eff(j))*temp/(K_inv_sum(j) + 1.0/h)).eval()(0, 0);
    }
    G1_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*logLambda.submat(which(H(_, j), 0), as_uvec(j))).eval()(0, 0);
    G2_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*create_ones(n_eff(j), n_eff(j))*logLambda.submat(which(H(_, j), 0), as_uvec(j))/(1.0*n_eff(j) + 1.0/h)).eval()(0, 0);
    loglambda_null(j) = -log(1.0*n_eff(j) + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0);
  }
  
  // MCMC
  for (it = 0; it < iter; it++) {
    // Update H
    if (zi) {
      for(j = 0; j < p; j++) {
        count_temp = 0;
        for(i = 0; i < n; i++) {
          if(Y(i, j) == 0) {
            if (nb) {
              prob_temp(0) = phi(j)*(log(phi(j)) - log(s(i)*exp(logLambda(i, j)) + phi(j))) + log(1 - pi(i));
            } else {
              prob_temp(0) = -s(i)*exp(logLambda(i, j)) + log(1 - pi(i));
            }
            prob_temp(1) = log(pi(i));
            max_temp = max(prob_temp);
            prob_temp(0) = prob_temp(0) - max_temp;
            prob_temp(1) = prob_temp(1) - max_temp;
            prob_temp(0) = exp(prob_temp(0));
            prob_temp(1) = exp(prob_temp(1));
            sum_temp = prob_temp(0) + prob_temp(1);
            prob_temp(0) = prob_temp(0)/sum_temp;
            prob_temp(1) = prob_temp(1)/sum_temp;
            H_temp = H(i, j);
            H(i, j) = rbinom(1, 1, prob_temp(1))(0);
            H_sum_temp(i) = H_sum_temp(i) - H_temp + H(i, j);
            if(H(i, j) != H_temp) {
              count_temp++;
            }
          }
          // pi(i) = 0.5;
        }
        if (count_temp > 0) {
          n_eff(j) = num(H(_, j), 0);
          if (gamma(j) == 1) {
            K = K_builder(dist, l(j));
            K = K.submat(which(H(_, j), 0), which(H(_, j), 0));
            log_det(log_K_det_temp, log_K_det_sign_temp, K);
            log_K_det(j) = log_K_det_temp*log_K_det_sign_temp;
            solve(temp, K, create_ones(n_eff(j), 1), solve_opts::likely_sympd + solve_opts::no_approx);
            K_inv_sum(j) = (create_ones(1, n_eff(j))*temp).eval()(0, 0);
            solve(temp, K, logLambda.submat(which(H(_, j), 0), as_uvec(j)), solve_opts::likely_sympd + solve_opts::no_approx);
            G1(j) = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*temp).eval()(0, 0);
            G2(j) = (temp.t()*create_ones(n_eff(j), n_eff(j))*temp/(K_inv_sum(j) + 1.0/h)).eval()(0, 0);
          }
          G1_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*logLambda.submat(which(H(_, j), 0), as_uvec(j))).eval()(0, 0);
          G2_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*create_ones(n_eff(j), n_eff(j))*logLambda.submat(which(H(_, j), 0), as_uvec(j))/(1.0*n_eff(j) + 1.0/h)).eval()(0, 0);
          loglambda_null(j) = -log(1.0*n_eff(j) + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0);
        }
      }
    }

    // Update phi
    if (nb) {
      for(j = 0; j < p; j++) {
        phi_temp = exp(r_truncnorm(log(phi(j)), tau_log_phi, log(1), log(100)));
        hastings = 0;
        for(i = 0; i < n; i++) {
          if (H(i, j) == 0) {
            hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + Y(i, j)) - (phi_temp + Y(i, j))*log(phi_temp + s(i)*exp(logLambda(i, j)));
            hastings = hastings - (phi(j)*log(phi(j)) - lgamma(phi(j)) + lgamma(phi(j) + Y(i, j)) - (phi(j) + Y(i, j))*log(phi(j) + s(i)*exp(logLambda(i, j))));
          }
        }
        hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
        hastings = hastings - ((a_phi - 1)*log(phi(j)) - b_phi*phi(j));
        if (it > burn) {
          try_phi++;
        }
        if(hastings >= log(double(rand()%10001)/10000)) {
          phi(j) = phi_temp;
          if (it > burn) {
            accept_phi++;
          }
        }
      }
    }
    
    // Update Lambda
    for (j = 0; j < p; j++) {
      if (gamma(j) == 1) {
        K = K_builder(dist, l(j));
      }
      for(f = 0; f < F; f++) {
        i = rand()%n;
        //if (Y(i, j) == 0)
        //{
        //  loglambda_temp = rnorm(1, 0, tau_log_lambda)(0);
        //}
        //else
        //{
          loglambda_temp = rnorm(1, logLambda(i, j), tau_log_lambda)(0);
        //}
        logLambda_temp = logLambda.col(j);
        logLambda_temp_0 = logLambda_temp.elem(which_stop(nei(i, _)));
        logLambda_temp(i) = loglambda_temp;
        logLambda_temp = logLambda_temp.elem(which_stop(nei(i, _)));
        n_temp = which_stop(nei(i, _)).size();
        hastings = 0;
        if (H(i, j) == 0) {
          if (nb) {
            hastings = hastings + (Y(i, j)*(log(s(i)) + loglambda_temp) - (phi(j) + Y(i, j))*log(phi(j) + s(i)*exp(loglambda_temp)));
            hastings = hastings - (Y(i, j)*(log(s(i)) + logLambda(i, j)) - (phi(j) + Y(i, j))*log(phi(j) + s(i)*exp(logLambda(i, j))));
          } else {
            hastings = hastings + (Y(i, j)*(log(s(i)) + loglambda_temp) - s(i)*exp(loglambda_temp));
            hastings = hastings - (Y(i, j)*(log(s(i)) + logLambda(i, j)) - s(i)*exp(logLambda(i, j)));
          }
        }
        //if (Y(i, j) == 0) {
          if (gamma(j) == 0) {
            G1_temp = (logLambda_temp.t()*logLambda_temp).eval()(0,0);
            G2_temp = (logLambda_temp.t()*create_ones(n_temp, n_temp)*logLambda_temp/(1.0*n_temp + 1.0/h)).eval()(0, 0);
            loglambda_null_temp = -log(1.0*n_temp + 1.0/h)/2.0 - (a_sigma + 1.0*n_temp/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0);
            hastings = hastings + loglambda_null_temp;
            G1_temp = (logLambda_temp_0.t()*logLambda_temp_0).eval()(0,0);
            G2_temp = (logLambda_temp_0.t()*create_ones(n_temp, n_temp)*logLambda_temp_0/(1.0*n_temp + 1.0/h)).eval()(0, 0);
            loglambda_null_temp = -log(1.0*n_temp + 1.0/h)/2.0 - (a_sigma + 1.0*n_temp/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0);
            hastings = hastings - loglambda_null_temp;
          } else {
            K_temp = K.submat(which_stop(nei(i, _)), which_stop(nei(i, _)));
            solve(temp, K_temp, create_ones(n_temp, 1), solve_opts::likely_sympd + solve_opts::no_approx);
            K_inv_sum_temp = (create_ones(1, n_temp)*temp).eval()(0, 0);
            solve(temp, K_temp, logLambda_temp, solve_opts::likely_sympd + solve_opts::no_approx);
            G1_temp = (logLambda_temp.t()*temp).eval()(0, 0);
            G2_temp = (temp.t()*create_ones(n_temp, n_temp)*temp/(K_inv_sum_temp + 1.0/h)).eval()(0, 0);
            hastings = hastings + (-(a_sigma + 1.0*n_temp/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0));
            solve(temp, K_temp, logLambda_temp_0, solve_opts::likely_sympd + solve_opts::no_approx);
            G1_temp = (logLambda_temp_0.t()*temp).eval()(0, 0);
            G2_temp = (temp.t()*create_ones(n_temp, n_temp)*temp/(K_inv_sum_temp + 1.0/h)).eval()(0, 0);
            hastings = hastings - (-(a_sigma + 1.0*n_temp/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0));
          }
        //}
        /*
        if (Y(i, j) == 0)
        {
          hastings = hastings + (-pow(logLambda(i, j), 2)/2.0/pow(tau_log_lambda, 2));
          hastings = hastings - (-pow(loglambda_temp, 2)/2.0/pow(tau_log_lambda, 2));
        }
         */
        if (it > burn) {
          try_lambda++;
        }
        if (hastings >= log(double(rand()%10001)/10000)) {
          logLambda(i, j) = loglambda_temp;
          if (it > burn) {
            accept_lambda++;
          }
        }
      }
      if (gamma(j) == 1) {
        K = K.submat(which(H(_, j), 0), which(H(_, j), 0));
        solve(temp, K, logLambda.submat(which(H(_, j), 0), as_uvec(j)), solve_opts::likely_sympd + solve_opts::no_approx);
        G1(j) = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*temp).eval()(0, 0);
        G2(j) = (temp.t()*create_ones(n_eff(j), n_eff(j))*temp/(K_inv_sum(j) + 1.0/h)).eval()(0, 0);
      }
      G1_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*logLambda.submat(which(H(_, j), 0), as_uvec(j))).eval()(0, 0);
      G2_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*create_ones(n_eff(j), n_eff(j))*logLambda.submat(which(H(_, j), 0), as_uvec(j))/(1.0*n_eff(j) + 1.0/h)).eval()(0, 0);
      loglambda_null(j) = -log(1.0*n_eff(j) + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0);
    }
    
    // Update gamma
    for(e = 0; e < E; e++) {
      j = rand()%p;
      gamma_temp = 1 - gamma(j);
      if(gamma_temp == 0) {// Delete
        hastings = 0;
        // Old
        hastings = -(-log_K_det(j)/2.0 - log(K_inv_sum(j) + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1(j)/2.0 - G2(j)/2.0));
        // New
        hastings = hastings + loglambda_null(j);
        logBF(it, j) = -hastings;
        // Proposal
        hastings = hastings + (-pow(log(l(j)) - log(l_start), 2)/2.0/pow(10.0*tau_log_l, 2));
        // Prior
        hastings = hastings - (-log(b_l - a_l));
        hastings = hastings + log(1 - omega) - log(omega);
      }
      else {// Add
        hastings = 0;
        // Old
        hastings = hastings - loglambda_null(j);
        // New
        l_temp = exp(r_truncnorm(log(l_start), 10.0*tau_log_l, log(a_l), log(b_l)));
        K = K_builder(dist, l_temp);
        K = K.submat(which(H(_, j), 0), which(H(_, j), 0));
        log_det(log_K_det_temp, log_K_det_sign_temp, K);
        solve(temp, K, create_ones(n_eff(j), 1), solve_opts::likely_sympd + solve_opts::no_approx);
        K_inv_sum_temp = (create_ones(1, n_eff(j))*temp).eval()(0,0);
        solve(temp, K, logLambda.submat(which(H(_, j), 0), as_uvec(j)), solve_opts::likely_sympd + solve_opts::no_approx);
        G1_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*temp).eval()(0,0);
        G2_temp = (temp.t()*create_ones(n_eff(j), n_eff(j))*temp/(K_inv_sum_temp + 1.0/h)).eval()(0,0);
        hastings = hastings + (-log_K_det_temp*log_K_det_sign_temp/2.0 - log(K_inv_sum_temp + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0));
        logBF(it, j) = hastings;
        // Proposal
        hastings = hastings - (-pow(log(l_temp) - log(l_start), 2)/2.0/pow(10.0*tau_log_l, 2));
        // Prior
        hastings = hastings + (-log(b_l - a_l));
        hastings = hastings - (log(1 - omega) - log(omega));
      }
      if (it > burn) {
        try_gamma++;
      }
      if (hastings >= log(double(rand()%10001)/10000)) {
        gamma(j) = gamma_temp;
        if(gamma_temp == 1) {// Add
          gamma_sum_temp++;
          l(j) = l_temp;
          log_K_det(j) = log_K_det_temp*log_K_det_sign_temp;
          K_inv_sum(j) = K_inv_sum_temp;
          G1(j) = G1_temp;
          G2(j) = G2_temp;
        }
        else {// Delete
          gamma_sum_temp--;
        }
        if(it > burn) {
          accept_gamma++;
        }
      }
    }
    // omega = rbeta(1, a_omega + gamma_sum_temp, b_omega + p - gamma_sum_temp)(0);
    
    // Update kernel parameter l
    for(j = 0; j < p; j++) {
      if(gamma(j) == 1) {
        hastings = 0;
        // Old
        hastings = hastings - (-log_K_det(j)/2.0 - log(K_inv_sum(j) + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1(j)/2.0 - G2(j)/2.0));
        //New
        l_temp = exp(r_truncnorm(log(l(j)), tau_log_l, log(a_l), log(b_l)));
        K = K_builder(dist, l_temp);
        K = K.submat(which(H(_, j), 0), which(H(_, j), 0));
        log_det(log_K_det_temp, log_K_det_sign_temp, K);
        solve(temp, K, create_ones(n_eff(j), 1), solve_opts::likely_sympd + solve_opts::no_approx);
        K_inv_sum_temp = (create_ones(1, n_eff(j))*temp).eval()(0,0);
        solve(temp, K, logLambda.submat(which(H(_, j), 0), as_uvec(j)), solve_opts::likely_sympd + solve_opts::no_approx);
        G1_temp = (logLambda.submat(which(H(_, j), 0), as_uvec(j)).t()*temp).eval()(0,0);
        G2_temp = (temp.t()*create_ones(n_eff(j), n_eff(j))*temp/(K_inv_sum_temp + 1.0/h)).eval()(0,0);
        hastings = hastings + (-log_K_det_temp*log_K_det_sign_temp/2.0 - log(K_inv_sum_temp + 1.0/h)/2.0 - (a_sigma + 1.0*n_eff(j)/2.0)*log(b_sigma + G1_temp/2.0 - G2_temp/2.0));
        if (it > burn) {
          try_l++;
        }
        if (hastings >= log(double(rand()%10001)/10000)) {
          l(j) = l_temp;
          log_K_det(j) = log_K_det_temp*log_K_det_sign_temp;
          K_inv_sum(j) = K_inv_sum_temp;
          G1(j) = G1_temp;
          G2(j) = G2_temp;
          if (it > burn) {
            accept_l++;
          }
        }
      }
      else {
        l(j) = 0;
      }
    }
    
    // Monitor the process
    if(it*100/iter == count) {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    H_sum(it) = 0;
    gamma_sum(it) = gamma_sum_temp;
    for(i = 0; i < n; i++) {
      pi_store(it, i) = pi(i);
      H_sum(it) = H_sum(it) + H_sum_temp(i);
    }
    for(j = 0; j < p; j++) {
      phi_store(it, j) = phi(j);
      gamma_store(it, j) = gamma(j);
      omega_store(it) = omega;
      l_store(it, j) = l(j);
      if(it >= burn) {
        gamma_ppi(j) = gamma_ppi(j) + gamma(j);
        for(i = 0; i < n; i++) {
          H_ppi(i, j) = H_ppi(i, j) + H(i, j);
        }
      }
    }
  }
  for(j = 0; j < p; j++) {
    gamma_ppi(j) = gamma_ppi(j)/(iter - burn);
    for(i = 0; i < n; i++) {
      H_ppi(i, j) = H_ppi(i, j)/(iter - burn);
    }
  }
  accept_phi = accept_phi/try_phi;
  accept_lambda = accept_lambda/try_lambda;
  accept_gamma = accept_gamma/try_gamma;
  accept_l = accept_l/try_l;

  return Rcpp::List::create(Rcpp::Named("logBF") = logBF, Rcpp::Named("H_sum") = H_sum, Rcpp::Named("H_ppi") = H_ppi, Rcpp::Named("pi") = pi_store, Rcpp::Named("phi") = phi_store, Rcpp::Named("logLambda") = logLambda, Rcpp::Named("gamma_sum") = gamma_sum, Rcpp::Named("gamma_ppi") = gamma_ppi, Rcpp::Named("gamma") = gamma_store, Rcpp::Named("omega") = omega_store, Rcpp::Named("l") = l_store, Rcpp::Named("accept_lambda") = accept_lambda, Rcpp::Named("accept_gamma") = accept_gamma, Rcpp::Named("accept_l") = accept_l);
}

arma::mat K_builder(arma::mat dist, double l) {
  int n = dist.n_rows;
  int i;
  arma::mat K(n, n);
  arma::vec eigval;
  arma::mat eigvec;
  K = exp(-dist % dist/2.0/l/l);
  eig_sym(eigval, eigvec, K);
  for (i = 0; i < n; i++) {
    if (eigval(i) < 1e-08) {
      eigval(i) = 1e-08;
    } else {
      break;
    }
  }
  K = eigvec * diagmat(eigval) * inv(eigvec);
  return (K);
}

arma::mat create_ones(int n, int p) {
  arma::mat ones(n, p);
  ones.ones();
  return (ones);
}

arma::mat create_eye(int n) {
  arma::mat I(n, n);
  I.eye();
  return (I);
}

int num(IntegerVector x, int c) {
  int n = x.size();
  int count = 0;
  int i;
  for (i = 0; i < n; i++) {
    if (x(i) == c) {
      count++;
    }
  }
  return count;
}

arma::uvec which(IntegerVector x, int c) {
  int n = x.size();
  int count = 0;
  int i;
  int m = num(x, c);
  arma::uvec index(m);
  for (i = 0; i < n; i++) {
    if (x(i) == c) {
      index(count) = i;
      count++;
    }
  }
  return index;
}

arma::uvec as_uvec(int j) {
  arma::uvec temp(1);
  temp(0) = j;
  return temp;
}

arma::uvec which_stop(IntegerVector x) {
  int n = x.size();
  int count = 0;
  arma::uvec index(n);
  int i;
  for (i = 0; i < n; i++) {
    if (x(i) == -1) {
      break;
    }
    index(count) = x(i);
    count++;
  }
  return index.subvec(0, count - 1);
}
  