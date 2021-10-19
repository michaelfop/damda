/*
 * Author: Michael Fop
 * Inductive conditional estimation for Gaussian multivariate distribution
 */

#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

/* Inductive conditional estimation */
Rcpp::List ice(arma::mat Y,
               arma::mat z,
               arma::mat mu,
               arma::cube sigma,
               arma::mat ybar,
               arma::cube O,
               arma::uvec obs, arma::uvec ext,
               arma::vec nc, int C)
{
  int Q = ext.size();
  int P = obs.size();
  int N = Y.n_rows;
  arma::cube sigmaout(Q,Q,C);
  arma::cube cross(P,Q,C);
  arma::mat muout(Q,C);
  arma::cube mucond(N,Q,C);
  arma::mat Yobs = Y.cols(obs);
  arma::uvec cvec = linspace<uvec>(0, C-1, C);
  arma::cube E(Q,Q,C);
  
  bool check;
  
  for ( int k = 0; k < C; k++ ) {
    arma::mat ss = sigma.slice(k);
    arma::mat W = O.slice(k);
    // arma::mat inv = inv_sympd( ss(obs,obs) );
    arma::mat inv = solve( ss(obs,obs), eye(P,P) );
    
    arma::mat A = inv.t()*W(obs,obs);
    arma::mat B = W(obs,ext).t()*inv;
    // arma::mat CC = inv_sympd( A*inv.t() ) * B.t();  // this sometimes produces inv_sympd not symmetric warning, not sure why - Solved? See below
    arma::mat CC = solve( A*inv.t(), eye(P,P) ) * B.t();  // this sometimes produces inv_sympd not symmetric warning, not sure why - Solved? See below
    arma::mat D = CC.t()*inv;
    E.slice(k) = 1/nc(k) * ( D*W(obs,obs)*D.t() - 2*(B*CC) + W(ext,ext) );
    E.slice(k) = (E.slice(k) + E.slice(k).t())/2;  // make sure obtained matrix is actually symmetric (not only numerically)
    cross.slice(k) = CC;
    sigmaout.slice(k) = E.slice(k) + D*CC;
    sigmaout.slice(k) = (sigmaout.slice(k) + sigmaout.slice(k).t())/2;
    
    arma::uvec kappa = find( cvec == k );
    arma::mat Ym = Yobs.each_row() - mu(obs,kappa).t();
    arma::mat Yz = Ym;
    Yz.each_col() %= z.col(k);
    
    mucond.slice(k) = Ym*D.t();
    muout.col(k) = ybar(ext,kappa) - sum(Yz*D.t(), 0).t()*1/nc(k);
    mucond.slice(k).each_row() += muout.col(k).t();
  }
  
  return Rcpp::List::create(Named("sigma") = sigmaout,
                            Named("cross") = cross,
                            Named("mu") = muout,
                            Named("sigmacond") = E,
                            Named("mucond") = mucond);
}