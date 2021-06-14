#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
int sim_null_realisation(int m, vec P) {
  int M = P.size();
  int ind = m;
  
  vec S = regspace(1, M);
  vec R(M, fill::zeros);
  R.subvec(0, m-1) = Rcpp::RcppArmadillo::sample(S, m, true, P);
  vec Q = unique(R); // unique elements drawn (plus 0)
  int qsize = Q.size() - 1; // number of unique elements drawn
  
  while(qsize < m) {
    R.subvec(ind, ind + m - qsize - 1) = Rcpp::RcppArmadillo::sample(S, m - qsize, true, P);
    Q = unique(R);
    ind = ind + m - qsize;
    qsize = Q.size() - 1;
  }
  
  return ind - m;
}

int mod(int a, int n)
{
  return a - floor(a/n)*n;
}   

// [[Rcpp::export]]
vec sim_null_dist(int m, int n, vec P) {
  vec S(n);
  for(int i = 0; i<n; i++) {
    S[i] = sim_null_realisation(m, P);
    if(mod(i, 1000) == 0) {
      if(mod(i,100000)==0) {
        Rcout << "#\n";
      }
      else{
        Rcout << "#";
      }
    }
  }
  return S;
}

