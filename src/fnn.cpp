#include <Rcpp.h>
#include <random>
#include <math.h>

// a function for a dual fnoise run
// mean 0.0 sd 1.0
RcppExport SEXP fnn(SEXP num, SEXP time, SEXP gamma, SEXP korrelaatio)
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  const double pii = 3.14159265359;
  int n, t, ss;
  double g, roo, spec, vaihe, freq;
  double ka, khaj;
  Rcpp::NumericVector numv(num);
  Rcpp::NumericVector timev(time);
  Rcpp::NumericVector gammav(gamma);
  Rcpp::NumericVector roov(korrelaatio);
  n = numv[0];
  t = timev[0];
  g = gammav[0];
  roo = roov[0];
  ss = t/2;
  double** w = new double*[n-1];
  double** x = new double*[n-1];
  for(int i = 0; i < n; ++i){
    w[i] = new double[t];
    x[i] = new double[t];
  }
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < t; ++j){
      w[i][j] = 0.0;
      x[i][j] = j*(2.0*pii/((double)t - 1.0));
    }
  }
  // phase generation here
  for(int i = 0; i < ss; ++i){
    freq = 1.0 + (double)i;
    spec = 1.0 / (pow(freq,g));
    for(int j = 0; j < n; ++j){
      vaihe = dist(mt) * 2.0 * pii + (1.0-roo) * dist(mt) * 2.0 * pii;
      for(int k = 0; k < t; ++k){
	w[j][k] = w[j][k] + spec*sin(x[j][k]*freq - vaihe);
      }
    }
  }
  Rcpp::NumericMatrix ulos(t,n);
  for(int j = 0; j < n; ++j){
    ka = 0.0;
    khaj = 0.0;
    for(int i = 0; i < t; ++i){
      ka += w[j][i];
    }
    ka = ka / t;
    for(int i = 0; i < t; ++i){
      khaj += (w[j][i]-ka)*(w[j][i]-ka);
    }
    khaj = sqrt(khaj/((double)t - 1.0));
    for(int i = 0; i < t; ++i){
      ulos(i,j) = (w[j][i]-ka)/khaj;
    }
  }
  for(int i = 0; i < n; ++i){
    delete [] w[i];
    delete [] x[i];
  }
  delete [] w;
  delete [] x;
  return(ulos);
}
