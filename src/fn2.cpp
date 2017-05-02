#include <Rcpp.h>
#include <random>
#include <vector>
#include <math.h>
#include "boost/multi_array.hpp"

typedef boost::multi_array<double,2> matriisi;
typedef matriisi::index indeksi;

// a function for a single fnoise run
// mean 0.0 sd 1.0
RcppExport SEXP fn2(SEXP time, SEXP gamma, SEXP korrelaatio)
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(1.0, 10.0);
  const double pii = 3.14159265359;
  double t, g, roo, ss, spec, vaihe, freq;
  Rcpp::NumericVector timev(time);
  Rcpp::NumericVector gammav(gamma);
  Rcpp::NumericVector roov(korrelaatio);
  t = timev[0];
  g = gammav[0];
  roo = roov[0];
  ss = t/2;
  std::vector<double> w, x;
  w.resize(t);
  x.resize(t);
  for(int i = 0; i < t; ++i){
    w[i] = 0.0;
    x[i] = i*(2.0*pii/((double)t - 1.0));
  }
  for(int i = 0; i < ss; ++i){
    vaihe = dist(mt) * 2.0 * pii;
    freq = 1.0 + (double)i;
    spec = 1.0 / (pow(freq,g));
    for(int k = 1; k < t; ++k){
      w[k] = w[k] + spec*sin(x[k]*freq-vaihe);
    }
  }
  double ka = 0.0;
  double khaj = 0.0;
  for(int i = 0; i < t; ++i){
    ka += w[i];
  }
  ka = ka / t;
  for(int i = 0; i < t; ++i){
    khaj += (w[i]-ka)*(w[i]-ka);
  }
  khaj = sqrt(khaj/((double)t - 1.0));
  Rcpp::NumericVector ulos(t);
  for(int i = 0; i < t; ++i){
    ulos[i] = (w[i]-ka)/khaj;
  }
  return(ulos);
}
