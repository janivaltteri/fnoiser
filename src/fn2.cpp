#include <Rcpp.h>
#include <random>
#include <math.h>

// a function for a dual fnoise run
// mean 0.0 sd 1.0
RcppExport SEXP fn2(SEXP time, SEXP gamma, SEXP korrelaatio)
{
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  const double pii = 3.14159265359;
  int t, ss;
  double g, roo, spec, vaihe1, vaihe2, freq;
  double ka1, ka2, khaj1, khaj2;
  Rcpp::NumericVector timev(time);
  Rcpp::NumericVector gammav(gamma);
  Rcpp::NumericVector roov(korrelaatio);
  t = timev[0];
  g = gammav[0];
  roo = roov[0];
  ss = t/2;
  double** w = new double*[2];
  double** x = new double*[2];
  for(int i = 0; i < 2; ++i){
    w[i] = new double[t];
    x[i] = new double[t];
  }
  for(int i = 0; i < 2; ++i){
    for(int j = 0; j < t; ++j){
      w[i][j] = 0.0;
      x[i][j] = j*(2.0*pii/((double)t - 1.0));
    }
  }
  for(int i = 0; i < ss; ++i){
    vaihe1 = dist(mt) * 2.0 * pii;
    vaihe2 = acos(roo) + vaihe1;
    freq = 1.0 + (double)i;
    spec = 1.0 / (pow(freq,g));
    for(int k = 1; k < t; ++k){
      w[0][k] = w[0][k] + spec*sin(x[0][k]*freq-vaihe1);
      w[1][k] = w[1][k] + spec*sin(x[1][k]*freq-vaihe2);
    }
  }
  ka1 = 0.0;
  khaj1 = 0.0;
  ka2 = 0.0;
  khaj2 = 0.0;
  for(int i = 0; i < t; ++i){
    ka1 += w[0][i];
    ka2 += w[1][i];
  }
  ka1 = ka1 / t;
  ka2 = ka2 / t;
  for(int i = 0; i < t; ++i){
    khaj1 += (w[0][i]-ka1)*(w[0][i]-ka1);
    khaj2 += (w[1][i]-ka2)*(w[1][i]-ka2);
  }
  khaj1 = sqrt(khaj1/((double)t - 1.0));
  khaj2 = sqrt(khaj2/((double)t - 1.0));
  Rcpp::NumericMatrix ulos(t,2);
  for(int i = 0; i < t; ++i){
    ulos(i,0) = (w[0][i]-ka1)/khaj1;
    ulos(i,1) = (w[1][i]-ka2)/khaj2;
  }
  for(int i = 0; i < 2; ++i){
    delete [] w[i];
    delete [] x[i];
  }
  delete [] w;
  delete [] x;
  return(ulos);
}
