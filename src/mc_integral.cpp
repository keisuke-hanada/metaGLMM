#include <Rcpp.h>
using namespace Rcpp;

// ---------- utilities ----------
inline double log1pexp(double x){
  if (x > 0) return x + std::log1p(std::exp(-x));
  return std::log1p(std::exp(x));
}

inline double inv_logit(double x){
  if (x >= 0) {
    double e = std::exp(-x);
    return 1.0/(1.0+e);
  } else {
    double e = std::exp(x);
    return e/(1.0+e);
  }
}

// linkinv: make.link compatible
inline double linkinv(double eta, int link_code){
  switch(link_code){
  case 1: return inv_logit(eta);                        // logit
  case 2: return R::pnorm5(eta, 0.0, 1.0, 1, 0);         // probit
  case 3: return 0.5 + std::atan(eta)/M_PI;             // cauchit
  case 4: {                                             // cloglog
    // 1 - exp(-exp(eta))
    if (eta > 30) return 1.0;
    double t = std::exp(eta);
    return 1.0 - std::exp(-t);
  }
  case 5: return eta;                                   // identity
  case 6: return std::exp(eta);                         // log
  case 7: return eta * eta;                             // sqrt
  case 8: return 1.0 / std::sqrt(eta);                  // 1/mu^2
  case 9: return 1.0 / eta;                             // inverse
  default: Rcpp::stop("Unsupported link_code");
  }
}

// log-mean-exp for one group vector z[0..B-1]
inline double log_mean_exp_row(const NumericMatrix& z, int g){
  int B = z.ncol();
  double m = z(g, 0);
  for (int b=1; b<B; ++b) if (z(g,b) > m) m = z(g,b);
  double s = 0.0;
  for (int b=0; b<B; ++b) s += std::exp(z(g,b) - m);
  return m + std::log(s) - std::log((double)B);
}

// [[Rcpp::export]]
List mc_nll_groups(
    const NumericVector& etak,
    const NumericVector& y,
    const NumericVector& factor,
    const IntegerVector& gid,     // 0..G-1
    const NumericVector& Z,       // 0/1 or 1
    const NumericVector& rstdnorm, // length B
    double tau2,
    double tau2_min,
    int fam_code,                 // 1 binom, 2 pois, 3 Gamma, 4 gauss
    int link_code                 // make.link code
){
  int N = etak.size();
  int B = rstdnorm.size();
  if (y.size() != N || factor.size() != N || gid.size() != N || Z.size() != N)
    Rcpp::stop("Input length mismatch.");

  int G = 0;
  for (int i=0; i<N; ++i) if (gid[i] + 1 > G) G = gid[i] + 1;
  if (G <= 0) Rcpp::stop("No groups detected.");

  double tau2_use = (tau2 > tau2_min ? tau2 : tau2_min);
  double sd = std::sqrt(tau2_use);

  NumericVector rst(B);
  for (int b=0; b<B; ++b) rst[b] = rstdnorm[b] * sd;

  // z(g,b): accumulate contributions for group g at MC point b
  NumericMatrix z(G, B); // initialized to 0

  const double eps = 1e-12;

  // One pass over observations
  for (int j=0; j<N; ++j){
    int g = gid[j];
    double etaj = etak[j];
    double yj   = y[j];
    double fj   = factor[j];
    double Zj   = Z[j];

    for (int b=0; b<B; ++b){
      double eta = etaj + Zj * rst[b];
      double mu  = linkinv(eta, link_code);

      double theta, btheta;

      if (fam_code == 1){                 // binomial
        mu = std::min(std::max(mu, eps), 1.0 - eps);
        theta  = std::log(mu/(1.0 - mu));
        btheta = log1pexp(theta);
      } else if (fam_code == 2){          // poisson
        mu = std::max(mu, eps);
        theta  = std::log(mu);
        btheta = std::exp(theta);
      } else if (fam_code == 3){          // Gamma (canonical)
        mu = std::max(mu, eps);
        theta  = -1.0/mu;
        btheta = -std::log(-theta);
      } else if (fam_code == 4){          // gaussian
        theta  = mu;
        btheta = 0.5 * theta * theta;
      } else {
        Rcpp::stop("Unsupported fam_code");
      }

      z(g, b) += fj * (yj * theta - btheta);
    }
  }

  NumericVector ll_vec(G);
  double total = 0.0;

  for (int g=0; g<G; ++g){
    double lme = log_mean_exp_row(z, g);
    double ll_g = -lme;     // minus log marginal
    ll_vec[g] = ll_g;
    total += ll_g;
  }

  return List::create(_["total"] = total, _["ll_vec"] = ll_vec);
}
