#include <Rcpp.h>
using namespace Rcpp;

// ---------- helpers ----------
static inline double log1pexp_stable(double x) {
  if (x > 35.0) return x;
  if (x < -35.0) return std::exp(x);
  return std::log1p(std::exp(x));
}

static inline double logSumExp_ptr(const double* x, int n) {
  double m = x[0];
  for (int i = 1; i < n; ++i) if (x[i] > m) m = x[i];
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += std::exp(x[i] - m);
  return m + std::log(s);
}

// family_id: 1 binomial, 2 poisson, 3 Gamma, 4 gaussian
static inline double clink_canonical(int family_id, double mu) {
  if (family_id == 1) { // binomial
    if (mu <= 0.0) mu = 1e-15;
    if (mu >= 1.0) mu = 1.0 - 1e-15;
    return std::log(mu / (1.0 - mu));
  } else if (family_id == 2) { // poisson
    if (mu <= 0.0) mu = 1e-15;
    return std::log(mu);
  } else if (family_id == 3) { // Gamma
    if (mu <= 0.0) mu = 1e-15;
    return -1.0 / mu;
  } else { // gaussian
    return mu;
  }
}

static inline double bfun_canonical(int family_id, double theta) {
  if (family_id == 1) { // binomial
    return log1pexp_stable(theta);
  } else if (family_id == 2) { // poisson
    return std::exp(theta);
  } else if (family_id == 3) { // Gamma
    if (theta >= -1e-15) theta = -1e-15;
    return -std::log(-theta);
  } else { // gaussian
    return 0.5 * theta * theta;
  }
}

// link_id (make.link): 1 logit, 2 probit, 3 cauchit, 4 cloglog, 5 identity, 6 log, 7 sqrt, 8 1/mu^2, 9 inverse
static inline double linkinv_std(int link_id, double eta) {
  if (link_id == 1) { // logit
    if (eta >= 0.0) {
      double e = std::exp(-eta);
      return 1.0 / (1.0 + e);
    } else {
      double e = std::exp(eta);
      return e / (1.0 + e);
    }
  } else if (link_id == 2) { // probit
    return R::pnorm5(eta, 0.0, 1.0, 1, 0);
  } else if (link_id == 3) { // cauchit
    return 0.5 + std::atan(eta) / M_PI;
  } else if (link_id == 4) { // cloglog
    double e = std::exp(eta);
    return 1.0 - std::exp(-e);
  } else if (link_id == 5) { // identity
    return eta;
  } else if (link_id == 6) { // log
    return std::exp(eta);
  } else if (link_id == 7) { // sqrt
    return eta * eta;
  } else if (link_id == 8) { // 1/mu^2
    if (eta <= 1e-15) eta = 1e-15;
    return 1.0 / std::sqrt(eta);
  } else { // inverse
    if (std::fabs(eta) <= 1e-15) eta = (eta >= 0 ? 1e-15 : -1e-15);
    return 1.0 / eta;
  }
}

// [[Rcpp::export]]
Rcpp::List ghq_nll_groups_fast(
    const NumericVector& etak,
    const NumericVector& y,
    const NumericVector& factor,
    const NumericVector& Z,
    const IntegerVector& gptr,       // length G+1, 0-based indices in arrays
    const double tau2,
    const int family_id,
    const int link_id,
    const NumericVector& ghq_x,       // length Q
    const NumericVector& ghq_w,       // length Q (for exp(-x^2) rule)
    Rcpp::Function linkinv_fun        // fallback only
) {
  const int n = etak.size();
  if (y.size() != n || factor.size() != n || Z.size() != n) {
    stop("etak/y/factor/Z must have the same length.");
  }
  if (tau2 <= 0.0) stop("tau2 must be positive inside ghq_nll_groups_fast.");

  const int G = gptr.size() - 1;
  if (G <= 0) stop("gptr must have length >= 2.");

  const int Q = ghq_x.size();
  if (ghq_w.size() != Q) stop("ghq_x and ghq_w must have same length.");
  if (Q <= 0) stop("Q must be positive.");

  NumericVector nll_g(G);
  const double s = std::sqrt(2.0 * tau2);
  const double log_sqrt_pi = 0.5 * std::log(M_PI);

  // work buffer for zq (log-domain)
  std::vector<double> zq(Q);

  for (int g = 0; g < G; ++g) {
    const int start = gptr[g];
    const int end   = gptr[g+1];
    if (start < 0 || end < start || end > n) stop("invalid gptr.");

    // initialize with log weights
    for (int q = 0; q < Q; ++q) {
      double w = ghq_w[q];
      if (w <= 0.0) w = 1e-300;
      zq[q] = std::log(w);
    }

    for (int j = start; j < end; ++j) {
      const double et0 = etak[j];
      const double fj  = factor[j];
      const double yj  = y[j];
      const double zj  = Z[j];

      for (int q = 0; q < Q; ++q) {
        const double vq  = s * ghq_x[q];
        const double eta = et0 + vq * zj;

        double mu;
        if (link_id >= 1 && link_id <= 9) {
          mu = linkinv_std(link_id, eta);
        } else {
          mu = as<double>(linkinv_fun(eta));
        }

        const double theta = clink_canonical(family_id, mu);
        const double bj    = bfun_canonical(family_id, theta);
        zq[q] += fj * ( yj * theta - bj );
      }
    }

    const double log_sum = logSumExp_ptr(zq.data(), Q);
    const double logLg   = log_sum - log_sqrt_pi; // 1/sqrt(pi) factor
    nll_g[g] = -logLg;
  }

  double total = 0.0;
  for (int g = 0; g < G; ++g) total += nll_g[g];

  return List::create(
    _["nll_total"] = total,
    _["nll_group"] = nll_g
  );
}
