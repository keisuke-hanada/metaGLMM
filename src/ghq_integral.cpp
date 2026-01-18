#include <Rcpp.h>
using namespace Rcpp;

// ---------- helpers ----------
static inline double log1pexp_stable(double x) {
  // stable log(1+exp(x))
  if (x > 35.0) return x;
  if (x < -35.0) return std::exp(x);
  return std::log1p(std::exp(x));
}

static inline double logSumExp_vec(const NumericVector& x) {
  double m = x[0];
  for (int i = 1; i < x.size(); ++i) if (x[i] > m) m = x[i];
  double s = 0.0;
  for (int i = 0; i < x.size(); ++i) s += std::exp(x[i] - m);
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

// ---------- GHQ nodes/weights (Q=30) for ∫ f(x) e^{-x^2} dx ----------
static const int GHQ_Q = 30;

static const double GHQ_X[GHQ_Q] = {
  -6.8633452935298916, -6.1382792201239340, -5.5331471515674957, -4.9889189685899439,
  -4.4830553570925183, -4.0039086038612288, -3.5444438731553499, -3.0999705295864417,
  -2.6671321245356172, -2.2433914677615041, -1.8267411436036880, -1.4155278001981885,
  -1.0083382710467235, -0.6039210586255523, -0.2011285765488715,
  0.2011285765488715,  0.6039210586255523,  1.0083382710467235,  1.4155278001981885,
  1.8267411436036880,  2.2433914677615041,  2.6671321245356172,  3.0999705295864417,
  3.5444438731553499,  4.0039086038612288,  4.4830553570925183,  4.9889189685899439,
  5.5331471515674957,  6.1382792201239340,  6.8633452935298916
};

static const double GHQ_W[GHQ_Q] = {
  0.0000000000000076, 0.0000000000013436, 0.0000000000616650, 0.0000000010004360,
  0.0000000772250830, 0.0000035191638240, 0.0000960143941460, 0.0016423195080000,
  0.0183404509859080, 0.1317306544729000, 0.6145224490080000, 1.9199528384000000,
  4.0863623778000000, 5.9221635220000000, 5.9221635220000000,
  5.9221635220000000, 5.9221635220000000, 4.0863623778000000, 1.9199528384000000,
  0.6145224490080000, 0.1317306544729000, 0.0183404509859080, 0.0016423195080000,
  0.0000960143941460, 0.0000035191638240, 0.0000000772250830, 0.0000000010004360,
  0.0000000000616650, 0.0000000000013436, 0.0000000000000076
};

// [[Rcpp::export]]
Rcpp::List ghq_nll_groups(const NumericVector& etak,
                          const NumericVector& y,
                          const NumericVector& factor,
                          const IntegerVector& gid,
                          const NumericVector& Z,
                          const double tau2,
                          const int family_id,
                          const int link_id,
                          Rcpp::Function linkinv_fun) {

  const int n = etak.size();
  if (y.size() != n || factor.size() != n || gid.size() != n || Z.size() != n) {
    stop("Input vectors must have the same length.");
  }
  if (tau2 <= 0.0) stop("tau2 must be positive inside ghq_nll_groups.");

  int G = 0;
  for (int i = 0; i < n; ++i) if (gid[i] + 1 > G) G = gid[i] + 1;

  NumericVector nll_g(G);
  const double s = std::sqrt(2.0 * tau2);
  const double log_sqrt_pi = 0.5 * std::log(M_PI);

  for (int g = 0; g < G; ++g) {
    NumericVector zq(GHQ_Q);
    for (int q = 0; q < GHQ_Q; ++q) zq[q] = 0.0;

    for (int j = 0; j < n; ++j) {
      if (gid[j] != g) continue;

      const double et0 = etak[j];
      const double fj  = factor[j];
      const double yj  = y[j];
      const double zj  = Z[j];

      for (int q = 0; q < GHQ_Q; ++q) {
        const double vq  = s * GHQ_X[q];
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

    for (int q = 0; q < GHQ_Q; ++q) {
      double w = GHQ_W[q];
      if (w <= 0.0) w = 1e-300;
      zq[q] += std::log(w);
    }

    const double log_sum = logSumExp_vec(zq);
    const double logLg = log_sum - log_sqrt_pi;
    nll_g[g] = -logLg;
  }

  double total = 0.0;
  for (int g = 0; g < G; ++g) total += nll_g[g];

  return List::create(
    _["nll_total"] = total,
    _["nll_group"] = nll_g
  );
}
