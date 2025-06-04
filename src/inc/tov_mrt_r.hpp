#pragma once

#include "constants.hpp"
#include "eos.hpp"

using d3 = std::array<double, 3>;

extern unsigned g_cnt;

struct TOV_mrt_r {
  EOS_epcs &eos;

  void operator()(double r, d3 const &pmy, d3 &dpmy_dr) {
    ++g_cnt;
    double p = pmy[0];
    auto [e, cs2] = eos.e_cs2_of_p(p);
    double kappa = 1. / cs2;
    tov_impl(pmy, dpmy_dr, r, e, kappa);
  }

  void tov_impl(d3 const &pmy, d3 &dpmy_dr, double r, double e, double kappa) {
    double p = pmy[0];
    double m = pmy[1];
    double y = pmy[2];

    double dpdr = -(e + p) * (m + C2 * r * r * r * p) / (r * (C1 * r - 2. * m));
    double dmdr = C2 * r * r * e;
    double dydr =
        (2. * m * y * (y + 2.) - C1 * r * (y * y + y - 6.) +
         C2 * r * r * r * (y * (e + 3. * p) - (e + p) * (3. + kappa))) /
        (r * (C1 * r - 2. * m));

    dpmy_dr[0] = dpdr;
    dpmy_dr[1] = dmdr;
    dpmy_dr[2] = dydr;
  }
};
