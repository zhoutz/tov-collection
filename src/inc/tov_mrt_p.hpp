#pragma once

#include "constants.hpp"
#include "eos.hpp"
#include <cmath>

using d3 = std::array<double, 3>;

extern unsigned g_cnt;

struct TOV_mrt_p {
  EOS &eos;

  void operator()(double p, d3 const &wmy, d3 &dwmy_dp) {
    ++g_cnt;
    auto [e, cs2] = eos.e_cs2_of_p(p);
    double kappa = 1. / cs2;
    tov_impl(wmy, dwmy_dp, p, e, kappa);
  }

  void tov_impl(d3 const &wmy, d3 &dwmy_dp, double p, double e, double kappa) {
    double w = wmy[0];
    double m = wmy[1];
    double y = wmy[2];
    double r = std::cbrt(w);

    double dwdp = -3. * w * (C1 * r - 2. * m) / ((e + p) * (m + C2 * w * p));
    double dmdp = C2 * e / 3. * dwdp;
    double dydp = -(2. * m * y * (y + 2.) - C1 * r * (y * y + y - 6.) +
                    C2 * w * (y * (e + 3. * p) - (e + p) * (3. + kappa))) /
                  ((e + p) * (m + C2 * w * p));

    dwmy_dp[0] = dwdp;
    dwmy_dp[1] = dmdp;
    dwmy_dp[2] = dydp;
  }
};
