#pragma once

#include "constants.hpp"
#include "eos.hpp"
#include <cmath>

using d3 = std::array<double, 3>;

extern unsigned g_cnt;

struct TOV_mrt_xi {
  EOS &eos;

  void operator()(double xi, d3 const &wmy, d3 &dwmy_dxi) {
    ++g_cnt;
    double p = std::exp(xi);
    auto [e, cs2] = eos.e_cs2_of_p(p);
    double kappa = 1. / cs2;
    tov_impl(wmy, dwmy_dxi, p, e, kappa);
  }

  void tov_impl(d3 const &wmy, d3 &dwmy_dxi, double p, double e, double kappa) {
    double w = wmy[0];
    double m = wmy[1];
    double y = wmy[2];
    double r = std::cbrt(w);

    double deno = 1. / ((e + p) * (m + C2 * w * p));
    double dwdxi = -3. * w * p * (C1 * r - 2. * m) * deno;
    double dmdxi = C2 * e / 3. * dwdxi;
    double dydxi = -p *
                   (2. * m * y * (y + 2.) - C1 * r * (y * y + y - 6.) +
                    C2 * w * (y * (e + 3. * p) - (e + p) * (3. + kappa))) *
                   deno;

    dwmy_dxi[0] = dwdxi;
    dwmy_dxi[1] = dmdxi;
    dwmy_dxi[2] = dydxi;
  }
};
