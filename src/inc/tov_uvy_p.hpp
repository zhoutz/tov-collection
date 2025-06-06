#pragma once

#include "constants.hpp"
#include "eos.hpp"

using d3 = std::array<double, 3>;

extern unsigned g_cnt;

struct TOV_uvy_p {
  EOS &eos;

  void operator()(double p, d3 const &uvy, d3 &duvy_dp) {
    ++g_cnt;
    auto [e, cs2] = eos.e_cs2_of_p(p);
    double kappa = 1. / cs2;
    tov_impl(uvy, duvy_dp, p, e, kappa);
  }

  void tov_impl(d3 const &uvy, d3 &duvy_dp, double p, double e, double kappa) {
    double u = uvy[0];
    double v = uvy[1];
    double y = uvy[2];

    double deno = 1. / ((e + p) * (v + C2 * u * p));
    double dudp = -2. * u * (C1 - 2. * v) * deno;
    double dvdp = (C1 - 2. * v) * (v - C2 * u * e) * deno;
    double dydp = (C2 * u * ((3. + kappa) * (e + p) - y * (e + 3. * p)) -
                   2. * v * y * (2. + y) + C1 * (y * y + y - 6.)) *
                  deno;

    duvy_dp[0] = dudp;
    duvy_dp[1] = dvdp;
    duvy_dp[2] = dydp;
  }
};
