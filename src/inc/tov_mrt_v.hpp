#pragma once

#include "constants.hpp"
#include "eos.hpp"
#include <cmath>

using d3 = std::array<double, 3>;

extern unsigned g_cnt;

struct TOV_mrt_v {
  EOS_epcs &eos;

  void operator()(double p, d3 const &vmy, d3 &dvmy_dp) {
    ++g_cnt;
    auto [e, cs2] = eos.e_cs2_of_p(p);
    double kappa = 1. / cs2;
    tov_impl(vmy, dvmy_dp, p, e, kappa);
  }

  void tov_impl(d3 const &vmy, d3 &dvmy_dp, double p, double e, double kappa) {
    double v = vmy[0];
    double m = vmy[1];
    double y = vmy[2];

    double r = std::cbrt(v);

    double C1r_2m = C1 * r - 2. * m;

    double dvdp_part1 = 3. * v * C1r_2m;
    double dvdp_part2 = e + p;
    double dvdp_part3 = m + C2 * v * p;

    double dmdv = C2 / 3. * e;

    double dvdp = -dvdp_part1 / (dvdp_part2 * dvdp_part3);
    double dmdp = dmdv * dvdp;

    double del = y - 2.;
    double eta = 2. * m / C1r_2m;

    double dydv_part1 = (-del * (del + 5.) + 2. * eta * (1. - del)) / (3. * v);
    double dydv_part2 =
        2. * m * y / v + C2 / 3. * (y * (e + 3. * p) - (e + p) * (3. + kappa));

    double dydv = dydv_part1 + dydv_part2 / C1r_2m;
    double dydp = dydv * dvdp;

    dvmy_dp[0] = dvdp;
    dvmy_dp[1] = dmdp;
    dvmy_dp[2] = dydp;
  }
};
