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

    double C2r2 = C2 * r * r;
    double C1r_2m = C1 * r - 2. * m;

    double dpdr_part1 = e + p;
    double dpdr_part2 = m + C2r2 * r * p;
    double dpdr_part3 = r * C1r_2m;

    double dpdr = -dpdr_part1 * dpdr_part2 / dpdr_part3;
    double dmdr = C2r2 * e;

    double del = y - 2.;
    double eta = 2. * m / C1r_2m;

    double dydr_part1 = (-del * (del + 5.) + 2. * eta * (1. - del)) / r;
    double dydr_part2 =
        6. * m * y / r + C2r2 * (y * (e + 3. * p) - (e + p) * (3. + kappa));

    double dydr = dydr_part1 + dydr_part2 / C1r_2m;

    dpmy_dr[0] = dpdr;
    dpmy_dr[1] = dmdr;
    dpmy_dr[2] = dydr;
  }
};
