#pragma once

#include "constants.hpp"
#include "eos.hpp"

using d3 = std::array<double, 3>;

extern unsigned g_cnt;

struct TOV_mrt_h {
  EOS &eos;

  void operator()(double h, d3 const &uvy, d3 &duvy_dh) {
    ++g_cnt;
    auto [e, p, cs2] = eos.e_p_cs2_of_h(h);
    double kappa = 1. / cs2;
    tov_impl(uvy, duvy_dh, p, e, kappa);
  }

  void tov_impl(d3 const &uvy, d3 &duvy_dh, double p, double e, double kappa) {
    double u = uvy[0];
    double v = uvy[1];
    double y = uvy[2];

    double dudh = -2. * u * (C1 - 2. * v) / (v + C2 * u * p);
    double dvdh = -(C1 - 2. * v) * (-v + C2 * u * e) / (v + C2 * u * p);
    double dydh = (C2 * u * ((3. + kappa) * (e + p) - y * (e + 3. * p)) -
                   2. * v * y * (2. + y) + C1 * (y * y + y - 6.)) /
                  (v + C2 * u * p);

    duvy_dh[0] = dudh;
    duvy_dh[1] = dvdh;
    duvy_dh[2] = dydh;
  }
};
