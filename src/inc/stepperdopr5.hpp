#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>

#define THROW(x)                                                               \
  {                                                                            \
    puts(x);                                                                   \
    std::exit(EXIT_FAILURE);                                                   \
  }

template <int n, class F> struct StepperDopr5 {
  using VecDoub = std::array<double, n>;
  static constexpr double EPS = std::numeric_limits<double>::epsilon();

  VecDoub dydx, dydxnew;
  VecDoub y, yout, yerr, ytemp;
  VecDoub k2, k3, k4, k5, k6;
  double atol, rtol;
  double x, hnext, errold;
  bool reject;

  StepperDopr5(double atol_, double rtol_) : atol(atol_), rtol(rtol_) {}

  void set_init(VecDoub const &y_init, VecDoub const &dydx_init, double x_init,
                double h_init) {
    y = y_init;
    dydx = dydx_init;
    x = x_init;
    hnext = h_init;
    reject = false;
    errold = 1.0e-4;
  }

  void integrate_to(double x_final, F &&derivs) {
    for (;;) {
      if ((x + hnext * 1.0001 - x_final) * hnext > 0.0) {
        hnext = x_final - x;
      }
      step(derivs);
      if ((x - x_final) * hnext >= 0.0) {
        break;
      }
    }
  }

  void step(F &&derivs) {
    double h = hnext;
    for (;;) {
      dy(h, derivs);
      double err = error();
      if (success(err, h))
        break;
      if (std::abs(h) <= std::abs(x) * EPS) {
        THROW("stepsize underflow in StepperDopr5");
      }
    }
    dydx = dydxnew;
    y = yout;
    x += h;
  }

  void dy(double h, F &&derivs) {
    static constexpr double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0,
                            a21 = 0.2, a31 = 3.0 / 40.0, a32 = 9.0 / 40.0,
                            a41 = 44.0 / 45.0, a42 = -56.0 / 15.0,
                            a43 = 32.0 / 9.0, a51 = 19372.0 / 6561.0,
                            a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0,
                            a54 = -212.0 / 729.0, a61 = 9017.0 / 3168.0,
                            a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0,
                            a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
                            a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0,
                            a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0,
                            a76 = 11.0 / 84.0, e1 = 71.0 / 57600.0,
                            e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
                            e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0,
                            e7 = -1.0 / 40.0;

    int i;
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * a21 * dydx[i];
    derivs(x + c2 * h, ytemp, k2);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
    derivs(x + c3 * h, ytemp, k3);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a41 * dydx[i] + a42 * k2[i] + a43 * k3[i]);
    derivs(x + c4 * h, ytemp, k4);
    for (i = 0; i < n; i++)
      ytemp[i] =
          y[i] + h * (a51 * dydx[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    derivs(x + c5 * h, ytemp, k5);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a61 * dydx[i] + a62 * k2[i] + a63 * k3[i] +
                             a64 * k4[i] + a65 * k5[i]);
    double xph = x + h;
    derivs(xph, ytemp, k6);
    for (i = 0; i < n; i++)
      yout[i] = y[i] + h * (a71 * dydx[i] + a73 * k3[i] + a74 * k4[i] +
                            a75 * k5[i] + a76 * k6[i]);
    derivs(xph, yout, dydxnew);
    for (i = 0; i < n; i++) {
      yerr[i] = h * (e1 * dydx[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] +
                     e6 * k6[i] + e7 * dydxnew[i]);
    }
  }

  double error() {
    double err = 0.0, sk;
    for (int i = 0; i < n; i++) {
      sk = atol + rtol * std::max(std::abs(y[i]), std::abs(yout[i]));
      double yerr2 = yerr[i] / sk;
      err += yerr2 * yerr2;
    }
    return std::sqrt(err / n);
  }

  bool success(double err, double &h) {
    static constexpr double beta = 0.4 / 5.0, alpha = 0.2 - beta * 0.75,
                            safe = 0.9, minscale = 0.2, maxscale = 10.0;
    double scale;
    if (err <= 1.0) {
      if (err == 0.0)
        scale = maxscale;
      else {
        scale = safe * std::pow(err, -alpha) * std::pow(errold, beta);
        if (scale < minscale)
          scale = minscale;
        if (scale > maxscale)
          scale = maxscale;
      }
      if (reject)
        hnext = h * std::min(scale, 1.0);
      else
        hnext = h * scale;
      errold = std::max(err, 1.0e-4);
      reject = false;
      return true;
    } else {
      scale = std::max(safe * std::pow(err, -alpha), minscale);
      h *= scale;
      reject = true;
#if 0
      printf("rejecting step, err = %g, scale = %g, h = %g, x=%g\n", err, scale,
             h, x);
#endif
      return false;
    }
  }
};

#undef THROW
