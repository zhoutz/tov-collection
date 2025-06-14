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

struct Dopr853_constants {
  static constexpr double c2 = 0.526001519587677318785587544488e-01;
  static constexpr double c3 = 0.789002279381515978178381316732e-01;
  static constexpr double c4 = 0.118350341907227396726757197510e+00;
  static constexpr double c5 = 0.281649658092772603273242802490e+00;
  static constexpr double c6 = 0.333333333333333333333333333333e+00;
  static constexpr double c7 = 0.25e+00;
  static constexpr double c8 = 0.307692307692307692307692307692e+00;
  static constexpr double c9 = 0.651282051282051282051282051282e+00;
  static constexpr double c10 = 0.6e+00;
  static constexpr double c11 = 0.857142857142857142857142857142e+00;
  static constexpr double c14 = 0.1e+00;
  static constexpr double c15 = 0.2e+00;
  static constexpr double c16 = 0.777777777777777777777777777778e+00;

  static constexpr double b1 = 5.42937341165687622380535766363e-2;
  static constexpr double b6 = 4.45031289275240888144113950566e0;
  static constexpr double b7 = 1.89151789931450038304281599044e0;
  static constexpr double b8 = -5.8012039600105847814672114227e0;
  static constexpr double b9 = 3.1116436695781989440891606237e-1;
  static constexpr double b10 = -1.52160949662516078556178806805e-1;
  static constexpr double b11 = 2.01365400804030348374776537501e-1;
  static constexpr double b12 = 4.47106157277725905176885569043e-2;

  static constexpr double bhh1 = 0.244094488188976377952755905512e+00;
  static constexpr double bhh2 = 0.733846688281611857341361741547e+00;
  static constexpr double bhh3 = 0.220588235294117647058823529412e-01;

  static constexpr double er1 = 0.1312004499419488073250102996e-01;
  static constexpr double er6 = -0.1225156446376204440720569753e+01;
  static constexpr double er7 = -0.4957589496572501915214079952e+00;
  static constexpr double er8 = 0.1664377182454986536961530415e+01;
  static constexpr double er9 = -0.3503288487499736816886487290e+00;
  static constexpr double er10 = 0.3341791187130174790297318841e+00;
  static constexpr double er11 = 0.8192320648511571246570742613e-01;
  static constexpr double er12 = -0.2235530786388629525884427845e-01;

  static constexpr double a21 = 5.26001519587677318785587544488e-2;
  static constexpr double a31 = 1.97250569845378994544595329183e-2;
  static constexpr double a32 = 5.91751709536136983633785987549e-2;
  static constexpr double a41 = 2.95875854768068491816892993775e-2;
  static constexpr double a43 = 8.87627564304205475450678981324e-2;
  static constexpr double a51 = 2.41365134159266685502369798665e-1;
  static constexpr double a53 = -8.84549479328286085344864962717e-1;
  static constexpr double a54 = 9.24834003261792003115737966543e-1;
  static constexpr double a61 = 3.7037037037037037037037037037e-2;
  static constexpr double a64 = 1.70828608729473871279604482173e-1;
  static constexpr double a65 = 1.25467687566822425016691814123e-1;
  static constexpr double a71 = 3.7109375e-2;
  static constexpr double a74 = 1.70252211019544039314978060272e-1;
  static constexpr double a75 = 6.02165389804559606850219397283e-2;
  static constexpr double a76 = -1.7578125e-2;

  static constexpr double a81 = 3.70920001185047927108779319836e-2;
  static constexpr double a84 = 1.70383925712239993810214054705e-1;
  static constexpr double a85 = 1.07262030446373284651809199168e-1;
  static constexpr double a86 = -1.53194377486244017527936158236e-2;
  static constexpr double a87 = 8.27378916381402288758473766002e-3;
  static constexpr double a91 = 6.24110958716075717114429577812e-1;
  static constexpr double a94 = -3.36089262944694129406857109825e0;
  static constexpr double a95 = -8.68219346841726006818189891453e-1;
  static constexpr double a96 = 2.75920996994467083049415600797e1;
  static constexpr double a97 = 2.01540675504778934086186788979e1;
  static constexpr double a98 = -4.34898841810699588477366255144e1;
  static constexpr double a101 = 4.77662536438264365890433908527e-1;
  static constexpr double a104 = -2.48811461997166764192642586468e0;
  static constexpr double a105 = -5.90290826836842996371446475743e-1;
  static constexpr double a106 = 2.12300514481811942347288949897e1;
  static constexpr double a107 = 1.52792336328824235832596922938e1;
  static constexpr double a108 = -3.32882109689848629194453265587e1;
  static constexpr double a109 = -2.03312017085086261358222928593e-2;

  static constexpr double a111 = -9.3714243008598732571704021658e-1;
  static constexpr double a114 = 5.18637242884406370830023853209e0;
  static constexpr double a115 = 1.09143734899672957818500254654e0;
  static constexpr double a116 = -8.14978701074692612513997267357e0;
  static constexpr double a117 = -1.85200656599969598641566180701e1;
  static constexpr double a118 = 2.27394870993505042818970056734e1;
  static constexpr double a119 = 2.49360555267965238987089396762e0;
  static constexpr double a1110 = -3.0467644718982195003823669022e0;
  static constexpr double a121 = 2.27331014751653820792359768449e0;
  static constexpr double a124 = -1.05344954667372501984066689879e1;
  static constexpr double a125 = -2.00087205822486249909675718444e0;
  static constexpr double a126 = -1.79589318631187989172765950534e1;
  static constexpr double a127 = 2.79488845294199600508499808837e1;
  static constexpr double a128 = -2.85899827713502369474065508674e0;
  static constexpr double a129 = -8.87285693353062954433549289258e0;
  static constexpr double a1210 = 1.23605671757943030647266201528e1;
  static constexpr double a1211 = 6.43392746015763530355970484046e-1;

  static constexpr double a141 = 5.61675022830479523392909219681e-2;
  static constexpr double a147 = 2.53500210216624811088794765333e-1;
  static constexpr double a148 = -2.46239037470802489917441475441e-1;
  static constexpr double a149 = -1.24191423263816360469010140626e-1;
  static constexpr double a1410 = 1.5329179827876569731206322685e-1;
  static constexpr double a1411 = 8.20105229563468988491666602057e-3;
  static constexpr double a1412 = 7.56789766054569976138603589584e-3;
  static constexpr double a1413 = -8.298e-3;

  static constexpr double a151 = 3.18346481635021405060768473261e-2;
  static constexpr double a156 = 2.83009096723667755288322961402e-2;
  static constexpr double a157 = 5.35419883074385676223797384372e-2;
  static constexpr double a158 = -5.49237485713909884646569340306e-2;
  static constexpr double a1511 = -1.08347328697249322858509316994e-4;
  static constexpr double a1512 = 3.82571090835658412954920192323e-4;
  static constexpr double a1513 = -3.40465008687404560802977114492e-4;
  static constexpr double a1514 = 1.41312443674632500278074618366e-1;
  static constexpr double a161 = -4.28896301583791923408573538692e-1;
  static constexpr double a166 = -4.69762141536116384314449447206e0;
  static constexpr double a167 = 7.68342119606259904184240953878e0;
  static constexpr double a168 = 4.06898981839711007970213554331e0;
  static constexpr double a169 = 3.56727187455281109270669543021e-1;
  static constexpr double a1613 = -1.39902416515901462129418009734e-3;
  static constexpr double a1614 = 2.9475147891527723389556272149e0;
  static constexpr double a1615 = -9.15095847217987001081870187138e0;

  static constexpr double d41 = -0.84289382761090128651353491142e+01;
  static constexpr double d46 = 0.56671495351937776962531783590e+00;
  static constexpr double d47 = -0.30689499459498916912797304727e+01;
  static constexpr double d48 = 0.23846676565120698287728149680e+01;
  static constexpr double d49 = 0.21170345824450282767155149946e+01;
  static constexpr double d410 = -0.87139158377797299206789907490e+00;
  static constexpr double d411 = 0.22404374302607882758541771650e+01;
  static constexpr double d412 = 0.63157877876946881815570249290e+00;
  static constexpr double d413 = -0.88990336451333310820698117400e-01;
  static constexpr double d414 = 0.18148505520854727256656404962e+02;
  static constexpr double d415 = -0.91946323924783554000451984436e+01;
  static constexpr double d416 = -0.44360363875948939664310572000e+01;

  static constexpr double d51 = 0.10427508642579134603413151009e+02;
  static constexpr double d56 = 0.24228349177525818288430175319e+03;
  static constexpr double d57 = 0.16520045171727028198505394887e+03;
  static constexpr double d58 = -0.37454675472269020279518312152e+03;
  static constexpr double d59 = -0.22113666853125306036270938578e+02;
  static constexpr double d510 = 0.77334326684722638389603898808e+01;
  static constexpr double d511 = -0.30674084731089398182061213626e+02;
  static constexpr double d512 = -0.93321305264302278729567221706e+01;
  static constexpr double d513 = 0.15697238121770843886131091075e+02;
  static constexpr double d514 = -0.31139403219565177677282850411e+02;
  static constexpr double d515 = -0.93529243588444783865713862664e+01;
  static constexpr double d516 = 0.35816841486394083752465898540e+02;

  static constexpr double d61 = 0.19985053242002433820987653617e+02;
  static constexpr double d66 = -0.38703730874935176555105901742e+03;
  static constexpr double d67 = -0.18917813819516756882830838328e+03;
  static constexpr double d68 = 0.52780815920542364900561016686e+03;
  static constexpr double d69 = -0.11573902539959630126141871134e+02;
  static constexpr double d610 = 0.68812326946963000169666922661e+01;
  static constexpr double d611 = -0.10006050966910838403183860980e+01;
  static constexpr double d612 = 0.77771377980534432092869265740e+00;
  static constexpr double d613 = -0.27782057523535084065932004339e+01;
  static constexpr double d614 = -0.60196695231264120758267380846e+02;
  static constexpr double d615 = 0.84320405506677161018159903784e+02;
  static constexpr double d616 = 0.11992291136182789328035130030e+02;

  static constexpr double d71 = -0.25693933462703749003312586129e+02;
  static constexpr double d76 = -0.15418974869023643374053993627e+03;
  static constexpr double d77 = -0.23152937917604549567536039109e+03;
  static constexpr double d78 = 0.35763911791061412378285349910e+03;
  static constexpr double d79 = 0.93405324183624310003907691704e+02;
  static constexpr double d710 = -0.37458323136451633156875139351e+02;
  static constexpr double d711 = 0.10409964950896230045147246184e+03;
  static constexpr double d712 = 0.29840293426660503123344363579e+02;
  static constexpr double d713 = -0.43533456590011143754432175058e+02;
  static constexpr double d714 = 0.96324553959188282948394950600e+02;
  static constexpr double d715 = -0.39177261675615439165231486172e+02;
  static constexpr double d716 = -0.14972683625798562581422125276e+03;
};

template <int n, class F> struct StepperDopr853 : Dopr853_constants {
  using VecDoub = std::array<double, n>;
  static constexpr double EPS = std::numeric_limits<double>::epsilon();

  VecDoub dydx;
  VecDoub y, yout, yerr, yerr2, ytemp;
  VecDoub k2, k3, k4, k5, k6, k7, k8, k9, k10;
  double atol, rtol;
  double x, hnext, errold;
  bool reject;

  StepperDopr853(double atol_, double rtol_) : atol(atol_), rtol(rtol_) {}

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
      double err = error(h);
      if (success(err, h))
        break;
      if (std::abs(h) <= std::abs(x) * EPS)
        THROW("stepsize underflow in StepperDopr853");
    }
    derivs(x + h, yout, dydx);
    y = yout;
    x += h;
  }

  void dy(double h, F &&derivs) {
    int i;
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * a21 * dydx[i];
    derivs(x + c2 * h, ytemp, k2);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
    derivs(x + c3 * h, ytemp, k3);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a41 * dydx[i] + a43 * k3[i]);
    derivs(x + c4 * h, ytemp, k4);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a51 * dydx[i] + a53 * k3[i] + a54 * k4[i]);
    derivs(x + c5 * h, ytemp, k5);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a61 * dydx[i] + a64 * k4[i] + a65 * k5[i]);
    derivs(x + c6 * h, ytemp, k6);
    for (i = 0; i < n; i++)
      ytemp[i] =
          y[i] + h * (a71 * dydx[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
    derivs(x + c7 * h, ytemp, k7);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a81 * dydx[i] + a84 * k4[i] + a85 * k5[i] +
                             a86 * k6[i] + a87 * k7[i]);
    derivs(x + c8 * h, ytemp, k8);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a91 * dydx[i] + a94 * k4[i] + a95 * k5[i] +
                             a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
    derivs(x + c9 * h, ytemp, k9);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a101 * dydx[i] + a104 * k4[i] + a105 * k5[i] +
                             a106 * k6[i] + a107 * k7[i] + a108 * k8[i] +
                             a109 * k9[i]);
    derivs(x + c10 * h, ytemp, k10);
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a111 * dydx[i] + a114 * k4[i] + a115 * k5[i] +
                             a116 * k6[i] + a117 * k7[i] + a118 * k8[i] +
                             a119 * k9[i] + a1110 * k10[i]);
    derivs(x + c11 * h, ytemp, k2);
    double xph = x + h;
    for (i = 0; i < n; i++)
      ytemp[i] = y[i] + h * (a121 * dydx[i] + a124 * k4[i] + a125 * k5[i] +
                             a126 * k6[i] + a127 * k7[i] + a128 * k8[i] +
                             a129 * k9[i] + a1210 * k10[i] + a1211 * k2[i]);
    derivs(xph, ytemp, k3);
    for (i = 0; i < n; i++) {
      k4[i] = b1 * dydx[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] +
              b10 * k10[i] + b11 * k2[i] + b12 * k3[i];
      yout[i] = y[i] + h * k4[i];
    }
    for (i = 0; i < n; i++) {
      yerr[i] = k4[i] - bhh1 * dydx[i] - bhh2 * k9[i] - bhh3 * k3[i];
      yerr2[i] = er1 * dydx[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] +
                 er9 * k9[i] + er10 * k10[i] + er11 * k2[i] + er12 * k3[i];
    }
  }

  double error(double h) {
    double err = 0.0, err2 = 0.0, sk, deno;
    for (int i = 0; i < n; i++) {
      sk = atol + rtol * std::max(std::abs(y[i]), std::abs(yout[i]));
      double t1 = yerr[i] / sk;
      err2 += t1 * t1;
      double t2 = yerr2[i] / sk;
      err += t2 * t2;
    }
    deno = err + 0.01 * err2;
    if (deno <= 0.0)
      deno = 1.0;
    return std::abs(h) * err * std::sqrt(1.0 / (n * deno));
  }

  bool success(const double err, double &h) {
    static const double beta = 0.4 / 8.0, alpha = 1.0 / 8.0 - beta * 0.2,
                        safe = 0.9, minscale = 0.333, maxscale = 6.0;
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
      return false;
    }
  }
};

#undef THROW
