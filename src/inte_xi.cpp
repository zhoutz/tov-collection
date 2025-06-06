#include "bench.hpp"
#include "constants.hpp"
#include "eos.hpp"
#include "stepperdopr5.hpp"
#include "tov_mrt_xi.hpp"
#include <array>

double rtol;

auto c_MR_point(TOV_mrt_xi &tov, double e_c) {
  double p_c = tov.eos.p_of_e(e_c);
  double y_c = 2.;

  double r0 = 1e-6;
  double w0 = r0 * r0 * r0;
  double m0 = C2 / 3. * w0 * e_c;
  double y0 = y_c;
  double xi0 = std::log(p_c);

  double dxi = -xi0 * 1e-14;
  double xi_boundary = std::log(p_boundary);

  d3 wmy{w0, m0, y0};
  d3 dwmy_dxi;
  tov(xi0, wmy, dwmy_dxi);

  StepperDopr5<3, decltype(tov)> stepper(0, rtol);
  stepper.set_init(wmy, dwmy_dxi, xi0, dxi);

  stepper.integrate_to(xi_boundary, tov);

  double ret_m = stepper.y[1];
  double ret_r = std::cbrt(stepper.y[0]);
  double ret_y = stepper.y[2];
  return std::array<double, 3>{ret_m, ret_r, ret_y};
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " in_fname out_fname rtol"
              << std::endl;
    exit(1);
  }
  std::string in_fname = argv[1];
  std::string out_fname = argv[2];
  rtol = std::stod(argv[3]);

  EOS eos;
  eos.read_file_natual(in_fname);
  TOV_mrt_xi tov{eos};
  bench_c_MR_point(tov, c_MR_point, out_fname);
}
