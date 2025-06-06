#include "bench.hpp"
#include "constants.hpp"
#include "eos.hpp"
#include "stepperdopr5.hpp"
#include "tov_uvy_p.hpp"
#include <array>

double rtol;

auto c_MR_point(TOV_uvy_p &tov, double e_c) {
  double r0 = 1e-6;
  double p0 = tov.eos.p_of_e(e_c);
  double m0 = C2 / 3. * r0 * r0 * r0 * e_c;
  double u0 = r0 * r0;
  double v0 = m0 / r0;
  double y0 = 2.;
  double dp0 = -p0 * 1e-8;

  d3 uvy{u0, v0, y0};
  d3 duvy_dp;
  tov(p0, uvy, duvy_dp);

  StepperDopr5<3, decltype(tov)> stepper(0, rtol);
  stepper.set_init(uvy, duvy_dp, p0, dp0);

  stepper.integrate_to(p_boundary, tov);

  double ret_r = std::sqrt(stepper.y[0]);
  double ret_m = ret_r * stepper.y[1];
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
  TOV_uvy_p tov{eos};
  bench_c_MR_point(tov, c_MR_point, out_fname);
}
