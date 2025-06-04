#include "bench.hpp"
#include "constants.hpp"
#include "eos.hpp"
#include "stepperdopr5.hpp"
#include "tov_mrt_h.hpp"
#include <array>

double rtol;

auto c_MR_point(TOV_mrt_h &tov, double e_c) {
  double r0 = 1e-6;
  double p0 = tov.eos.p_of_e(e_c);
  double h0 = tov.eos.h_of_p(p0);
  double m0 = C2 / 3. * r0 * r0 * r0 * e_c;
  double u0 = r0 * r0;
  double v0 = m0 / r0;
  double y0 = 2.;
  double dh0 = -h0 * 1e-8;
  double h_boundary = tov.eos.h_of_p(p_boundary);

  d3 uvy{u0, v0, y0};
  d3 duvy_dh;
  tov(h0, uvy, duvy_dh);

  StepperDopr5<3, decltype(tov)> stepper(0, rtol);
  stepper.set_init(uvy, duvy_dh, h0, dh0);

  stepper.integrate_to(h_boundary, tov);

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
  TOV_mrt_h tov{eos};
  bench_c_MR_point(tov, c_MR_point, out_fname);
}
