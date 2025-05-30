#include "bench.hpp"
#include "constants.hpp"
#include "eos.hpp"
#include "stepperdopr5.hpp"
#include "tov_mrt_r.hpp"
#include <array>

double rtol;

auto c_MR_point(TOV_mrt_r &tov, double e_c) {
  double r0 = 1e-6;
  double dr = 1e-4;
  double p_c = tov.eos.p_of_e(e_c);
  double m_c = C2 / 3. * r0 * r0 * r0 * e_c;
  double y_c = 2.;

  d3 pmy{p_c, m_c, y_c};
  d3 dpmy_dr;
  tov(r0, pmy, dpmy_dr);

  StepperDopr5<3, decltype(tov)> stepper(0, rtol);
  stepper.set_init(pmy, dpmy_dr, r0, dr);

  while (stepper.y[0] > p_boundary) {
    stepper.step(tov);
  }

  int n_newton_raphson = 1;
  while (n_newton_raphson--) {
    dr = -(stepper.y[0] - p_boundary) / stepper.dydx[0];
    stepper.hnext = dr;
    stepper.integrate_to(stepper.x + dr, tov);
  }

  double ret_m = stepper.y[1];
  double ret_r = stepper.x;
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

  EOS_epcs eos;
  eos.read_file_natual(in_fname);
  TOV_mrt_r tov{eos};
  bench_c_MR_point(tov, c_MR_point, out_fname);
}
