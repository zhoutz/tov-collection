#pragma once

#include "constants.hpp"
#include "vec.hpp"
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

inline unsigned g_cnt;

template <typename TOV_type, typename Func>
void bench_c_MR_point(TOV_type &&tov, Func &&c_MR_point,
                      std::string const &out_fname) {
  std::vector<double> e_c_vec = geomspace(e_c_lb, e_c_ub, n_e_c);

  using mr_type = decltype(std::function{c_MR_point})::result_type;
  std::vector<mr_type> vec;

  g_cnt = 0;
  auto time_start = std::chrono::high_resolution_clock::now();
  for (double e_c : e_c_vec) {
    vec.emplace_back(c_MR_point(tov, e_c));
    // break; // Remove this line to run the full benchmark
  }
  auto time_end = std::chrono::high_resolution_clock::now();
  double t_us =
      std::chrono::duration<double, std::micro>(time_end - time_start).count();

  {
    std::ofstream of(out_fname);
    if (!of) {
      std::cerr << "Error opening file" << std::endl;
      exit(1);
    }
    of << std::setprecision(15);
    of << t_us << " " << g_cnt << "\n";
    for (auto const &mrt : vec) {
      for (auto const &m : mrt) {
        of << m << " ";
      }
      of << "\n";
    }

    of.close();
    std::cout << ">> " << out_fname << std::endl;
  }
}
