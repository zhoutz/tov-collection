#pragma once

#include <vector>

inline std::vector<double> linspace(double start, double end, int num_points) {
  std::vector<double> vec(num_points);
  double step = (end - start) / (num_points - 1);
  for (int i = 0; i < num_points; ++i) {
    vec[i] = start + i * step;
  }
  return vec;
}

inline std::vector<double> geomspace(double start, double end, int num_points) {
  std::vector<double> vec(num_points);
  double log_start = std::log(start);
  double log_end = std::log(end);
  double step = (log_end - log_start) / (num_points - 1);
  for (int i = 0; i < num_points; ++i) {
    vec[i] = std::exp(log_start + i * step);
  }
  return vec;
}
