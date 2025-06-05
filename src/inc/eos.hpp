#pragma once

#include "linear_interp.hpp"
#include <fstream>
#include <iostream>
#include <vector>

struct EOS {
  std::vector<double> e, p, cs2, h;
  LinearInterp ee, pp, hh;

  void read_file_natual(const std::string &filename) {
    std::ifstream fin(filename);
    if (!fin) {
      std::cerr << "Error: cannot open file " << filename << std::endl;
      std::exit(1);
    }

    double e_, p_, cs2_, h_;
    while (fin >> e_ >> p_ >> cs2_ >> h_) {
      e.push_back(e_);
      p.push_back(p_);
      cs2.push_back(cs2_);
      h.push_back(h_);
    }
    fin.close();
    e.shrink_to_fit();
    p.shrink_to_fit();
    cs2.shrink_to_fit();
    h.shrink_to_fit();

    std::cout << "<< " << filename << " (" << e.size() << ")\n";

    ee.reset(e.data(), e.size());
    pp.reset(p.data(), p.size());
    hh.reset(h.data(), h.size());
  }

  double p_of_e(double e_) {
    int i = ee.hunt(e_);
    return ee.lin_interp(e_, e[i], e[i + 1], p[i], p[i + 1]);
  }

  auto e_cs2_of_p(double p_) {
    int i = pp.hunt(p_);
    return std::make_pair(
        pp.lin_interp(p_, p[i], p[i + 1], e[i], e[i + 1]),
        pp.lin_interp(p_, p[i], p[i + 1], cs2[i], cs2[i + 1]));
  }
  double h_of_p(double p_) {
    int i = pp.hunt(p_);
    return pp.lin_interp(p_, p[i], p[i + 1], h[i], h[i + 1]);
  }

  auto e_p_cs2_of_h(double h_) {
    int i = hh.hunt(h_);
    return std::make_tuple(
        hh.lin_interp(h_, h[i], h[i + 1], e[i], e[i + 1]),
        hh.lin_interp(h_, h[i], h[i + 1], p[i], p[i + 1]),
        hh.lin_interp(h_, h[i], h[i + 1], cs2[i], cs2[i + 1]));
  }
};
