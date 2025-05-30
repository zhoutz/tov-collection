#pragma once

#include "linear_interp.hpp"
#include <fstream>
#include <iostream>
#include <vector>

struct EOS_epcs {
  std::vector<double> e, p, cs2;
  LinearInterp ee, pp;

  void read_file_natual(const std::string &filename) {
    std::ifstream fin(filename);
    if (!fin) {
      std::cerr << "Error: cannot open file " << filename << std::endl;
      std::exit(1);
    }

    double e_, p_, cs2_;
    while (fin >> e_ >> p_ >> cs2_) {
      e.push_back(e_);
      p.push_back(p_);
      cs2.push_back(cs2_);
    }
    fin.close();
    e.shrink_to_fit();
    p.shrink_to_fit();
    cs2.shrink_to_fit();

    std::cout << "<< " << filename << " (" << e.size() << ")\n";

    ee.reset(e.data(), e.size());
    pp.reset(p.data(), p.size());
  }

  double p_of_e(double e_) {
    int i = ee.hunt(e_);
    return ee.linterp(e_, e[i], e[i + 1], p[i], p[i + 1]);
  }
  auto e_cs2_of_p(double p_) {
    int i = pp.hunt(p_);
    return std::make_pair(pp.linterp(p_, p[i], p[i + 1], e[i], e[i + 1]),
                          pp.linterp(p_, p[i], p[i + 1], cs2[i], cs2[i + 1]));
  }
};
