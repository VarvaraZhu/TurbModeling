#include<iostream>
#include<fstream>
#include<cmath>
#include <vector>

void SolvePrandtl(const std::vector<double> & Y, std::vector<double> &U, std::vector<double> &nuT,\
  const double & H, const double & dy, const double & dt, const double & rho, const double & nu,\
    const size_t & nIter, const double & pDrop, const double & Uw);


void SolveWilcox(const std::vector<double> & Y, std::vector<double> &U, std::vector<double> &nuT,\
      std::vector<double> &K, std::vector<double> &W,\
      const double & H, const double & dy, const double & dt, const double & rho, const double & nu,\
        const size_t & nIter, const double & pDrop, const double & Uw);
