#include<iostream>
#include<fstream>

#include <vector>

void SolveLam(const std::vector<double> & Y, std::vector<double> &U, \
  const double & dy, const double & dt, const double & rho, const double & nu,\
    const size_t & nIter, const double & pDrop, const double & Uw);
