#include<iostream>
#include<fstream>

#include <vector>

void Solve(const std::vector<double> & Y, std::vector<double> &U,\
  const double & dy, const double & dt, const double & pho, const double & nu,\
    const size_t & nIter, const double & pDrop);
