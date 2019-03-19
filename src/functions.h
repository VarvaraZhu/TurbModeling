#include<iostream>
#include <vector>

double max(const double & a, const double & b);

double min(const double & a, const double & b);

double vecMaxAbs(const std::vector<double> & Y);

double vecAverage(const std::vector<double> & Y);

void solveMatrix(const std::vector<double> & A, \
      const std::vector<double> & B, std::vector<double> & C,\
       std::vector<double> & F, std::vector<double> & X);
