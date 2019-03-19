#include <vector>
#include <cmath>

#include "functions.h"

double max(const double & a, const double & b){
  if (a > b)
    return a;
  else
    return b;
}

double min(const double & a, const double & b){
  if (a < b)
    return a;
  else
    return b;
}


double vecMaxAbs(const std::vector<double> & Y){
    double M = Y[0];
    for(std::size_t k = 0; k < Y.size(); k++){
        if (std::abs(Y[k]) > M)
            M = std::abs(Y[k]);
    }
    return M;
}

double vecAverage(const std::vector<double> & Y){
    double sum = 0;
    for(std::size_t k = 0; k < Y.size(); k++)
        sum += Y[k];
    return (sum / Y.size());
}


void solveMatrix(const std::vector<double> & A, \
      const std::vector<double> & B, std::vector<double> & C,\
      std::vector<double> & F, std::vector<double> & X){
    std::size_t N = A.size();
    double m;

    for (std::size_t i = 1; i < N; ++i) {

        m = A[i] / C[i - 1];
        C[i] -= m * B[i - 1];
        F[i] -= m * F[i - 1];
    }


    X[N - 1] = F[N - 1] / C[N - 1];

    for (int i = (N - 2); i >= 0; --i) {
        X[i] = (F[i] - B[i] * X[i + 1]) / C[i];
    }

    return;
}
