#include "solver.h"

#include "functions.h"

void Solve(const std::vector<double> & Y, std::vector<double> &U,\
  const double & dy, const double & dt, const double & rho, const double & nu, \
    const size_t & nIter, const double & pDrop){

    std::vector<double> res(U.size());

    for (size_t j = 0; j != res.size(); j++)
        res[j] = 0;

    std::ofstream  residuals("res.csv");

    if (residuals.is_open()) {

        for(size_t i = 0; i < nIter; i++){
            for(size_t j = 1; j < U.size() - 1; j++){
                res[j] = dt * (-pDrop / rho + nu  / dy / dy * (U[j + 1] + U[j - 1] - 2 * U[j]));

                U[j] += res[j];
            }
            double maxRes = vecMaxAbs(res);
            std::cout << i << " " << maxRes << std::endl;
            residuals << i << " " << maxRes << "\n";

            if (maxRes < 1e-14) break;
        }
        residuals.close();
      }

      else {
          std::cout << "Invalid output file" << std::endl;
          return;
      }

      return;

    }
