#include "solver.h"

#include "functions.h"

void SolveLam(const std::vector<double> & Y, std::vector<double> &U, \
  const double & dy, const double & dt, const double & rho, const double & nu,\
    const size_t & nIter, const double & pDrop, const double & Uw){

    std::size_t N = Y.size();

    std::vector<double> Un(N, 0.0);

    std::vector<double> A(N, 0.0);
    std::vector<double> B(N, 0.0);
    std::vector<double> C(N, 0.0);
    std::vector<double> F(N, 0.0);


    std::vector<double> resU(N, 0.0); //Velocity residuals
    std::ofstream  residuals("res.dat");  //Residuals output file

    if (residuals.is_open()) {

        for(size_t i = 0; i < nIter; i++){

            A[0] = 0;
            C[0] = 1;
            B[0] = 0;
            F[0] = 0;

            A[N - 1] = 0;
            C[N - 1] = 1;
            B[N - 1] = 0;
            F[N - 1] = Uw;

            for(size_t j = 1; j < N - 1; j++){
              A[j] = - nu / dy / dy;
              C[j] = 1 / dt + 2 * nu / dy / dy;
              B[j] =  - nu / dy / dy;
              F[j] = 1 / dt * U[j] - pDrop / rho;
            }


            solveMatrix(A, B, C, F, Un);

            for(size_t j = 0; j < N; j++){
              resU[j] = (U[j] - Un[j]) / max(Un[j], 1e-14);
              U[j] = Un[j];
            }

            double maxResU = vecMaxAbs(resU);

            std::cout << i << " " << maxResU << std::endl;
            residuals << i << " " << maxResU << "\n";


            if (maxResU < 1e-14) break;
        }
        //std::cout << i << " " << maxResU << std::endl;

        residuals.close();
      }

      else {
          std::cout << "Invalid output file" << std::endl;
          return;
      }

      return;

    }
