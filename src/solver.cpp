#include "solver.h"

#include "functions.h"

void SolveLam(const std::vector<double> & Y, std::vector<double> &U, std::vector<double> &nuT,\
  const double & H, const double & dy, const double & dt, const double & rho, const double & nu,\
    const size_t & nIter, const double & pDrop, const double & Uw){

    std::size_t N = Y.size();

    std::vector<double> Un(N);

    std::vector<double> A(N);
    std::vector<double> B(N);
    std::vector<double> C(N);
    std::vector<double> F(N);


    double kapa = 0.41;
    double y_plus_w, y_plus, y_star, U_star;

    std::vector<double> resU(N, 0.0); //Velocity residuals
    std::ofstream  residuals("res.plt");  //Residuals output file

    if (residuals.is_open()) {
        for(size_t i = 0; i < nIter; i++){
          U_star = std::pow(nu * (U[1] - U[0]) / rho / dy, 0.5);
          y_plus_w = U_star * dy / nu;
          std::cout << y_plus_w << " " << U_star << std::endl;

            nuT[0] = 0;
            nuT[N - 1] = 0;

            for(size_t j = 1; j < N - 1; ++j){
                y_star = min(Y[j], H - Y[j]);
                y_plus = y_star * U_star / nu;
                nuT[j] = std::pow(kapa * y_star, 2) * std::abs((U[j + 1] - U[j - 1])) /2.0 / dy * std::pow((1 - std::exp(-y_plus / 26)), 2);
//                std::cout << j << std::pow(kapa * y_star, 2) * std::abs((U[j + 1] - U[j - 1])) /2.0 / dy * std::pow((1 - std::exp(-y_plus / 26)), 2)  << " " << std::abs((U[j] - U[j - 1])) << std::endl;
                  //nuT[j] = 0;
            }
//Implicit

            A[0] = 0;
            C[0] = 1;
            B[0] = 0;
            F[0] = 0;

            A[N - 1] = 0;
            C[N - 1] = 1;
            B[N - 1] = 0;
            F[N - 1] = Uw;

            for(size_t j = 1; j < N - 1; j++){
              A[j] = - (2 * nu + nuT[j] + nuT[j - 1]) / 2 / dy / dy;
              //C[j] = 1 / dt + (2 * nu + nuT[j + 1] + nuT[j - 1])/ dy / dy;
              C[j] = 1 / dt + (4 * nu + nuT[j + 1] + 2 * nuT[j] + nuT[j - 1]) / 2 / dy / dy;

              B[j] =  - (2 * nu + nuT[j + 1] + nuT[j]) /2 / dy / dy;;
              F[j] = 1 / dt * U[j] - pDrop / rho;
              //F[j] = - pDrop / rho;

            }

            solveMatrix(A, B, C, F, Un);

            for(size_t j = 0; j < N; j++){
              resU[j] = (U[j] - Un[j]) / max(Un[j], 1e-10);
              U[j] = Un[j];
            }

//Explicit
      /*      U[0] = 0;
            U[N - 1] = Uw;
            for(size_t j = 1; j < N - 1; j++){
              resU[j] = dt * (-std::pow(U_star, 2) * 2 / H + std::pow(dy, -2) * ((nu + nuT[j + 1]) * (U[j + 1] - U[j]) - (nu + nuT[j]) * (U[j] - U[j - 1])));
              U[j] += resU[j];
            }
*/            double maxResU = vecMaxAbs(resU);

            std::cout << i << " " << maxResU <<  " " << vecMaxAbs(nuT) << std::endl;
            residuals << i << " " << maxResU << "\n";



            if (maxResU < 1e-15) break;
        }

        residuals.close();
      }

      else {
          std::cout << "Invalid output file" << std::endl;
          return;
      }

      return;

    }
