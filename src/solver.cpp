#include "solver.h"

#include "functions.h"

void SolvePrandtl(const std::vector<double> & Y, std::vector<double> &U, std::vector<double> &nuT,\
    const double & H, const double & dy, const double & dt, const double & rho, const double & nu,\
        const size_t & nIter, const double & pDrop, const double & Uw){

    std::size_t N = Y.size();

    std::vector<double> Un(N);

    std::vector<double> A(N);
    std::vector<double> B(N);
    std::vector<double> C(N);
    std::vector<double> F(N);


    double kapa = 0.41;
    double a = 26;

    double y_plus_w, y_plus, y_tilde, U_star;

    std::vector<double> resU(N, 0.0); //Velocity residuals
    std::ofstream  residuals("res_Pr.plt");  //Residuals output file

    if (residuals.is_open()) {
        for(size_t i = 0; i != nIter; i++){
            U_star = std::pow(nu * (- U[2] + 4 * U[1] - 3 * U[0]) / 2 / dy, 0.5);
            y_plus_w = U_star * (Y[1] - Y[0]) / nu;

            nuT[0] = 0;
            nuT[N - 1] = 0;

            for(size_t j = 1; j != N - 1; ++j){
                y_tilde = min(Y[j], H - Y[j]);
                y_plus = y_tilde * U_star / nu;
                nuT[j] = std::pow(kapa * y_tilde, 2) * std::pow((1 - std::exp(-y_plus / a)), 2) *\
                            std::abs((U[j + 1] - U[j - 1]) / 2 / dy);
            }

            A[0] = 0;
            C[0] = 1;
            B[0] = 0;
            F[0] = 0;

            A[N - 1] = 0;
            C[N - 1] = 1;
            B[N - 1] = 0;
            F[N - 1] = Uw;

            for(size_t j = 1; j != N - 1; j++){
                A[j] = - (2 * nu + nuT[j] + nuT[j - 1]) / 2 / dy / dy;
                B[j] = - (2 * nu + nuT[j + 1] + nuT[j]) / 2 / dy / dy;

                C[j] = 1 / dt - A[j] - B[j];

                F[j] = 1 / dt * U[j] - pDrop;
            }

            solveMatrix(A, B, C, F, Un);

            for(size_t j = 0; j != N; j++){
                resU[j] = (U[j] - Un[j]) / max(Un[j], 1e-10);
                U[j] = Un[j];
            }

            double maxResU = vecMaxAbs(resU);

            residuals << i << " " << maxResU << "\n";
            //std::cout << i << " " << maxResU << "\n";

            if (maxResU < 1e-12) break;
          }

          residuals.close();
    }

    else {
          std::cout << "Invalid output file" << std::endl;
          return;
    }

    return;

}




void SolveWilcox(const std::vector<double> &Y, std::vector<double> &U, std::vector<double> &nuT,\
      std::vector<double> &K, std::vector<double> &W,\
          const double & H, const double & dy, const double & dt, const double & rho, const double & nu,\
            const size_t & nIter, const double & pDrop, const double & Uw){

    std::size_t N = U.size();

    std::vector<double> Un(N);
    std::vector<double> Kn(N);
    std::vector<double> Wn(N);

    std::vector<double> A(N);
    std::vector<double> B(N);
    std::vector<double> C(N);
    std::vector<double> F(N);

    std::vector<double> resU(N, 0.0); //Velocity residuals
    std::vector<double> resK(N, 0.0); //K residuals
    std::vector<double> resW(N, 0.0); //W residuals

    std::vector<double> alpha_star(N, 0.0);
    std::vector<double> beta_star(N, 0.0);
    std::vector<double> gamma(N, 0.0);

    //Model Constants
    double sigma_k = 0.6;
    double sigma_w = 0.5;
    double sigma_d;
    double dkdw;

    double c_lim = 7.0 / 8.0;

    double beta_0 = 0.0708;
    double beta_star_0 = 0.09;
    double beta = beta_0;

    double r_beta = 8.0;
    double r_k = 6.0;
    double r_w = 2.61;
    double Re_t;

    double alpha_0 = 1.0 / 9.0;
    double alpha_star_0 = beta_0 / 3.0;

    double omega_circum;

    std::ofstream  residuals("res.plt");  //Residuals output file

    if (residuals.is_open()) {
        for(size_t i = 0; i != nIter; i++){
            //Solve Velocity Eq
            A[0] = 0;
            C[0] = 1;
            B[0] = 0;
            F[0] = 0;

            A[N - 1] = 0;
            C[N - 1] = 1;
            B[N - 1] = 0;
            F[N - 1] = Uw;

            for(size_t j = 1; j != N - 1; j++){
                A[j] = - (2 * nu + nuT[j] + nuT[j - 1]) / 2 / dy / dy;
                B[j] = - (2 * nu + nuT[j + 1] + nuT[j]) / 2 / dy / dy;

                C[j] = 1 / dt - A[j] - B[j];

                F[j] = 1 / dt * U[j] - pDrop;
            }

            solveMatrix(A, B, C, F, Un);

            for(size_t j = 0; j != N; j++){
                resU[j] = (U[j] - Un[j]) / max(Un[j], 1e-10);
                U[j] = Un[j];
            }

            for(size_t j = 0; j < N; j++){
                Re_t = K[j] / W[j] / nu;
                alpha_star[j] = (alpha_star_0 + Re_t / r_k) / (1 + Re_t / r_k);
                //alpha_star[j] = 1;
                beta_star[j] = 0.09 * (100.0 * beta_0 / 27.0 + pow(Re_t / r_beta, 4)) / (1 + pow(Re_t / r_beta, 4));
                //beta_star[j] = 0.09;
                gamma[j] = 13.0 / 25.0 * (alpha_0 + Re_t / r_w) / (1 + Re_t / r_w) / alpha_star[j];
                //gamma[j] = 13.0 / 25.0;
            }


                //Solve K Eq
            A[0] = 0;
            C[0] = 1;
            B[0] = 0;
            F[0] = 0;

            A[N - 1] = 0;
            C[N - 1] = 1;
            B[N - 1] = 0;
            F[N - 1] = 0;

            for(size_t j = 1; j != N - 1; j++){
                B[j] = - (2 * nu + sigma_k * alpha_star[j + 1] * K[j + 1] / W[j + 1] \
                                 + sigma_k * alpha_star[j] * K[j] / W[j]) / 2 / dy / dy;

                A[j] = - (2 * nu + sigma_k * alpha_star[j - 1] * K[j - 1] / W[j - 1] \
                                 + sigma_k * alpha_star[j] * K[j] / W[j]) / 2 / dy / dy;

                C[j] = 1 / dt - A[j] - B[j] + beta_star[j] * W[j];

                F[j] = 1 / dt * K[j] + 2 * nuT[j] * pow((U[j + 1] - U[j - 1]) / 2 / dy, 2);
            }

            solveMatrix(A, B, C, F, Kn);

            for(size_t j = 0; j != N; j++){
                resK[j] = (K[j] - Kn[j]) / max(Kn[j], 1e-10);
                K[j] = Kn[j];
            }

                  //Solve W Eq
            A[0] = 0;
            C[0] = 1;
            B[0] = 0;
            F[0] =  60 * nu / 0.075 / dy / dy;

            A[N - 1] = 0;
            C[N - 1] = 1;
            B[N - 1] = 0;
            F[N - 1] = 60 * nu / 0.075 / dy / dy;

            for(size_t j = 1; j != N - 1; j++){
                dkdw = (K[j + 1] - K[j - 1]) * (W[j + 1] - W[j - 1]) / 4 / dy / dy;
                if(dkdw > 0 )
                    sigma_d = 1.0 / 8.0;
                else
                    sigma_d = 0.0;

                B[j] = - (2 * nu + sigma_w * alpha_star[j + 1] * K[j + 1] / W[j + 1] \
                                 + sigma_w * alpha_star[j] * K[j] / W[j]) / 2 / dy / dy;

                A[j] = - (2 * nu + sigma_w * alpha_star[j - 1] * K[j - 1] / W[j - 1]\
                                 + sigma_w * alpha_star[j] * K[j] / W[j]) / 2 / dy / dy;

                C[j] = 1 / dt - A[j] - B[j] + beta * W[j];

                F[j] = 1 / dt * W[j] + 2 * nuT[j] * gamma[j] * W[j] / K[j] * pow((U[j + 1] - U[j - 1]) / 2 / dy, 2)\
                     + sigma_d / W[j] * dkdw;

              }

              solveMatrix(A, B, C, F, Wn);


              for(size_t j = 0; j != N; j++){
                  resW[j] = (W[j] - Wn[j]) / max(Wn[j], 1e-10);
                  W[j] = Wn[j];
              }

              nuT[0] = 0;
              nuT[N - 1] = 0;

              for(size_t j = 1; j != N - 1; j++){
                  omega_circum = max(W[j], c_lim * std::abs(U[j + 1] - U[j - 1]) \
                                  / 2 / dy * std::pow(alpha_star[j] / beta_star_0, 0.5));
                  nuT[j] = alpha_star[j] * K[j] / omega_circum;
              }

              double maxResU = vecMaxAbs(resU);
              double maxResK = vecMaxAbs(resK);
              double maxResW = vecMaxAbs(resW);

             //std::cout << i << " " << maxResU <<  " " << maxResK << " " << maxResW << std::endl;
             residuals << i << " " << maxResU <<  " " << maxResK << " " << maxResW << "\n";

              if (max(maxResU, max(maxResK, maxResW)) < 1e-10) break;
        }

        residuals.close();
    }

    else {
              std::cout << "Invalid output file" << std::endl;
              return;
          }

    return;

}
