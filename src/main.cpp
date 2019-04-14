
#include<iostream>
#include<fstream>

#include <vector>
#include <string>

#include "solver.h"
#include "functions.h"

int main(void) {

    double H; //High of channel
    double rho; //Fluid Density
    double nu;  //Fluid viscosity

    double Uw;  //Velocity of upper wall
    double Re_tau; //Velocity of upper wall

    size_t NJ;  //Number of nodes along j-direction

    double dt;  //Time step size
    size_t nIter; //Max number of iteration
    std::string outFile;
//********************** Read Input Data **************************************

    std::fstream  inData("inData.dat");
    std::cout << "Reading input data from inData.dat " << std::endl;
    if (inData.is_open()) {

        inData >> H;
        inData >> rho;
        inData >> nu;
        inData >> Uw;
        inData >> Re_tau;
        inData >> NJ;
        inData >> dt;
        inData >> nIter;
        inData >> outFile;
        inData.close();
    }
    else {
        std::cout << "Invalid input file" << std::endl;
        return 0;
    }

//********************** Calculate patameters of flow and solution **************************************

    double U_tau = 2 * nu * Re_tau / H;

    std::cout << "U_tau " << U_tau << std::endl;

    double Re; //Reynolds number
    Re = 14.64 * std::pow(Re_tau, 8.0 / 7.0);

    std::cout << "Re " << Re << std::endl;

    double pDrop; //Pressure drop per length unit

    if (Uw > 0)
        pDrop = 0;
    else
        pDrop = -2 * U_tau * U_tau / H;

    std::cout << "Pressure drop per length unis, Pa/m " << pDrop << std::endl;

    std::cout << "Velocity of upper wall, m/s " << Uw << std::endl;

    double dy;  //Length of cell
    dy = H / (NJ - 1);

//********************** Create mesh **************************************
    std::cout << "Creating mesh..." << std::endl;

    std::vector<double> Y(NJ);  //Y-coordinates of mesh nodes

    for (size_t j = 0; j != NJ; j++)
        Y[j] = j * dy;

//********************** Allocate & Initiate All Arrays  **********************************
    //std::cout << "Initiating fields..." << std::endl;

    std::vector<double> U(NJ, 0.0); //X component of Velocity

    U[0] = 0;
    U[NJ - 1] = Uw;

    std::vector<double> nuT(NJ, 0.0); //Turbulence viscosity

    std::vector<double> K(NJ, 0.0); //Turbulence kinetic energy

    std::vector<double> W(NJ, 0.0); //Dissipation velocity

//********************** Calclulate Initial Approach ***************************
    std::cout << "Calculating initial approach..." << std::endl;

    SolvePrandtl(Y, U, nuT, H, dy, dt, rho, nu, nIter, pDrop, Uw);

    K[0] = 0;
    K[NJ - 1] = 0;

    W[0] = 60 * nu / 0.075 / std::pow(dy, 2);
    W[NJ - 1] = 60 * nu / 0.075 / std::pow(dy, 2);

    for (size_t j = 1; j != NJ - 1; j++){
        K[j] = nuT[j] / 0.3 * std::abs((U[j + 1] - U[j - 1]) / 2 / dy);
        W[j] = K[j] / nuT[j];
    }

    std::ofstream  out("outPr.plt");

    std::cout << "Writing approach fields in file..." << std::endl;

    if (out.is_open()) {
        for(size_t j = 0; j != NJ; j++){
            out << Y[j] << " " << U[j]  << " " << nuT[j] << " " << K[j]  << " " << W[j] << "\n";
        }
        out.close();
    }
    else {
        std::cout << "Invalid output file" << std::endl;
        return 0;
    }


//********************** Solve equations ***************************

    std::cout << "Solving Reynold's equations..." << std::endl;

    SolveWilcox(Y, U, nuT, K, W, H, dy, dt, rho, nu, nIter, pDrop, Uw);

    double U_star = std::pow(nu * (- U[2] + 4 * U[1] - 3 * U[0]) / 2 / dy, 0.5);
    double y_plus_w = U_star * dy / nu;
    std::cout << "U_tau " << U_star << std::endl;

    std::cout << "Y+ " << y_plus_w << std::endl;

//********************** Output Fields *****************************************
    std::ofstream  outData(outFile);

    std::cout << "Writing fields in file..." << std::endl;

    if (outData.is_open()) {
        for(size_t j = 0; j != NJ; j++){
            outData << Y[j] << " " << U[j]  << " " << nuT[j] << " " << K[j]  << " " << W[j] << "\n";
        }
        outData.close();
    }
    else {
        std::cout << "Invalid output file" << std::endl;
        return 0;
    }

    return 0;
}
