
#include<iostream>
#include<fstream>

#include <vector>

#include "solver.h"

int main(void) {

    double H; //High of channel
    double rho; //Fluid Density
    double nu;  //Fluid viscosity

    double Uw;  //Velocity of upper wall
    double Uref; //Velocity of upper wall

    size_t NJ;  //Number of nodes along j-direction

    double dt;  //Time step size
    size_t nIter; //Max number of iteration

//********************** Read Input Data **************************************
    std::fstream  inData("inData.dat");
    std::cout << "Reading input data from inData.dat " << std::endl;
    if (inData.is_open()) {

        inData >> H;
        inData >> rho;
        inData >> nu;
        inData >> Uw;
        inData >> Uref;
        inData >> NJ;
        inData >> dt;
        inData >> nIter;

        inData.close();
    }
    else {
        std::cout << "Invalid input file" << std::endl;
        return 0;
    }

//********************** Calculate patameters of flow and solution **************************************

    double Re; //Reynolds number
    Re = Uref * H * rho / nu;

    std::cout << "Re " << Re << std::endl;

    double pDrop; //Pressure drop per length unit
    pDrop = -24 / Re / H * rho * Uref * Uref / 2;

    double dy;  //Length of cell
    dy = H / (NJ - 1);

    double CFL;
    CFL = Uref * dt / dy;
    std::cout << "CFL " << CFL << std::endl;

//********************** Create mesh **************************************
    std::cout << "Create mesh..." << std::endl;

    std::vector<double> Y(NJ);  //Y-coordinates of mesh nodes

    for (size_t j = 0; j != NJ; j++)
        Y[j] = j * dy;

//********************** Allocate & Initiate All Arrays  **********************************
    std::cout << "Initiate Fields..." << std::endl;

    std::vector<double> U(NJ); //X component of Velocity

    U[0] = 0 ;
    U[NJ - 1] = Uw;

    for (size_t j = 1; j != NJ - 1; j++)
        U[j] = Uref;

//********************** CalCFLlate Pressure Gradient ******************************
    std::cout << "Solving equations..." << std::endl;

    Solve(Y, U, dy, dt, rho, nu, nIter, pDrop);

//********************** Output Fields *****************************************
    std::ofstream  outData("out.dat");

    std::cout << "Output fields in file " << std::endl;

    if (outData.is_open()) {
        for(size_t j = 0; j != NJ; j++){
            outData << Y[j] << " " << U[j] << " " << pDrop / nu / 2 * ((Y[j] - 0.5 * H) * (Y[j] - 0.5 * H) - H * H / 4) << "\n";

        }
        outData.close();
    }
    else {
        std::cout << "Invalid output file" << std::endl;
        return 0;
    }


    return 0;
}
