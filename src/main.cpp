
#include<iostream>
#include<fstream>

#include <vector>

#include "solver.h"

int main(void) {

  double H;     //Heigh of channel
  H = 1;

  double rho;   //Fluid Density
  rho = 1;

  double nu;    //Fluid viscosity
  nu = 1e-3;

  double Uw;    //Velocity of upper wall
  Uw = 0;

  double Uref;    //Velocity of upper wall
  Uref = 1;

  double Re;
  Re = Uref * H * rho / nu;

  std::cout << "Re " << Re << std::endl;

  double pDrop; //Pressure drop per length unit

  pDrop = -24 / Re / H * rho * Uref * Uref / 2;

  size_t NJ;    //Number of nodes along j-direction
  NJ = 201;

  double dy;
  dy = H / (NJ - 1);

  double dt;
  dt = 0.005;

  double CFL;
  CFL = Uref * dt / dy;
  std::cout << "CFL " << CFL << std::endl;

  size_t nIter; //Max number of iteration
  nIter = 1000000;
  //********************** Read Input Data **************************************
/*
  if (inputData.is_open()) {

    inputData.close();
  }
  else {
    std::cout << "Invalid input file" << std::endl;
    return 0;
  }*/

  //********************** Create mesh **************************************
  std::cout << "Create mesh..." << std::endl;

  std::vector<double> Y(NJ);

  for (size_t j = 0; j != NJ; j++)
          Y[j] = j * dy;

  //********************** Allocate & Initiate All Arrays  **********************************
  std::cout << "Initiate Fields..." << std::endl;

  std::vector<double> U(NJ);

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
