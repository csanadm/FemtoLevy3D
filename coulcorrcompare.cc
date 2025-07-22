#include <iostream>
#include <fstream>
#include <cmath>
#include "Levy3D_CoulCalc.h"

using namespace std;

const double hbarcgev = 0.1973269631;

int main()
{
  Levy3D_CoulCalc* ccc = new Levy3D_CoulCalc();
  ccc->SetParticleMass(0.13957039);

  double alpha = 1.2;
  double Ro = 5.7;
  double Rs = 5.0;
  double Rl = 6.8;
  double lambda = 0.8;
  double betaT = 0.9;

  for(double Q=0.001; Q<0.150; Q+=0.001)
  {
    double QoPCMS = Q*sqrt(1-betaT*betaT);
    double RoPCMS = Ro/sqrt(1-betaT*betaT);
    double RPCMS = sqrt((RoPCMS*RoPCMS + Rs*Rs + Rl*Rl)/3);
		
    // Cfull (3D Coulomb)
    double Cfull_qout  = ccc->Full3DCorrFuncValue(alpha, RoPCMS, Rs, Rl, lambda, QoPCMS, 0, 0);
    double Cfull_qside = ccc->Full3DCorrFuncValue(alpha, RoPCMS, Rs, Rl, lambda, 0,      Q, 0);
    double Cfull_qlong = ccc->Full3DCorrFuncValue(alpha, RoPCMS, Rs, Rl, lambda, 0,      0, Q);
		
    // Capprox (spherical Coulomb)
    double Capprox_qout  = 1 - lambda + lambda * ccc->Full3DCoulCorrValue(alpha, RPCMS, RPCMS, RPCMS, 1.0, QoPCMS, 0, 0) * (1.0 + exp( -pow((Q*Q*Ro*Ro)/hbarcgev/hbarcgev, alpha/2.0) ));
    double Capprox_qside = 1 - lambda + lambda * ccc->Full3DCoulCorrValue(alpha, RPCMS, RPCMS, RPCMS, 1.0, 0,      Q, 0) * (1.0 + exp( -pow((Q*Q*Rs*Rs)/hbarcgev/hbarcgev, alpha/2.0) ));
    double Capprox_qlong = 1 - lambda + lambda * ccc->Full3DCoulCorrValue(alpha, RPCMS, RPCMS, RPCMS, 1.0, 0,      0, Q) * (1.0 + exp( -pow((Q*Q*Rl*Rl)/hbarcgev/hbarcgev, alpha/2.0) ));

    cout << Q << "\t"
         << Cfull_qout << "\t" << Capprox_qout << "\t"
         << Cfull_qside << "\t" << Capprox_qside << "\t"
         << Cfull_qlong << "\t" << Capprox_qlong << "\n";
    cerr << ".";
  }
  cerr << endl;

  delete ccc;
  return 0;
}
