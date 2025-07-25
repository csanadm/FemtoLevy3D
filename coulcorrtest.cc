#include <iostream>
#include "Levy3D_CoulCalc.h"

using namespace std;

int main()
{
  // Create an instance of the CoulCorrCalc class
  Levy3D_CoulCalc* cccinstance = new Levy3D_CoulCalc();
	cccinstance->SetParticleMass(0.13957039);
  // Declare parameter variables
  double alpha = 1.2;
  double Ro = 5.3;
  double Rs = 5.1;
  double Rl = 5.8;
  double lambda = 0.8;
  // Start loop for calculating the correlation function
  for(double Q=0.001; Q<0.2; Q+=0.001)
  {
    // Full correlation function, including the Coulomb effect, with the specified lambda value
    //double Full3DCorrFuncValueDiag = cccinstance->Full3DCorrFuncValue(alpha, Ro, Rs, Rl, lambda, Q, Q, Q);
    double Full3DCorrFuncValueOut  = cccinstance->Full3DCorrFuncValue(alpha, Ro, Rs, Rl, lambda, Q, 0, 0);
    double Full3DCorrFuncValueSide = cccinstance->Full3DCorrFuncValue(alpha, Ro, Rs, Rl, lambda, 0, Q, 0);
    double Full3DCorrFuncValueLong = cccinstance->Full3DCorrFuncValue(alpha, Ro, Rs, Rl, lambda, 0, 0, Q);
    // Printout
    cout << Q << "\t" << Full3DCorrFuncValueOut << "\t" << Full3DCorrFuncValueSide << "\t" << Full3DCorrFuncValueLong << endl;
  }
  // Delete CoulCorrCalc instance and return
  delete cccinstance;
  return 0;
}

