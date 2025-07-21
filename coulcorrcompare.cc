#include <iostream>
#include "Levy3D_CoulCalc.h"
//#include "CoulCorrCalc.h"

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
	double betaT = 0.90;
  // Start loop for calculating the correlation function
  for(double Qo=0.001; Qo<0.10; Qo+=0.01)
    for(double Qs=0.001; Qs<0.10; Qs+=0.01)
      for(double Ql=0.001; Ql<0.10; Ql+=0.01)
      {
				cerr << Qo << "\t" << Qs << "\t" << Ql; 
      	// Full corr. func. value in 3D
        double Full3DCorrFuncValue = cccinstance->Full3DCorrFuncValue(alpha, Ro, Rs, Rl, lambda, Qo, Qs, Ql);
      	
      	// Approximating a spherical source for the Coulomb correction
      	double Qinv = sqrt((1.0-betaT*betaT)*Qo*Qo + Qs*Qs + Ql*Ql);
        double RPCMS = sqrt((Ro*Ro/(1-betaT*betaT) + Rs*Rs + Rl*Rl)/3);
        //double coulcorrvalue = cccinstanceold->CoulCorrValue(alpha, RPCMS, Qinv);
        double coulcorrvalue = cccinstance->Full3DCoulCorrValue(alpha, RPCMS, RPCMS, RPCMS, 1.0, Qinv/sqrt(3.), Qinv/sqrt(3.), Qinv/sqrt(3.));
        double purecorrfunc = 1.0 + exp(-pow( (Qo*Qo*Ro*Ro+Qs*Qs*Rs*Rs+Ql*Ql*Rl*Rl)/HBARC/HBARC, alpha/2.0));
        double Cqapprox = 1-lambda+lambda*(coulcorrvalue*purecorrfunc);
      	cerr << "\t" << Full3DCorrFuncValue << "\t" << Cqapprox << endl;
      	//cerr << Qo << "\t" << Qs << "\t" << Ql << "\t" << Full3DCorrFuncValue << "\t" << Cqapprox << endl;
      }
  // Delete CoulCorrCalc instance and return
  delete cccinstance;
  return 0;
}

