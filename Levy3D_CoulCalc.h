#ifndef _Levy3D_CoulCalc_h_
#define _Levy3D_CoulCalc_h_

//#include <algorithm>
//#include <functional>
//#include <iostream>
//#include <cmath>
//#include <limits>
#include "basics.h"
#include "HypCalculator.h"
#include "Levy3D_CoulCalc.h"

class Levy3D_CoulCalc
{
 public:
  Levy3D_CoulCalc();
  ~Levy3D_CoulCalc();
  
  
  // Set integration properties
  void SetIntegrationProperties(int NMaxIter, double epsTolerance);
	
	// Set particle mass (not reduced mass) in GeV/c^2
  void SetParticleMass(const double _Mc2); 
  
  // Get number of function calls in last calculation
  int GetNFuncCalls() { return NFuncCalls; }

  // Calculate full 3D correlation function, Ro, Rs, Rl: HBT radii in fm, Qo, Qs, Ql: momentum differences in PCMS, GeV/c
	double Full3DCorrFuncValue(double alpha, double Ro, double Rs, double Rl, double lambda, double Qo, double Qs, double Ql);
	
  // Calculate full 3D Coulomb correction, Ro, Rs, Rl: HBT radii in fm, Qo, Qs, Ql: momentum differences in PCMS, GeV/c
	double Full3DCoulCorrValue(double alpha, double Ro, double Rs, double Rl, double lambda, double Qo, double Qs, double Ql);

 private:
  // Private functions, explained in the source code
  double f_s_aa_y_int(const double y, const double q, const double Rcc_LCMS, const double alpha, const double betat);
  double f_s(const double q, const double R, const double alpha, const double betat);
  double A_1_s_wo_int(const double x, const double k, const double R, const double alpha, const double betat);
  double A_2_s_wo_int(const double x, const double k, const double R, const double alpha, const double betat);
  double A1_int(const double k, const double R, const double alpha, const double betat);
  double A2_int(const double k, const double R, const double alpha, const double betat);

  inline double f_3d(const double _qout, const double _qside, const double _qlong); // pars must be already initialized!
  inline double f1_3d(const double phi, const double a, const double beta);
  inline double f2_3d(const double phi, const double b, const double y);
  double f1_3d_phiintegrated(const double xi, const double a); // beta = xi^2/(1.-xi)
  double f2_3d_phiintegrated(const double y, const double b);
  double f1_3d_phixi_integrated(const double a);
  double f2_3d_phiy_integrated(const double b);

  double f1_A1_x_integrand(const double x);
  double get_A1();

  double f2_3d_phixi_transformed_integrand(const double xi, const double x, const double atan1x, const double b);
  double f2_3d_phixi_transformed_integrated(const double x, const double atan1x, const double b);
  double f2_A2P_x_integrand(const double x);
  double get_A2P();

  double f2_A2P_xxi_integrand(const double xi, const double x); // x -> x^2, y -> x * tan(xi * atan(1./x)); as in note
  double f2_A2d_1y_integrand(const double y); // y \in [-1,1]
  double get_A2d();

  void init_3d_sourcepars(const double _R2oo, const double _R2ss, const double _R2ll, const double _alpha);

  void set_kvector(const double _kx, const double _ky, const double _kz);
  double get_eta() const { return eta; }

  void calc_qvector_for_A1(const double a, const double beta, const double phi);
  void calc_qvector_for_A2(const double b, const double y, const double phi);

  void calc_eta(const double k); // k: MeV/c, Mc2: in MeV, must be already set by set_particle_mass()
	
  double get_C2Pi_3d_full(const double _alpha, const double _R2oo, const double _R2ss, const double _R2ll, const double Qo, const double Qs, const double Ql);

  // An instance of the hypergeometric 2F1 calculator
  HypCalculator* HypCalculatorInstance;
  
  // Imaginary unit (1i since c++17)
  complex<double> I = complex<double>(0., 1.);
  complex<double> complex_zero = complex<double>(0., 0.);

  double Mc2;
  double eta;
  double Ychi; // = exp(2*pi*eta);
  double Gamow_factor;

  double f_3d_0;
  double f_3d_2k;

  double k_abs;
  double kt;
  double kx;
  double ky;
  double kz;
  double Theta;
  double Phi;
  double CosTheta;
  double SinTheta;
  double CosPhi;
  double SinPhi;

  double R2oo;
  double R2ss;
  double R2ll;
  double alpha;

  double epx_x;  double epy_x;  double epz_x;
  double epx_y;  double epy_y;  double epz_y;
  double epx_z;  double epy_z;  double epz_z;
  
  double qx;
  double qy;
  double qz;

  complex<double> Fm_x2_x2;
  complex<double> Fm_2x_up;
  complex<double> Fm_2x_dw;
  complex<double> Fm_XX;

  // Integration properties
  static const unsigned int NGaussKronrod = 15; // this has to be constant, to be set at compile time
  unsigned int NMaxIter;
  double epsTolerance;
  int NFuncCalls;
};

#endif // _Levy3D_CoulCalc_h_

