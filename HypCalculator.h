#ifndef _HypCalculator_h_
#define _HypCalculator_h_

#include "functions.h"

class HypCalculator
{
 public:
  HypCalculator();
  HypCalculator(const HypCalculator& myHypCalculator);
  ~HypCalculator() {}
  
  // In the following, 2F1(a,b,c,z) is the Gaussian hypergeometric function. For details, see the NIST DLMF: https://dlmf.nist.gov/15, especially https://dlmf.nist.gov/15.8
  
  // Initializetion of eta, useful if 2F1 is calculated for the same eta for many z values, but at the same eta (from which a, b, c are derived)
  void initialize_eta(const double _eta);

  // 2F1(i*eta, 1+i*eta, 1, x-i0), where x is real and positive; use after calling initialize_eta() 
  complex<double> Fplus(const double x);
  
  // 2F1(-i*eta, 1-i*eta, 2, x pm i*0), where x is real positive; use after calling initialize_eta(); is_abovecut true for x+i0, false for x-i0
  complex<double> Fminus(const double x, const bool is_abovecut);
  
// private:
  // 2F1(_a,_b,_c,z), using the power series around x=0
  complex<double> Gauss_2F1_series0(const complex<double> _a, const complex<double> _b, const complex<double> _c, const double _x);

  // 2F1(i*eta,1+i*eta,1,x), using the power series around x=0; use after calling initialize_eta() 
  complex<double> Gauss_2F1_series0(const double x);
  
  // 2F1(-i*eta,1-i*eta,2,x), using the power series around x=0; use after calling initialize_eta()  
  complex<double> Gauss_2F1_series0_minus(const double x);
  
  // 2F1(i*eta,1+i*eta;1;x), with the formula for 1-x argument; use after calling initialize_eta(); assuming x-i0, i.e. below branch cut
  complex<double> Gauss_2F1_1minusx_noint(const double x);
  
  // 2F1(-i*eta,1-i*eta,2,x), with the power series in 1-x; use after calling initialize_eta(); is_abovecut true for x+i0, false for x-i0
  complex<double> Gauss_2F1_1minusx_noint_minus(const double x, const bool is_abovecut);
  
  // 2F1(i*eta,1+i*eta,1,x), using the series (incl. l'Hospital) in 1/x; use after calling initialize_eta(); assuming x-i0, i.e. below branch cut
  complex<double> Gauss_2F1_1overx_spec_c1(const double x);
  
  // 2F1(-i*eta,1-i*eta,2,x), using the series (incl. l'Hospital) in 1/x; use after calling initialize_eta(); is_abovecut true for x+i0, false for x-i0
  complex<double> Gauss_2F1_1overx_spec_c2(const double x, const bool is_abovecut);

  complex<double> gplus;

// private:
  
  // Internal constants and variables
  bool is_eta_initialized;

  complex<double> sinabc;
  complex<double> sinab;
  complex<double> Gamma_c;
  complex<double> Gamma_ab1c;
  complex<double> Gamma_ca;
  complex<double> Gamma_cb;
  complex<double> Gamma_a;
  complex<double> Gamma_b;
  complex<double> Gamma_c1ab;
  complex<double> A;
  complex<double> B;
  complex<double> C;

  complex<double> sinabc_minus;
  complex<double> sinab_minus;
  complex<double> Gamma_c_minus;
  complex<double> Gamma_ab1c_minus;
  complex<double> Gamma_ca_minus;
  complex<double> Gamma_cb_minus;
  complex<double> Gamma_a_minus;
  complex<double> Gamma_b_minus;
  complex<double> Gamma_c1ab_minus;
  complex<double> A_minus;
  complex<double> B_minus;
  complex<double> C_minus;

  double eta;
  double pieta_cth_pieta;
  complex<double> sinhpieta_pieta;
  complex<double> tilde_s0eta;
  complex<double> tilde_s0eta_c2;
  complex<double> tilde_nu0;
};

#endif // _HypCalculator_h_
