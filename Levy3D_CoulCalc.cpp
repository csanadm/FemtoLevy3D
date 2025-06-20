#include "Levy3D_CoulCalc.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>

// Constructor
Levy3D_CoulCalc::Levy3D_CoulCalc()
{
  Mc2 = SQR(139.57039); // to start with the pion mass... better to be set when starting calculations anyway!
  eta = 0.;
  Ychi = 1.;
  Gamow_factor = 1.;

  k_abs = 0.;
  kt = 0.;
  kx = 0.;
  ky = 0.;
  kz = 0.;
  Theta = 0.;
  Phi = 0.;
  CosTheta = 0.;
  SinTheta = 0.;
  CosPhi = 0.;
  SinPhi = 0.;
  
  epx_x = 0.;  epy_x = 0.;  epz_x = 0.;
  epx_y = 0.;  epy_y = 0.;  epz_y = 0.;
  epx_z = 0.;  epy_z = 0.;  epz_z = 0.;
  
  qx = 0.;
  qy = 0.;
  qz = 0.;

  f_3d_0 = 0.;
  f_3d_2k = 0.;

  R2oo = 0.;
  R2ss = 0.;
  R2ll = 0.;
  alpha = 0.;

  Fm_x2_x2 = complex_zero;
  Fm_2x_up = complex_zero;
  Fm_2x_dw = complex_zero;
  Fm_XX = complex_zero;

  HypCalculatorInstance = new HypCalculator();
  NMaxIter = 3;
  epsTolerance = 1e-2;
}

// Destructor
Levy3D_CoulCalc::~Levy3D_CoulCalc()
{
  delete HypCalculatorInstance;
}

// Set integration properties
void Levy3D_CoulCalc::SetIntegrationProperties(int _NMaxIter, double _epsTolerance)
{
  NMaxIter = _NMaxIter;
  epsTolerance = _epsTolerance;
  SetParticleMass(Mass_Pi);
}

void Levy3D_CoulCalc::SetParticleMass(const double _Mc2) // Mc2 is the mass of the particle (not the reduced mass) IN GEV/C^2
{
  Mc2 = _Mc2*1000; // GeV -> MeV conversion
}

// Calculating the Sommerfeld parameter eta
void Levy3D_CoulCalc::calc_eta(const double k)
{
  eta = FINESTRUCTURE_CONSTANT * Mc2 * 0.5 / k;
  Ychi = exp(2.*M_PI * eta);
  Gamow_factor = (2.*M_PI * eta) / (Ychi - 1.);
}

double Levy3D_CoulCalc::f_3d(const double _qout, const double _qside, const double _qlong) // pars must be already initialized!
{
  return exp(-0.5 * pow((_qout * _qout * R2oo + _qside * _qside * R2ss + _qlong * _qlong * R2ll) / SQR(HBARC) , alpha/2.));
}

double Levy3D_CoulCalc::f1_3d(const double phi, const double a, const double beta ) // never called when beta==0 !!!
{
  calc_qvector_for_A1(a, beta, phi);
  return (f_3d(qx, qy, qz) - f_3d_0) / beta;
}

double Levy3D_CoulCalc::f2_3d(const double phi, const double b, const double y)
{
  calc_qvector_for_A2(b, y, phi);
  return (f_3d(qx, qy, qz) - f_3d_2k) / (1.+y);
}

double Levy3D_CoulCalc::f1_3d_phiintegrated(const double xi, const double a)
{
  if(xi == 0)
    return 0.; // except if alpha <= 0.5; DONT GO THERE!!!!
  double error = 0.;
  double beta = xi * xi / (1.-xi);
  auto func = bind(&Levy3D_CoulCalc::f1_3d, this, placeholders::_1, a, beta);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, -M_PI, M_PI, NMaxIter, epsTolerance, &error);
  return result / (a*a + beta*beta) * xi * (2. - xi) / (1.-xi) / (1.-xi);
}

double Levy3D_CoulCalc::f1_3d_phixi_integrated(const double a)
{
  double error = 0.;
  auto func = bind(&Levy3D_CoulCalc::f1_3d_phiintegrated, this, placeholders::_1, a);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, 0., 1., NMaxIter, epsTolerance, &error);
  return result;
}

double Levy3D_CoulCalc::f2_3d_phiintegrated(const double y, const double b)
{
  double error = 0.;
  auto func = bind(&Levy3D_CoulCalc::f2_3d, this, placeholders::_1, b, y);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, -M_PI, M_PI, NMaxIter, epsTolerance, &error);
  return result * (1.-y) / (1. + b*y*y);
}

double Levy3D_CoulCalc::f2_3d_phixi_transformed_integrand(const double xi, const double x, const double atan1x, const double b)
{
  double cosxiatan = cos(xi * atan1x);
  return f2_3d_phiintegrated(x * tan(xi * atan1x), b) / SQR(cosxiatan);
}

double Levy3D_CoulCalc::f2_3d_phixi_transformed_integrated(const double x, const double atan1x, const double b)
{
  double error = 0.;
  auto func = bind(&Levy3D_CoulCalc::f2_3d_phixi_transformed_integrand, this, placeholders::_1, x, atan1x, b);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, -1., 1., NMaxIter, epsTolerance, &error);
  return result;
}

double Levy3D_CoulCalc::f2_A2d_1y_integrand(const double y) // y \in [-1,1]
{
  return f2_3d_phiintegrated(y, 1.);
}

double Levy3D_CoulCalc::get_A2d()
{
  double error = 0.;
  auto func = bind(&Levy3D_CoulCalc::f2_A2d_1y_integrand, this, placeholders::_1);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, -1.0, 1.0, NMaxIter, epsTolerance, &error);
  return 2. * imag(HypCalculatorInstance->gplus) * 1./Gamow_factor * result;
}

double Levy3D_CoulCalc::f1_A1_x_integrand(const double x)
{
  double f1_3d_xx = f1_3d_phixi_integrated(x);
  double f1_3d_1x = f1_3d_phixi_integrated(1./x);
  double f1_3d_11 = f1_3d_phixi_integrated(1.);
  complex<double> gplusc = conj(HypCalculatorInstance->gplus);
  complex<double> xpow = pow(x, 1.+2.*I*eta);
  Fm_x2_x2 = HypCalculatorInstance->Fminus(x*x, false) / SQR(x);
  Fm_2x_up = HypCalculatorInstance->Fminus(1./(x*x), true)  * xpow;
  Fm_2x_dw = HypCalculatorInstance->Fminus(1./(x*x), false) * xpow;
  complex<double> result = pow( 1./(1.+x), 1.+2.*I*eta) * ( 4.*I*gplusc * (1./SQR(x+I) - xpow / SQR(x-I)   ) * f1_3d_11 +        Fm_2x_up * f1_3d_xx + Fm_x2_x2 * f1_3d_1x)
                          +pow( 1./(1.-x), 1.+2.*I*eta) * ( 4.*I*gplusc * (Ychi*xpow/SQR(x+I) + 1./SQR(x-I)) * f1_3d_11 - Ychi * Fm_2x_dw * f1_3d_xx + Fm_x2_x2 * f1_3d_1x);
  return real(result);
}

double Levy3D_CoulCalc::get_A1()
{
  double error = 0.;
  auto func = bind(&Levy3D_CoulCalc::f1_A1_x_integrand, this, placeholders::_1);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, 0.0, 1.0, NMaxIter, epsTolerance, &error);
  return -1./M_PI * result;
}

double Levy3D_CoulCalc::f2_A2P_x_integrand(const double x)
{
  double X = x*x; 
  double atan1x = atan(1./x);
  double f2_3d_xx = f2_3d_phixi_transformed_integrated(x, atan1x, X);
  double f2_3d_1x = f2_3d_phixi_transformed_integrated(x, atan1x, 1./X);
  double f2_3d_11 = f2_3d_phixi_transformed_integrated(x, atan1x, 1.);
  complex<double> gplusc = conj(HypCalculatorInstance->gplus);
  Fm_XX = HypCalculatorInstance->Fminus( 4.*X  / SQR(1. + X), false);
  complex<double> result = 1./(1.-X) * pow( (1.+X)/(1.-X), 2.*I*eta) * (Fm_XX * (Ychi * f2_3d_xx - 1./X * f2_3d_1x) - 4. * gplusc * (Ychi - 1.) * f2_3d_11 / (1.+X));
  return real(result) * 2.*x*x*atan1x;
}

double Levy3D_CoulCalc::get_A2P()
{
  double error = 0.;
  auto func = bind(&Levy3D_CoulCalc::f2_A2P_x_integrand, this, placeholders::_1);
  double result = boost::math::quadrature::gauss_kronrod<double, NGaussKronrod>::integrate(func, 0.0, 1.0, NMaxIter, epsTolerance, &error);
  return 1./M_PI * result;
}

double Levy3D_CoulCalc::f2_A2P_xxi_integrand(const double xi, const double x) // x -> x^2, y -> x * tan(xi * atan(1./x)); as in note
{
  double X = x*x; 
  double dXx = 1.; 
  double atan1x = atan(1./x);
  double xiatan = xi * atan1x;
  double cosxiatan = cos(xiatan);
  double f2_3d_xx = f2_3d_phiintegrated(x * tan(xiatan), X)    / SQR(cosxiatan);;
  double f2_3d_1x = f2_3d_phiintegrated(x * tan(xiatan), 1./X) / SQR(cosxiatan);;
  double f2_3d_11 = f2_3d_phiintegrated(x * tan(xiatan), 1.)   / SQR(cosxiatan);;
  complex<double> gplusc = conj(HypCalculatorInstance->gplus);
  Fm_XX = HypCalculatorInstance->Fminus( 4.*X  / SQR(1. + X), false);
  complex<double> result = 1./(1.-X) * pow( (1.+X)/(1.-X), 2.*I*eta) * (Fm_XX * (Ychi * f2_3d_xx - 1./X * f2_3d_1x) - 4. * gplusc * (Ychi - 1.) * f2_3d_11 / (1.+X));
  return dXx * real(result) * 2.*X*atan1x;
}

double Levy3D_CoulCalc::get_C2Pi_3d_full(const double _alpha, const double _R2oo, const double _R2ss, const double _R2ll, const double Qo, const double Qs, const double Ql)
{
  init_3d_sourcepars(_R2oo, _R2ss, _R2ll, _alpha);
  set_kvector(Qo / 2.0, Qs / 2.0, Ql / 2.0);
  return Gamow_factor * ( f_3d_0 + f_3d_2k +  eta / M_PI * (get_A1() + get_A2P() + get_A2d()) );
}

// This uses the regular radii and GeV units!!!!
double Levy3D_CoulCalc::Full3DCorrFuncValue(double alpha, double Ro, double Rs, double Rl, double lambda, double Qo, double Qs, double Ql)
{
	double Rocc = Ro*pow(2.,1./alpha);
	double Rscc = Rs*pow(2.,1./alpha);
	double Rlcc = Rl*pow(2.,1./alpha);
	return 1. - lambda + lambda * get_C2Pi_3d_full(alpha, Rocc*Rocc, Rscc*Rscc, Rlcc*Rlcc, Qo*1000., Qs*1000, Ql*1000);
}

void Levy3D_CoulCalc::init_3d_sourcepars(const double _R2oo, const double _R2ss, const double _R2ll, const double _alpha)
{
  R2oo = _R2oo;
  R2ss = _R2ss;
  R2ll = _R2ll;
  alpha = _alpha;
  f_3d_0 = f_3d(0., 0., 0.);
}

void Levy3D_CoulCalc::set_kvector(const double _kx, const double _ky, const double _kz)
{
  kx = _kx;
  ky = _ky;
  kz = _kz;
  k_abs = sqrt(SQR(kx) + SQR(ky) + SQR(kz));
  kt = sqrt(SQR(kx) + SQR(ky));
  if(k_abs > 1e-15)
  {
    CosTheta = kz / k_abs;
    SinTheta = kt / k_abs;
    Theta = acos(CosTheta);
    if(kt > 1e-15)
    {
      CosPhi = kx / kt;
      SinPhi = ky / kt;
      Phi = atan2(ky, kx);
    }
  }

  epx_x = -CosTheta * CosPhi;  epy_x = SinPhi;   epz_x = SinTheta * CosPhi;
  epx_y = -CosTheta * SinPhi;  epy_y = -CosPhi;  epz_y = SinTheta * SinPhi;
  epx_z = SinTheta;            epy_z = 0.;       epz_z = CosTheta;

  f_3d_2k = f_3d(2.*kx, 2.*ky, 2.*kz);
  calc_eta(k_abs);
  HypCalculatorInstance->initialize_eta(eta);
}

void Levy3D_CoulCalc::calc_qvector_for_A1(const double a, const double beta, const double phi)
{
  qx = 0.;
  qy = 0.;
  qz = 0.;
  if( (beta >= 0) && ( (fabs(a) > 1e-15) || fabs(beta) > 1e-15) )
  {
    double norm = 2. * k_abs * a * beta / (SQR(a) + SQR(beta));
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    qx = norm * ( a * cosphi * epx_x + a * sinphi * epy_x + beta * epz_x); 
    qy = norm * ( a * cosphi * epx_y + a * sinphi * epy_y + beta * epz_y); 
    qz = norm * ( a * cosphi * epx_z + a * sinphi * epy_z + beta * epz_z); 
  }
}

void Levy3D_CoulCalc::calc_qvector_for_A2(const double b, const double y, const double phi)
{
  qx = 0.;
  qy = 0.;
  qz = 0.;
  if( (b >= 0) && (fabs(y) <= 1.) )
  {
    double norm = k_abs * (1.-y) / (1. + b*y*y);
    double term_xy = sqrt(b) * (1.+y);
    double term_z  = (b*y-1.);
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    qx = norm * ( term_xy * cosphi * epx_x + term_xy * sinphi * epy_x + term_z * epz_x); 
    qy = norm * ( term_xy * cosphi * epx_y + term_xy * sinphi * epy_y + term_z * epz_y); 
    qz = norm * ( term_xy * cosphi * epx_z + term_xy * sinphi * epy_z + term_z * epz_z); 
  }
}
