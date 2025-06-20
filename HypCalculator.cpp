#include "HypCalculator.h"
const complex<double> I = complex<double>(0., 1.);
const double EPSILON = 1e-12;
const double EPSILON_ABC = 1e-6;
const double euler_gamma = 0.5772156649;

HypCalculator::HypCalculator()
{
  is_eta_initialized = false;
  eta = 0.;

  gplus       = complex<double>(0., 0.);
  sinabc      = complex<double>(0., 0.);
  sinab       = complex<double>(0., 0.);
  Gamma_c     = complex<double>(0., 0.);
  Gamma_ab1c  = complex<double>(0., 0.);
  Gamma_ca    = complex<double>(0., 0.);
  Gamma_cb    = complex<double>(0., 0.);
  Gamma_a     = complex<double>(0., 0.);
  Gamma_b     = complex<double>(0., 0.);
  Gamma_c1ab  = complex<double>(0., 0.);

  A = complex<double>(0., 0.);
  B = complex<double>(0., 0.);
  C = complex<double>(0., 0.);

  sinabc_minus     = complex<double>(0., 0.); 
  sinab_minus      = complex<double>(0., 0.);
  Gamma_c_minus    = complex<double>(0., 0.);
  Gamma_ab1c_minus = complex<double>(0., 0.);
  Gamma_ca_minus   = complex<double>(0., 0.);
  Gamma_cb_minus   = complex<double>(0., 0.);
  Gamma_a_minus    = complex<double>(0., 0.);
  Gamma_b_minus    = complex<double>(0., 0.);
  Gamma_c1ab_minus = complex<double>(0., 0.);
  A_minus = complex<double>(0., 0.); 
  B_minus = complex<double>(0., 0.);
  C_minus = complex<double>(0., 0.);

  pieta_cth_pieta = 0.;
  tilde_s0eta = complex<double>(0., 0.);
  tilde_s0eta_c2 = complex<double>(0., 0.);
  tilde_nu0 = complex<double>(0., 0.);
  sinhpieta_pieta = complex<double>(0., 0.);

};

HypCalculator::HypCalculator(const HypCalculator& myHypCalculator)
{
  is_eta_initialized = myHypCalculator.is_eta_initialized;
}

complex<double> HypCalculator::Gauss_2F1_series0(const complex<double> _a, const complex<double> _b, const complex<double> _c, const double x)
{
  complex<double> result(0., 0);
  if(fabs(x) < 0.9999)
  {
    complex<double> n(0., 0.);
    complex<double> term(1., 0);
    result = term;
    while(abs(term) > EPSILON)
    {
      term = term * x / (n + 1.) * (_a + n) * (_b + n) / (_c + n);
      result += term;
      n += 1.;
    }
  }
  return result;
}

void HypCalculator::initialize_eta(const double _eta)
{
  eta = _eta;
  
  A = complex<double>(0., eta);
  B = complex<double>(1., eta);
  C = complex<double>(1., 0.);
  A_minus = complex<double>(0., -eta);
  B_minus = complex<double>(1., -eta);
  C_minus = complex<double>(2., 0.);
  
  gplus      = Gamma(0.5-I*eta) / Gamma(1.-I*eta) / sqrt(M_PI) * pow(2.,-1.-2.*I*eta);
  sinabc     = sin(M_PI * (C-A-B));  
  sinab      = sin(M_PI * (A-B));  
  Gamma_c    = Gamma(C);
  Gamma_ab1c = Gamma(A+B+1.-C); 
  Gamma_ca   = Gamma(C-A); 
  Gamma_cb   = Gamma(C-B); 
  Gamma_a    = Gamma(A); 
  Gamma_b    = Gamma(B); 
  Gamma_c1ab = Gamma(C+1.-A-B); 
  sinabc_minus     = sin(M_PI * (C_minus-A_minus-B_minus));  
  sinab_minus      = sin(M_PI * (A_minus-B_minus));  
  Gamma_c_minus    = Gamma(C_minus);
  Gamma_ab1c_minus = Gamma(A_minus+B_minus+1.-C_minus); 
  Gamma_ca_minus   = Gamma(C_minus-A_minus); 
  Gamma_cb_minus   = Gamma(C_minus-B_minus); 
  Gamma_a_minus    = Gamma(A_minus); 
  Gamma_b_minus    = Gamma(B_minus); 
  Gamma_c1ab_minus = Gamma(C_minus+1.-A_minus-B_minus); 
  
  if(eta == 0.)
  {
    pieta_cth_pieta = 1.;
    sinhpieta_pieta = 1.;
  }
  else
  {
    pieta_cth_pieta = M_PI * eta / tanh(M_PI * eta);
    sinhpieta_pieta = sinh(M_PI * eta) / (M_PI * eta);
  }
  
  complex<double> s0(2. * eta * euler_gamma - eta , -1. * pieta_cth_pieta);
  tilde_s0eta = s0;
  complex<double> ieta1(1., eta);
  tilde_s0eta += 2. * eta * digamma(ieta1);

  complex<double> s0_c2(1. - pieta_cth_pieta, 2. * eta * euler_gamma - eta);
  tilde_s0eta_c2 = s0_c2;
  complex<double> mieta1(1., -eta);
  tilde_s0eta_c2 += 2. * I * eta * digamma(mieta1);

  tilde_nu0 = -1. * eta;
  
  is_eta_initialized = true;
}


complex<double> HypCalculator::Gauss_2F1_series0(const double x)
{
  complex<double> result(-9999., 0);
  if(fabs(x) < 0.9999)
  {
    complex<double> n(0., 0.);
    complex<double> term(1., 0);
    result = term;
    while(abs(term) > EPSILON)
    {
      term = term * x / (n + 1.) * (A + n) * (B + n) / (C + n);
      result += term;
      n += 1. ;
    }
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_series0_minus(const double x)
{
  complex<double> result(-9999., 0);
  if(fabs(x) < 0.9999)
  {
    complex<double> n(0., 0.);
    complex<double> term(1., 0);
    result = term;
    while(abs(term) > EPSILON)
    {
      term = term * x / (n + 1.) * (A_minus + n) * (B_minus + n) / (C_minus + n);
      result += term;
      n += 1. ;
    }
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1minusx_noint(const double x)
{
  complex<double> result(0., 0);
  if((fabs(1. - x) < 0.9999) && (x != 1.))
    if(abs(sinabc) > EPSILON_ABC)
      result = M_PI / sinabc * Gamma_c * ( Gauss_2F1_series0(A, B, A+B+1.-C, 1.-x) / Gamma_ab1c / Gamma_ca / Gamma_cb
                                          -Gauss_2F1_series0(C-A, C-B, C+1.-A-B, 1.-x) / Gamma_a / Gamma_b / Gamma_c1ab / pow(1.-x, A+B-C));
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1minusx_noint_minus(const double x, const bool is_abovecut)
{
  complex<double> result(0., 0);
  complex<double> aux = (is_abovecut ? complex<double>(x-1., 0.) : complex<double>(1.-x, 0.) );
  if(fabs(1. - x) < 0.9999)
  {
    if(abs(sinabc_minus) > EPSILON_ABC)
      result = M_PI / sinabc_minus * Gamma_c_minus * ( Gauss_2F1_series0(A_minus, B_minus, A_minus+B_minus+1.-C_minus, 1.-x) / Gamma_ab1c_minus / Gamma_ca_minus
                                         / Gamma_cb_minus
                                          -Gauss_2F1_series0(C_minus-A_minus, C_minus-B_minus, C_minus+1.-A_minus-B_minus, 1.-x) 
                        / Gamma_a_minus / Gamma_b_minus / Gamma_c1ab_minus / pow( ( is_abovecut ? -1. : 1. ) * aux , A_minus+B_minus-C_minus));
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1overx_spec_c1(const double x)
{
  complex<double> result(0., 0.);
  complex<double> aux(-x, 0.);
  if(fabs(1./x) < 0.9999)
  {
    complex<double> s = - eta * log(aux) + tilde_s0eta; // ide beírtam egy mínuszt
    complex<double> nu = tilde_nu0 / x;
    double n = 0.;
    complex<double> sum = nu * s;
    while(abs(nu * s) > EPSILON)
    {
      nu = nu * (1. + I * eta / (n + 1.)) * (1. + (I * eta - 1.)/(n + 2.)) / x;
      s +=  2. * eta /(n + 1. + I * eta) - eta * (2. * n + 3.)/(n * n + 3. * n + 2.);
      sum += (nu * s);
      n += 1.;
    }   
    result = sinhpieta_pieta * pow(aux, -I * eta) * (1. + sum);
  }
  return result;
}

complex<double> HypCalculator::Gauss_2F1_1overx_spec_c2(const double x, const bool is_abovecut)
{
  complex<double> result(0., 0.);
  complex<double> aux = (is_abovecut ? complex<double>(x, 0.) : complex<double>(-x, 0.) );
  if(fabs(1./x) < 0.9999)
  {
    complex<double> s = - I * eta * log( ( is_abovecut ? -1. : 1. ) * aux) + tilde_s0eta_c2;
    complex<double> nu = 1. / x;
    double n = 0.;
    complex<double> sum = nu * s;
    while(abs(nu * s) > EPSILON)
    {
      nu = nu * (1. - (I * eta + 1.) / (n + 1.)) * (1. - (I * eta + 1.)/(n + 2.)) / x;
      s +=  I * eta /(n + 1. - I * eta) + I *eta / (n - I * eta) - I * eta * (2. * n + 3.)/(n * n + 3. * n + 2.);
      sum += (nu * s);
      n += 1.;
    }   
    result = sinhpieta_pieta * pow( ( is_abovecut ? -1. : 1. ) * aux, I * eta) * ( 1. / (1. + I*eta) + sum);
  }
  return result;
}

complex<double> HypCalculator::Fplus(const double x) // -> Fplus
{
  complex<double> result(0., 0.);
  if(x >= 0. && is_eta_initialized)
  {
    if(x <= 0.6)
      result = Gauss_2F1_series0(x);
    else if(x <= 1.6)
      result = Gauss_2F1_1minusx_noint(x);
    else
      result = Gauss_2F1_1overx_spec_c1(x);
  }
  return result;

}

complex<double> HypCalculator::Fminus(const double x, const bool is_abovecut) // -> Fminus, 
{
  complex<double> result(0., 0.);
  if(x >= 0. && is_eta_initialized)
  {
    if(x <= 0.6)
      result = Gauss_2F1_series0_minus(x);
    else if(x <= 1.6)
      result = Gauss_2F1_1minusx_noint_minus(x, is_abovecut);
    else
      result = Gauss_2F1_1overx_spec_c2(x, is_abovecut);
  }
  return (1. + I * eta) * result;
}
