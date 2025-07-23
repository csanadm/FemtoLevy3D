#include <iostream>
#include <cmath>
#include "TH2D.h"
#include "TFile.h"
#include "Levy3D_CoulCalc.h"

using namespace std;

const double hbarcgev = 0.1973269631;

const double alpha = 1.2;
const double Ro = 5.7;
const double Rs = 5.0;
const double Rl = 6.8;
const double lambda = 0.8;
const double betaT = 0.95;

Levy3D_CoulCalc* ccc;

void FillHists(const int i1, const int i2, const double qout, const double qside, const double qlong, TH2D* hCfull, TH2D* hCappr, TH2D* hCdiff)
{
  double QoPCMS = qout * sqrt(1 - betaT * betaT);
  double RoPCMS = Ro / sqrt(1 - betaT * betaT);
  double RPCMS = sqrt((RoPCMS * RoPCMS + Rs * Rs + Rl * Rl) / 3.0);
  
  double Cfull = ccc->Full3DCorrFuncValue(alpha, RoPCMS, Rs, Rl, lambda, QoPCMS, qside, qlong);
  double K2 = ccc->Full3DCoulCorrValue(alpha, RPCMS, RPCMS, RPCMS, 1.0, QoPCMS, qside, qlong);
  double qRq = qout * qout * Ro * Ro + qside * qside * Rs * Rs + qlong * qlong * Rl * Rl;
  double Capprox = 1 - lambda + lambda * K2 * (1.0 + exp(-pow(qRq/(hbarcgev*hbarcgev), alpha/2.0)));
  
  hCfull->SetBinContent(i1, i2, Cfull);
  hCappr->SetBinContent(i1, i2, Capprox);
  hCdiff->SetBinContent(i1, i2, abs(Capprox-Cfull));
  return ;
}

int main()
{
  TFile* f = new TFile("coulcorrcompare.2D.root", "RECREATE");
  
  ccc = new Levy3D_CoulCalc();
  ccc->SetParticleMass(0.13957039);

  const int Nbins = 100;
  double qmax = 0.10;

  TH2D* hCfullqOqL = new TH2D("Cfull2DqOqL",            "Full C(q_{O},q_{L});q_{O} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCapprqOqL = new TH2D("Cappr2DqOqL",   "Approximative C(q_{O},q_{L});q_{O} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCdiffqOqL = new TH2D("Cdiff2DqOqL", "Difference C_{full}-C_{approx};q_{O} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  
  TH2D* hCfullqOqS = new TH2D("Cfull2DqOqS",            "Full C(q_{O},q_{S});q_{O} [GeV/c];q_{S} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCapprqOqS = new TH2D("Cappr2DqOqS",   "Approximative C(q_{O},q_{S});q_{O} [GeV/c];q_{S} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCdiffqOqS = new TH2D("Cdiff2DqOqS", "Difference C_{full}-C_{approx};q_{O} [GeV/c];q_{S} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  
  TH2D* hCfullqSqL = new TH2D("Cfull2DqSqL",            "Full C(q_{S},q_{L});q_{S} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCapprqSqL = new TH2D("Cappr2DqSqL",   "Approximative C(q_{S},q_{L});q_{S} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCdiffqSqL = new TH2D("Cdiff2DqSqL", "Difference C_{full}-C_{approx};q_{S} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  
  TH2D* hCfullqTqL = new TH2D("Cfull2DqTqL",            "Full C(q_{T},q_{L});q_{T} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCapprqTqL = new TH2D("Cappr2DqTqL",   "Approximative C(q_{T},q_{L});q_{T} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCdiffqTqL = new TH2D("Cdiff2DqTqL", "Difference C_{full}-C_{approx};q_{T} [GeV/c];q_{L} [GeV/c]", Nbins, 0, qmax, Nbins, 0, qmax);

  for (int i1 = 1; i1 <= Nbins; ++i1)
  {
    for (int i2 = 1; i2 <= Nbins; ++i2)
    {
      double q1 = (i1-0.5)*(qmax/Nbins);
      double q2 = (i2-0.5)*(qmax/Nbins);
      double q3 = q1/sqrt(2.);
      
      FillHists(i1, i2, q1,  0, q2, hCfullqOqL, hCapprqOqL, hCdiffqOqL);
      FillHists(i1, i2, q1, q2,  0, hCfullqOqS, hCapprqOqS, hCdiffqOqS);
      FillHists(i1, i2, 0,  q1, q2, hCfullqSqL, hCapprqSqL, hCdiffqSqL);
      FillHists(i1, i2, q3, q3, q2, hCfullqTqL, hCapprqTqL, hCdiffqTqL);
    }
    cerr << ".";
  }
  cerr << endl;

  hCfullqOqL->Write(); 
  hCapprqOqL->Write(); 
  hCdiffqOqL->Write();

  hCfullqOqS->Write();  
  hCapprqOqS->Write();  
  hCdiffqOqS->Write(); 

  hCfullqSqL->Write();  
  hCapprqSqL->Write();  
  hCdiffqSqL->Write();

  hCfullqTqL->Write();  
  hCapprqTqL->Write();  
  hCdiffqTqL->Write();
  
  f->Write();
  f->Close();

  return 0;
}
