#include <iostream>
#include <cmath>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
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
  double betaT = 0.95;

  const int Nbins = 100;
  double qmax = 0.10;

  TH2D* hCfull = new TH2D("Cfull2D", "Full C(q_{T},q_{L});q_{T} [GeV];q_{L} [GeV]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCapprox = new TH2D("Capprox2D", "Approximative C(q_{T},q_{L});q_{T} [GeV];q_{L} [GeV]", Nbins, 0, qmax, Nbins, 0, qmax);
  TH2D* hCdiff = new TH2D("Cdiff2D", "Difference C_{full}-C_{approx};q_{T} [GeV];q_{L} [GeV]", Nbins, 0, qmax, Nbins, 0, qmax);

  for (int iT = 1; iT <= Nbins; ++iT)
  {
    double qT = hCfull->GetXaxis()->GetBinCenter(iT);
    for (int iL = 1; iL <= Nbins; ++iL)
    {
      double qL = hCfull->GetYaxis()->GetBinCenter(iL);

      // Map to 3D qOut = qT, qSide = 0
      double qout = qT/sqrt(2.);
      double qside = qT/sqrt(2.);
      double qlong = qL;

      double QoPCMS = qout * sqrt(1 - betaT * betaT);
      double RoPCMS = Ro / sqrt(1 - betaT * betaT);
      double RPCMS = sqrt((RoPCMS * RoPCMS + Rs * Rs + Rl * Rl) / 3.0);

      double Cfull = ccc->Full3DCorrFuncValue(alpha, RoPCMS, Rs, Rl, lambda, QoPCMS, qside, qlong);
      double K2 = ccc->Full3DCoulCorrValue(alpha, RPCMS, RPCMS, RPCMS, 1.0, qout, qside, qlong);
      double qRq = qout * qout * Ro * Ro + qside * qside * Rs * Rs + qlong * qlong * Rl * Rl;
      double Capprox = 1 - lambda + lambda * K2 * (1.0 + exp(-pow(qRq/(hbarcgev*hbarcgev), alpha/2.0)));

      hCfull->SetBinContent(iT, iL, Cfull);
      hCapprox->SetBinContent(iT, iL, Capprox);
      hCdiff->SetBinContent(iT, iL, (Capprox-Cfull)/Cfull);
    }
    cerr << ".";
  }
  cerr << endl;

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1", "Full", 800, 600);
  hCfull->Draw("COLZ");
  c1->Print("coulcorrcompare.2D.Cfull.png");
  hCapprox->Draw("COLZ");
  c1->Print("coulcorrcompare.2D.Capprox.png");
  hCdiff->Draw("COLZ");
  c1->SetLogz(1);
  c1->Print("coulcorrcompare.2D.Cdiff.png");

  return 0;
}
