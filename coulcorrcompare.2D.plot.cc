#include <iostream>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"

using namespace std;

const double alpha = 1.2;
const double Ro = 5.7;
const double Rs = 5.0;
const double Rl = 6.8;
const double lambda = 0.8;
const double betaT = 0.95;

void PrintHists(TCanvas *c1, TH2D* hCfull, TH2D* hCappr, TH2D* hCdiff, const char* tag)
{
  hCfull->SetTitle(Form("#alpha=%.1f, R_{o}=%.1f, R_{s}=%.1f, R_{l}=%.1f, #lambda=%.1f, #beta_{T}=%.1f",alpha,Ro,Rs,Rl,lambda,betaT));
  hCfull->Draw("COLZ");
  c1->Print(Form("coulcorrcompare.2D.Cfull.%s.png",tag));
  
  hCappr->SetTitle(Form("#alpha=%.1f, R_{o}=%.1f, R_{s}=%.1f, R_{l}=%.1f, #lambda=%.1f, #beta_{T}=%.1f",alpha,Ro,Rs,Rl,lambda,betaT));
  hCappr->Draw("COLZ");
  c1->Print(Form("coulcorrcompare.2D.Capprox.%s.png",tag));
  
  hCdiff->SetTitle(Form("#alpha=%.1f, R_{o}=%.1f, R_{s}=%.1f, R_{l}=%.1f, #lambda=%.1f, #beta_{T}=%.1f",alpha,Ro,Rs,Rl,lambda,betaT));
  hCdiff->Draw("COLZ");
  hCdiff->GetZaxis()->SetRangeUser(2.e-8,8.e-1);
  c1->SetLogz(1);
  c1->Print(Form("coulcorrcompare.2D.Cdiff.%s.png",tag));
  return; 
}

int main()
{
  TFile* f = new TFile("coulcorrcompare.2D.root");

  TH2D* hCfullqOqL = (TH2D*)f->Get("Cfull2DqOqL");
  TH2D* hCapprqOqL = (TH2D*)f->Get("Cappr2DqOqL");
  TH2D* hCdiffqOqL = (TH2D*)f->Get("Cdiff2DqOqL");
  
  TH2D* hCfullqOqS = (TH2D*)f->Get("Cfull2DqOqS");
  TH2D* hCapprqOqS = (TH2D*)f->Get("Cappr2DqOqS");
  TH2D* hCdiffqOqS = (TH2D*)f->Get("Cdiff2DqOqS");
  
  TH2D* hCfullqSqL = (TH2D*)f->Get("Cfull2DqSqL");
  TH2D* hCapprqSqL = (TH2D*)f->Get("Cappr2DqSqL");
  TH2D* hCdiffqSqL = (TH2D*)f->Get("Cdiff2DqSqL");
  
  TH2D* hCfullqTqL = (TH2D*)f->Get("Cfull2DqTqL");
  TH2D* hCapprqTqL = (TH2D*)f->Get("Cappr2DqTqL");
  TH2D* hCdiffqTqL = (TH2D*)f->Get("Cdiff2DqTqL");
  
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1", "Full", 800, 600);
  PrintHists(c1, hCfullqOqL, hCapprqOqL, hCdiffqOqL, "qOqL");
  PrintHists(c1, hCfullqOqS, hCapprqOqS, hCdiffqOqS, "qOqS");
  PrintHists(c1, hCfullqSqL, hCapprqSqL, hCdiffqSqL, "qSqL");
  PrintHists(c1, hCfullqTqL, hCapprqTqL, hCdiffqTqL, "qTqL");

  return 0;
}
