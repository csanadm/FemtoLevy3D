#include <iostream>
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TPaletteAxis.h"

void Print2x2(const char* what,
              TH2D* h1, TH2D* h2,
              TH2D* h3, TH2D* h4)
{
  TCanvas* c = new TCanvas(Form("c_%s",what),what,800,800);

  // Define pads manually (tight, no spacing)
  double margin = 0.17;
  double midX = 0.50;
  double midY = 0.50;
  double maxX = 1.00;
  double maxY = 1.00;

  TPad* pad1 = new TPad("pad1","",0.00, midY, midX, maxY); // top-left
  TPad* pad2 = new TPad("pad2","",midX, midY, maxX, maxY); // top-right
  TPad* pad3 = new TPad("pad3","",0.00, 0.00, midX, midY); // bottom-left
  TPad* pad4 = new TPad("pad4","",midX, 0.00, maxX, midY); // bottom-right

  pad1->Draw(); pad2->Draw(); pad3->Draw(); pad4->Draw();

  TH2D* hists[4] = {h1,h2,h3,h4};
  TPad* pads[4]  = {pad1,pad2,pad3,pad4};

  for (int i=0;i<4;++i) {
    pads[i]->cd();
    pads[i]->SetRightMargin((i%2==1) ? margin : 0.01);
    pads[i]->SetLeftMargin((i%2==0) ? margin : 0.01);
    pads[i]->SetBottomMargin((i>=2) ? margin : 0.01);
    pads[i]->SetTopMargin((i<2) ? margin : 0.01);

    if (strcmp(what,"Cdiff")==0) {
      pads[i]->SetLogz(1);
      hists[i]->GetZaxis()->SetRangeUser(2.e-5,5.e-2);
    }

    double qmax = 0.095;
    hists[i]->GetXaxis()->SetRangeUser(0,qmax);
    hists[i]->GetYaxis()->SetRangeUser(0,qmax);

    // Hide inner axes
    if (i%2==1) hists[i]->GetYaxis()->SetLabelSize(0);
    if (i<2)    hists[i]->GetXaxis()->SetLabelSize(0);

    hists[i]->SetTitle("");
    // draw hist
    if (i%2==1) {
      hists[i]->Draw("COLZ");
      gPad->Modified();
      gPad->Update();
      TPaletteAxis *palette = (TPaletteAxis*)hists[i]->GetListOfFunctions()->FindObject("palette");
      palette->SetX1NDC(0.87);
      palette->SetX2NDC(0.92);
      palette->SetY1NDC(0.01+(i-1)*0.077);
      palette->SetY2NDC(0.99+(i-3)*0.077);
      gPad->Modified();
      gPad->Update();
    }
    else        hists[i]->Draw("COL");
    
    // now get pad coordinates
    double xlow = 0;
    double xup  = qmax;
    double ylow = 0;
    double yup  = qmax;
    
    // add top axes
    TGaxis *topAxis = new TGaxis(xlow, yup, xup, yup, xlow, xup, 510, "-U");
    if (i < 2) {
      topAxis->SetLabelSize(hists[i]->GetXaxis()->GetLabelSize());
      topAxis->SetLabelFont(hists[i]->GetXaxis()->GetLabelFont());
      topAxis->SetTitle(hists[i]->GetXaxis()->GetTitle());
      topAxis->SetTitleSize(hists[i]->GetXaxis()->GetTitleSize());
      topAxis->SetTitleOffset(0.4);
    }
    topAxis->Draw();
    
    // add right axes
    TGaxis *rightAxis = new TGaxis(xup, ylow, xup, yup, ylow, yup, 510, "+U");
    if (i % 2 == 1) {
      rightAxis->SetLabelSize(hists[i]->GetYaxis()->GetLabelSize());
      rightAxis->SetLabelFont(hists[i]->GetYaxis()->GetLabelFont());
      rightAxis->SetTitle(hists[i]->GetYaxis()->GetTitle());
      rightAxis->SetTitleSize(hists[i]->GetYaxis()->GetTitleSize());
      rightAxis->SetTitleOffset(0.5);
    }
    rightAxis->Draw();
  }

  c->Print(Form("coulcorrcompare.2D.%s.pdf",what));
}

int main()
{
  TFile* f = new TFile("coulcorrcompare.2D.root");

  // full
  TH2D* hCfullqOqL = (TH2D*)f->Get("Cfull2DqOqL");
  TH2D* hCfullqOqS = (TH2D*)f->Get("Cfull2DqOqS");
  TH2D* hCfullqSqL = (TH2D*)f->Get("Cfull2DqSqL");
  TH2D* hCfullqTqL = (TH2D*)f->Get("Cfull2DqTqL");

  // approx
  TH2D* hCapprqOqL = (TH2D*)f->Get("Cappr2DqOqL");
  TH2D* hCapprqOqS = (TH2D*)f->Get("Cappr2DqOqS");
  TH2D* hCapprqSqL = (TH2D*)f->Get("Cappr2DqSqL");
  TH2D* hCapprqTqL = (TH2D*)f->Get("Cappr2DqTqL");

  // diff
  TH2D* hCdiffqOqL = (TH2D*)f->Get("Cdiff2DqOqL");
  TH2D* hCdiffqOqS = (TH2D*)f->Get("Cdiff2DqOqS");
  TH2D* hCdiffqSqL = (TH2D*)f->Get("Cdiff2DqSqL");
  TH2D* hCdiffqTqL = (TH2D*)f->Get("Cdiff2DqTqL");

  gStyle->SetOptStat(0);

  Print2x2("Cfull", hCfullqOqL,hCfullqOqS,hCfullqSqL,hCfullqTqL);
  Print2x2("Cappr", hCapprqOqL,hCapprqOqS,hCapprqSqL,hCapprqTqL);
  Print2x2("Cdiff", hCdiffqOqL,hCdiffqOqS,hCdiffqSqL,hCdiffqTqL);

  return 0;
}
