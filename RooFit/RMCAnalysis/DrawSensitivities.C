#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TMath.h"

using namespace TMath;

void DrawSensitivities(){

  Bool_t UseAllTargets=1;

  // Only External
  //Double_t S  = 1.56168*Power(10,-15);
  //Double_t Ca = 8.1176*Power(10,-16);
  //Double_t Ti = 3.25893*Power(10,-15);
  //Double_t Cr = 1.73621*Power(10,-15);
  //Double_t Fe = 2.81146*Power(10,-15);
  //Double_t Ni = 6.61892*Power(10,-15);
  //Double_t Zn = 5.8104*Power(10,-15);
  //Double_t Ge = 6.11824*Power(10,-15);

  // Include Internal
  //Double_t S  = 1.10766*Power(10,-15); // Z=16
  Double_t S  = 1.00695*Power(10,-15); // Z=16
  Double_t Ca = 7.91461*Power(10,-16); // Z=20
  //Double_t Ti = 2.15718*Power(10,-15); // Z=22
  Double_t Ti = 2.07354*Power(10,-15); // Z=22
  Double_t Cr = 1.68553*Power(10,-15); // Z=24
  Double_t Fe = 2.59838*Power(10,-15); // Z=26
  Double_t Ni = 6.38783*Power(10,-15); // Z=28
  Double_t Zn = 5.59349*Power(10,-15); // Z=30
  //Double_t Ge = 5.40889*Power(10,-15); // Z=32
  Double_t Ge = 5.37784*Power(10,-15); // Z=32

  const Int_t n = 8;
  Double_t Z[8]   = {16,20,22,24,26,28,30,32};
  Double_t ses[8] = {S,Ca,Ti,Cr,Fe,Ni,Zn,Ge};

  //const Int_t n = 4;
  //Double_t Z[n]   = {20,22,26,28};
  //Double_t ses[n] = {Ca,Ti,Fe,Ni};

  TCanvas *c = new TCanvas("canvas","canvas",800,800);
  TGraph * gr = new TGraph(n,Z,ses);
  gr->SetMarkerStyle(33);
  gr->SetMarkerSize(2.7);
  gr->SetMarkerColor(kBlue);

  Double_t xMin=14.5;
  Double_t xMax=34.5;
  Double_t textSize=0.04;

  TLine *line = new TLine(xMin,Ca,xMax,Ca);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);
  TLatex *latexS  = new TLatex(gr->GetX()[0]-gr->GetX()[1]*0.02, gr->GetY()[0]+gr->GetY()[1]*0.15," #font[42]{^{32}S}");
  latexS->SetTextSize(textSize);
  TLatex *latexCa = new TLatex(gr->GetX()[1]-gr->GetX()[1]*0.02, gr->GetY()[1]+gr->GetY()[1]*0.15," #font[42]{^{40}Ca}");
  latexCa->SetTextSize(textSize);
  TLatex *latexTi = new TLatex(gr->GetX()[2]-gr->GetX()[1]*0.02, gr->GetY()[2]+gr->GetY()[1]*0.15," #font[42]{^{48}Ti}");
  latexTi->SetTextSize(textSize);
  TLatex *latexCr = new TLatex(gr->GetX()[3]-gr->GetX()[1]*0.02, gr->GetY()[3]+gr->GetY()[1]*0.15," #font[42]{^{50}Cr}");
  latexCr->SetTextSize(textSize);
  TLatex *latexFe = new TLatex(gr->GetX()[4]-gr->GetX()[1]*0.02, gr->GetY()[4]+gr->GetY()[1]*0.15," #font[42]{^{54}Fe}");
  latexFe->SetTextSize(textSize);
  TLatex *latexNi = new TLatex(gr->GetX()[5]-gr->GetX()[1]*0.02, gr->GetY()[5]+gr->GetY()[1]*0.15," #font[42]{^{58}Ni}");
  latexNi->SetTextSize(textSize);
  TLatex *latexZn = new TLatex(gr->GetX()[6]-gr->GetX()[1]*0.02, gr->GetY()[6]+gr->GetY()[1]*0.15," #font[42]{^{64}Zn}");
  latexZn->SetTextSize(textSize);
  TLatex *latexGe = new TLatex(gr->GetX()[7]-gr->GetX()[1]*0.02, gr->GetY()[7]+gr->GetY()[1]*0.15," #font[42]{^{70}Ge}");
  latexGe->SetTextSize(textSize);
  TLatex *latexLine = new TLatex(27, gr->GetY()[1]+gr->GetY()[1]*0.1," #font[42]{7.9 #times 10^{-16}}");
  latexLine->SetTextSize(textSize);
  latexLine->SetTextColor(kRed);
  
  gr->GetListOfFunctions()->Add(line);
  gr->GetListOfFunctions()->Add(latexS);
  gr->GetListOfFunctions()->Add(latexCa);
  gr->GetListOfFunctions()->Add(latexTi);
  gr->GetListOfFunctions()->Add(latexCr);
  gr->GetListOfFunctions()->Add(latexFe);
  gr->GetListOfFunctions()->Add(latexNi);
  gr->GetListOfFunctions()->Add(latexZn);
  gr->GetListOfFunctions()->Add(latexGe);
  gr->GetListOfFunctions()->Add(latexLine);

  gr->GetXaxis()->SetLimits(xMin,xMax);
  gr->GetXaxis()->SetTitle("Atomic Number");
  gr->GetXaxis()->SetTitleOffset(1.13);
  gr->GetYaxis()->SetTitle("Sensitivity (90% C.L.)");
  gr->SetMaximum(7.5*Power(10,-15));
  gr->SetTitle("");

  gr->Draw("AP");
}
