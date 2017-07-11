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
  //Double_t S  = 1.995*Power(10,-15);
  Double_t Ca = 7.91461*Power(10,-16);
  Double_t Ti = 2.15718*Power(10,-15);
  //Double_t Cr = 1.74832*Power(10,-15);
  Double_t Fe = 2.59838*Power(10,-15);
  Double_t Ni = 6.38783*Power(10,-15);
  //Double_t Zn = 5.88725*Power(10,-15);
  //Double_t Ge = 6.50173*Power(10,-15);

  const Int_t n = 4;
  //TString ele[8] = {"S","Ca","Ti","Cr","Fe","Ni","Ge"};
  //Double_t ele[8] = {1,2,3,4,5,6,7,8};
  //Double_t Z[8]   = {16,20,22,24,26,28,30,32};
  //Double_t ses[8] = {S,Ca,Ti,Cr,Fe,Ni,Zn,Ge};

  Double_t Z[n]   = {20,22,26,28};
  Double_t ses[n] = {Ca,Ti,Fe,Ni};

  TCanvas *c = new TCanvas("canvas","canvas",600,600);
  TGraph * gr = new TGraph(n,Z,ses);
  gr->SetMarkerStyle(33);
  gr->SetMarkerSize(2.7);
  gr->SetMarkerColor(kBlue);

  Double_t xMin=18.;
  Double_t xMax=30.;

  TLine *line = new TLine(xMin,Ca,xMax,Ca);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);
  //TLatex *latexS  = new TLatex(gr->GetX()[0], gr->GetY()[0]+gr->GetY()[1]*0.05," ^{32}S");
  //latexS->SetTextSize(0.03);
  TLatex *latexCa = new TLatex(gr->GetX()[0], gr->GetY()[0]+gr->GetY()[1]*0.05," ^{40}Ca");
  latexCa->SetTextSize(0.047);
  TLatex *latexTi = new TLatex(gr->GetX()[1], gr->GetY()[1]+gr->GetY()[1]*0.05," ^{48}Ti");
  latexTi->SetTextSize(0.047);
  //TLatex *latexCr = new TLatex(gr->GetX()[3], gr->GetY()[3]+gr->GetY()[1]*0.05," ^{50}Cr");
  //latexCr->SetTextSize(0.03);
  TLatex *latexFe = new TLatex(gr->GetX()[2], gr->GetY()[2]+gr->GetY()[1]*0.05," ^{54}Fe");
  latexFe->SetTextSize(0.047);
  TLatex *latexNi = new TLatex(gr->GetX()[3], gr->GetY()[3]+gr->GetY()[1]*0.05," ^{58}Ni");
  latexNi->SetTextSize(0.047);
  //TLatex *latexZn = new TLatex(gr->GetX()[6], gr->GetY()[6]+gr->GetY()[1]*0.05," ^{64}Zn");
  //latexZn->SetTextSize(0.03);
  //TLatex *latexGe = new TLatex(gr->GetX()[7], gr->GetY()[7]+gr->GetY()[1]*0.05," ^{70}Ge");
  //latexGe->SetTextSize(0.03);
  TLatex *latexLine = new TLatex(26, gr->GetY()[0]+gr->GetY()[0]*0.1," 7.9 #times 10^{-16}");
  latexLine->SetTextSize(0.047);
  latexLine->SetTextColor(kRed);
  
  gr->GetListOfFunctions()->Add(line);
  //gr->GetListOfFunctions()->Add(latexS);
  gr->GetListOfFunctions()->Add(latexCa);
  gr->GetListOfFunctions()->Add(latexTi);
  //gr->GetListOfFunctions()->Add(latexCr);
  gr->GetListOfFunctions()->Add(latexFe);
  gr->GetListOfFunctions()->Add(latexNi);
  //gr->GetListOfFunctions()->Add(latexZn);
  //gr->GetListOfFunctions()->Add(latexGe);
  gr->GetListOfFunctions()->Add(latexLine);

  gr->GetXaxis()->SetLimits(xMin,xMax);
  //gr->GetXaxis()->SetLabelSize(0);
  gr->GetXaxis()->SetTitle("Atomic Number");
  gr->GetXaxis()->SetTitleOffset(1.13);
  gr->GetYaxis()->SetTitle("Sensitivity (90% C.L.)");
  gr->SetMaximum(7.5*Power(10,-15));
  gr->SetTitle("");

  gr->Draw("AP");
}
