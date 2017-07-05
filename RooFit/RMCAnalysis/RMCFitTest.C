#include "TROOT.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"

using namespace TMath;

Double_t fitf(Double_t *x, Double_t *par)
{
  //Double_t arg = 0;
  //if (par[2] != 0) arg = (x[0] - par[1])/par[2];
  //Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  //Double_t fitval = (par[0]+par[1]*x[0]*x[0]+par[2]*Power(x[0],4))*(par[3]-x[0]);

  //Double_t fitval = (par[0]+par[1]*x[0]+par[2]*x[0]*x[0])*(par[3]-x[0]);
  Double_t fitval = //exp(par[0]*(x[0]-par[4]))
    (
     //exp(par[0]*(x[0]-par[6]))
     //par[0]+par[1]*pow(par[6]-x[0],1)
     par[0]*(pow((par[2]-x[0]),par[1]))//+par[2]*pow((par[4]-x[0]),par[3])
     //+par[3]*pow((par[6]-x[0]),4)
     //+par[4]*pow((par[6]-x[0]),9)
     //+par[5]*pow((par[6]-x[0]),5)
     //+par[2]*log(pow(par[4]-x[0],par[3])+1)
     );

  return fitval;
}
void RMCFitTest(string element)
{
  gStyle->SetOptFit(111);
  Double_t lowB;
  Double_t upB;

  if (element=="Al"){
    lowB=90.30;
    upB=101.34;
  }
  else if (element=="S"){
    lowB=96;
    upB=102.03;
  }
  else if (element=="Ca"){
    lowB=96;
    upB=102.06;
  }
  else if (element=="Ti"){
    lowB=92;
    upB=99.17;
  }
  else if (element=="Cr"){
    lowB=96;
    upB=101.86;
  }
  else if (element=="Fe"){
    lowB=96;
    upB=101.93;
  }
  else if (element=="Ni"){
    lowB=96;
    upB=101.95;
  }
  else if (element=="Zn"){
    lowB=96;
    upB=101.43;
  }
  else if (element=="Ge"){
    lowB=96;
    upB=100.02;
  }
  else return;


  TString filename = "Simplified_"+TString(element)+"_Bulk.root";
  TFile *f = new TFile(filename);
  TTree *t = (TTree*)f->Get("trdata");
  TH1D  *h = new TH1D("E","E",100,lowB,upB);
  
  Double_t Pairep_genTrE;
  Bool_t ifPairVertexAtTarget;
  t->SetBranchAddress("Pairep_genTrE",&Pairep_genTrE);
  t->SetBranchAddress("ifPairVertexAtTarget", &ifPairVertexAtTarget);
  
  for (Int_t i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if (Pairep_genTrE>lowB && ifPairVertexAtTarget==1){
      //h->Fill(Pairep_genTrE-lowB);
      h->Fill(Pairep_genTrE);
    }
  }

  TCanvas *c1 = new TCanvas("c1","the fit canvas",600,600);

  TF1 *func = new TF1("fitf",fitf,0.,upB,7);
  //TF1 *func = new TF1("fitf",fitf,0.,upB,3);
  func->SetParNames("par0","par1","par2","par3","par4");
  func->SetParLimits(0,0.,1000.);
  func->SetParLimits(1,0.,1000.);
  //func->SetParLimits(2,0.,1000.);
  //func->SetParLimits(3,0.,1000.);
  func->FixParameter(2,upB);
  //func->FixParameter(2,upB);
  
  h->Fit("fitf","r");

}
