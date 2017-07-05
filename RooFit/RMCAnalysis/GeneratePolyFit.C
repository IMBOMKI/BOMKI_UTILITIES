#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooLandau.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooBinning.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooChi2Var.h"
#include "RooAbsArg.h"
#include "RooFitResult.h"
#include "RooFit.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TMath.h"

//using namespace TMath;
using namespace RooFit;

void GeneratePolyFit(){ 
  RooRealVar x("x","x",0.,120.);
  RooRealVar par0("par0","par0",0.,0.,100000.);
  RooRealVar par1("par1","par1",0.,0.,10000.);
  RooRealVar par2("par2","par2",0.,0.,100.);
  //RooRealVar par3("par3","par3",0.,-100.,120.);
  //RooRealVar par4("par4","par4",0.,-100.,120.);

  //RooAbsPdf* polyFit =  RooClassFactory::makePdfInstance("polyFit","(par0+par1*Power(x,2)+par2*Power(x,4))*(par3-x)",RooArgSet(x,par0,par1,par2,par3));
  RooAbsPdf* polyFit = RooClassFactory::makePdfInstance("polyFit",
     "(par0*(TMath::Power((par2-x),par1)))",
     RooArgSet(x,par0,par1,par2)
     );

  gROOT->ProcessLine(".L ./RooPolyFit.cxx+");
  return;
}
