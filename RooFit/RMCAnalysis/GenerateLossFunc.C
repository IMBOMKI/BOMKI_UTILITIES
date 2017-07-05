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

using namespace RooFit;

void GenerateLossFunc(){ 
  RooRealVar x("x","x",0.,1.);
  RooRealVar landauMean("landauMean","mean of Landau",0.,0.,10.);
  RooRealVar landauVar("landauVar","variance of Landau", 0.,0.,10.);
  RooAbsPdf* energyLoss =  RooClassFactory::makePdfInstance("eLoss","TMath::Landau(-x,landauMean,landauVar)",RooArgSet(x,landauMean,landauVar));
  gROOT->ProcessLine(".L ./RooELossPdf.cxx+");
  return;
}
