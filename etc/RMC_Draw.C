#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include <iostream>
#include <sstream>

void RMC_Draw(){

  TCanvas *c = new TCanvas("c","c",400,400);
  
  TF1 *rmc = new TF1("rmc","(1-2*x+2*x*x)*x*(1-x)*(1-x)",0,1);
  rmc->Draw();

 
}
