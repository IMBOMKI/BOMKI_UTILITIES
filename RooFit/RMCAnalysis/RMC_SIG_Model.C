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
#include "RooTFnBinding.h"
#include "RooPlot.h"
#include "RooClassFactory.h"
#include "RooELossPdf.h"
#include "Math/DistFunc.h"
#include "RooRandom.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <TMath.h>
#include <map>
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"
#include <iomanip>
#include <fstream>

using namespace RooFit;
using namespace RooStats;
using namespace TMath;

////////////////////////////////////////////////
//
//  root> gSystem->Load("RooELossPdf_cxx.so");
//  root> .x RMC_SIG_Model.C("ele")
//
////////////////////////////////////////////////

Double_t fitf(Double_t *x, Double_t *par)
{
  Double_t fitval = par[0]*(pow((par[2]-x[0]),par[1]));
  return fitval;
}

std::pair<double,double> FindClosestKey(std::map<double,double> map, double b){
  double key=-1;
  double val;
  double diff=abs(key-b);
  for (std::map<double,double>::iterator it=map.begin(); it!=map.end(); ++it){
    if (abs(it->first-b)<diff){
      key  = it->first;
      val  = it->second; 
      diff = abs(it->first-b);
    }
  }
  return std::make_pair(key,val);
}

void RMC_SIG_Model(string element){
  gSystem->Load("RooELossPdf_cxx.so");
  gSystem->Load("RooPolyFitPdf_cxx.so");

  ////////////////////////////////////////////
  // Set SIG MODEL and Find Momentum Window //
  ////////////////////////////////////////////

  Bool_t ifUseInternalRMC=1;

  int A; //Atomic Mass
  Double_t par0;
  Double_t par1;
  Double_t par2;
  Double_t par3;

  Double_t lmean;
  Double_t lvar;
  Double_t sigE;
  Double_t rmc_mean;
  Double_t rmc_var;
  Double_t BR_rmc;
  Double_t rmc_end;
  Double_t NumOfStoppedMu=pow(10,18);
  Double_t AcceptanceAl=0.01;
  Double_t TimeFOM_Al=0.191;
  Double_t TimeFOM;
  Double_t Acceptance;
  Double_t lowB;
  Double_t upB=120; 
  //if (element=="Al") upB=94;
  if (element=="Al") upB=120;
  Double_t ProbGND=0.9;
  Double_t fcap;
  Int_t BinNumber;

  Int_t N=48;
  if (element=="Al") N=2;
  Double_t RMCAcceptance[N];
  Double_t SigAcceptance[N];
  Double_t Sensitivity[N];
  Double_t Nup[N];
  Double_t b[N];
  Double_t OptWindowMin[N];
  Double_t OptWindowMax[N];
  
  std::map<double, double> BkgToSig; 

  // Sensitivity Map
  std::string fileDir = "/home/bomki/workspace/BOMKI_UTILITIES/Mathematica/";
  std::string fileNames[4] = {"Integer50.dat","Integer300.dat","Fine.dat","SuperFine.dat"};
  for (int i=0; i<sizeof(fileNames)/sizeof(fileNames[0]); i++){
    std::ifstream file;
    std::string fullFileName=fileDir+fileNames[i];
    file.open(fullFileName);
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream ss(line);
      double Nb, Nup;
      ss >> Nb >> Nup; 
      BkgToSig[Nb]=Nup;
      //std::cout << Nb << " " << Nup << std::endl;
    }
  }

  // Values in the paper of feldman & cousin 
  /*
  BkgToSig[0.]= 2.44;  
  BkgToSig[0.5]= 2.86;
  BkgToSig[1.]= 3.28;
  BkgToSig[1.5]= 3.62;
  BkgToSig[2.]= 3.94;
  BkgToSig[2.5]= 4.20;
  BkgToSig[3.]= 4.42;
  BkgToSig[4.]= 4.83;
  BkgToSig[5.]= 5.18;
  BkgToSig[6.]= 5.53;
  BkgToSig[7.]= 5.90;
  BkgToSig[8.]= 6.18;
  BkgToSig[9.]= 6.49;
  BkgToSig[10.]= 6.76;
  BkgToSig[11.]= 7.02;
  BkgToSig[12.]= 7.28;
  BkgToSig[13.]= 7.51;
  BkgToSig[14.]= 7.75;
  BkgToSig[15.]= 7.99;
  */    

  // Set pre-determined values
  if (element=="Al"){
    lmean=0.433262;
    lvar=0.10388;
    sigE = 92.3;
    lowB=90.30;
    //BinNumber = 30;
    //BinNumber = (upB-90)/0.5
    BR_rmc = 6.21751*pow(10,-7);      
    par0=1.1024*pow(10,-1);
    par1=3.41745;
    rmc_end=101.34;
    A=27;
    fcap=0.61;
  }
  else if (element=="S"){
    lmean = 0.338;
    lvar  = 0.0683;
    sigE  = 101.8;
    //rmc_mean = -1.2533;
    //rmc_var  = 1.83163;
    rmc_mean  = -1.88921;
    rmc_var  = 2.01315;
    lowB  = 96.0;
    //BinNumber = (upB-90)/0.5
    BinNumber = 58;
    BR_rmc = 1.38693*pow(10,-7);      
    TimeFOM=0.142;
    fcap = 0.75;
    par0=1.66554*pow(10,-1);
    par1=3.68714;
    rmc_end=102.03;
    A=32;
  }
  else if (element=="Ca"){
    lmean = 0.375;
    lvar  = 0.082;
    sigE  = 103.55;
    //rmc_mean  = -1.10043;
    //rmc_var  = 1.79044;
    rmc_mean  = -1.48143;
    rmc_var  = 1.92208;   
    lowB  = 96.0;
    BinNumber = 58;
    BR_rmc = 1.40288*pow(10,-7);
    TimeFOM=0.078;
    fcap = 0.8507;
    par0=2.34185*pow(10,-1);
    par1=3.50945;
    rmc_end=102.06;
    A=40;
  }
  else if (element=="Ti"){ //Bad Fit
    lmean = 0.484;
    lvar  = 0.136;
    sigE  = 98.89;
    //rmc_mean = -1.33979;
    //rmc_var  = 2.1369;
    rmc_mean  = -1.95096;
    rmc_var  = 2.3384;
    lowB  = 92.0;
    BinNumber = 66;
    BR_rmc = 1.59901*pow(10,-7);
    TimeFOM=0.076;
    fcap = 0.8529;
    par0=3.35934*pow(10,-1);
    par1=3.51422;
    rmc_end=99.17;
    A=48;
  }
  else if (element=="Cr"){
    lmean = 0.572;
    lvar  = 0.225;
    sigE  = 104.06;
    rmc_mean  = -1.30725;
    rmc_var  = 1.80533;
    lowB  = 96.0;
    BinNumber = 58;
    BR_rmc = 8.076*pow(10,-8);
    TimeFOM=0.0378;
    fcap = 0.8939;
    par0=4.06216*pow(10,-1);
    par1=3.56581;
    rmc_end=101.86;
    A=50;
  }
  else if (element=="Fe"){
    lmean = 0.610;
    lvar  = 0.223;
    sigE  = 103.3;
    rmc_mean  = -1.41425;
    rmc_var  = 1.84893;
    lowB  = 96.0;
    BinNumber = 58;
    BR_rmc = 1.09183*pow(10,-7);
    TimeFOM=0.0269;
    fcap = 0.9085;
    par0=4.09697*pow(10,-1);
    par1=3.64415;
    rmc_end=101.93;
    A=54;
  }
  else if (element=="Ni"){ 
    lmean = 0.672;
    lvar  = 0.212;
    sigE  = 104.25;
    rmc_mean  = -1.48623;
    rmc_var  = 1.87357;
    lowB  = 96.0;
    BinNumber = 58;
    BR_rmc = 0.952345*pow(10,-7);
    TimeFOM=0.0093;
    fcap=0.9306;
    par0=4.29729*pow(10,-1);
    par1=3.62841;
    rmc_end=101.95;
    A=58;
  }
  else if (element=="Zn"){
    lmean = 0.567;
    lvar  = 0.206;
    sigE  = 103.1;
    rmc_mean  = -1.24976;
    rmc_var  = 1.67151;
    lowB  = 96.0;
    BinNumber = 58;
    BR_rmc = 7.70664*pow(10,-8);
    TimeFOM=0.0112;
    fcap=0.9299;
    par0=3.80561*pow(10,-1);
    par1=3.685;
    rmc_end=101.43;
    A=64;
  }
  else if (element=="Ge"){
    lmean = 0.499;
    lvar  = 0.155;
    sigE  = 100.67;
    rmc_mean  = -1.05651;
    rmc_var  = 1.26261;
    lowB  = 96.0;
    BinNumber = 58;
    BR_rmc = 3.81926*pow(10,-8);
    TimeFOM=0.0133;
    fcap = 0.9272;    
    par0=6.19571*pow(10,-1);
    par1=3.80235;
    par2=2.36536;
    rmc_end=100.02;
    A=70;
  }
  else return;

  if (element!="Al") Acceptance = AcceptanceAl*TimeFOM/TimeFOM_Al;
  else if (element=="Al") Acceptance=AcceptanceAl;

  Double_t lowBshift=20;
  if (element=="Al") lowBshift=0.3;

  BinNumber = (upB-lowB+lowBshift)/0.5;

  RooRealVar lowerBound("lowerBound","lowerBound", lowB-lowBshift);
  RooRealVar LOWB("LOWB","LOWB", lowB);
  RooRealVar upperBound("lowerBound","lowerBound", upB);
  RooRealVar x("x","x",lowerBound.getValV(),upperBound.getValV());
  x.setBins(6000);
  if (element=="Al") x.setBins(100);

  // Theoretical Model (Intrinsic)
  RooRealVar landauMean("landauMean","mean of Landau",lmean);
  RooRealVar landauVar("landauVar","variance of Landau", lvar);
  RooRealVar Par0("Par0","Par0",par0);
  RooRealVar Par1("Par1","Par1",par1);
  RooRealVar Par2("Par2","Par2",rmc_end);

  RooRealVar sigEVal("sigEVal","signal energy", sigE);
  RooGenericPdf sigPdf_g4("sigG4","sigG4","Landau(-x+sigEVal,landauMean,landauVar)",RooArgSet(x,sigEVal,landauMean,landauVar));  

  // Smearing from Detector Response 
  RooRealVar gausMean("gausMean","gausMean",0);
  RooRealVar gausVar("gausVar","gausVar",0.2);
  RooGaussian smearPdf("smear","smear",x,gausMean,gausVar);
 
  // Convolution
  RooFFTConvPdf sigPdf_detResp("sig_DetResp","sig_DetResp",x,sigPdf_g4,smearPdf);

  // Draw the models
  TCanvas* c1 = new TCanvas("Sig_Model","Sig_Model",600,600);  
  RooPlot* modelFrame = x.frame(Title(""));
  modelFrame->SetTitle("");
  sigPdf_g4.plotOn(modelFrame);
  sigPdf_detResp.plotOn(modelFrame,LineColor(kRed));
  //sigCdf_detResp->plotOn(modelFrame,LineColor(kGreen));
  modelFrame->Draw();    

  ///////////////////
  // RMC Parameter //
  ///////////////////

  TCanvas* c2 = new TCanvas("RMC_Model","RMC_Model",500,1000);  
  c2->Divide(1,2);
  c2->cd(1);
  
  TString filename = "Simplified_"+TString(element)+".root";
  TFile *f =  TFile::Open(filename);
  TTree *t = (TTree*)f->Get("trdata");

  Double_t Pairep_genTrE;
  Bool_t ifPairProdOccurs;
  Bool_t ifPairVertexAtTarget;
  t->SetBranchAddress("Pairep_genTrE",&Pairep_genTrE);
  t->SetBranchAddress("ifPairProdOccurs",&ifPairProdOccurs);
  t->SetBranchAddress("ifPairVertexAtTarget",&ifPairVertexAtTarget);

  int N_pair = t->Draw("Pairep_genTrE","ifPairProdOccurs==1"); 
  int N_VatT = t->Draw("Pairep_genTrE","ifPairProdOccurs==1 && ifPairVertexAtTarget==1"); 
  int N_mom  = t->Draw("Pairep_genTrE",Form("ifPairProdOccurs==1 && ifPairVertexAtTarget==1 && Pairep_genTrE>%f",lowB)); 


  ///////////////
  // RMC Model //
  ///////////////

  c2->cd(2);

  RooRealVar gausMean_rmc("gausMean_rmc", "gausMean_rmc",rmc_mean);
  RooRealVar gausVar_rmc("gausVar_rmc", "gausVar_rmc",rmc_var);
  RooFormulaVar gausMean_shifted("gausMean_shifted", "gausMean_rmc+LOWB",RooArgSet(gausMean_rmc,LOWB));
  RooGaussian rmcPdf_shifted("rmc","rmc",x,gausMean_shifted,gausVar_rmc);

  // Loss Function
  RooELossPdf eLoss("eLoss","eLoss",x, landauMean,landauVar);    

  RooFFTConvPdf rmcPdf_landau("rmcPdf_landau","rmcPdf_landau",x,rmcPdf_shifted,eLoss);    
  RooFFTConvPdf rmcPdf_detResp("rmcPdf_detResp","rmcPdf_detResp",x,rmcPdf_landau,smearPdf);    
  RooPlot* rmcFrame = x.frame(Title(""));
  rmcFrame->SetTitle("RMC_shifted");  

  rmcPdf_shifted.plotOn(rmcFrame,LineColor(kBlue));
  rmcPdf_landau.plotOn(rmcFrame,LineColor(kRed));
  rmcPdf_detResp.plotOn(rmcFrame,LineColor(kGreen));
  rmcFrame->Draw();

  TCanvas * c4=new TCanvas("poly","poly",600,600);

  c4->cd();

  RooPlot* rmcFrame2 = x.frame(Title("")); 

  Par2.setConstant(kTRUE);
  RooPolyFitPdf polyFit("polyFit","polyFit",x,Par0,Par1,Par2);
  polyFit.plotOn(rmcFrame2);

  RooNumConvPdf LandauGaus("landauGaus","landauGaus",x,eLoss,smearPdf);  
  RooFFTConvPdf poly_detResp("rmcPdf_detResp","rmcPdf_detResp",x,polyFit,LandauGaus);
  poly_detResp.plotOn(rmcFrame2,LineColor(kBlack));
  //RooPolyFitPdf poly_detResp("polyFit","polyFit",x,Par0,Par1,Par2);  


  rmcFrame2->Draw();

  // Find Energy window with 90 % acceptance
  Double_t Window_ELow;
  Double_t stepSize=0.01;
  Int_t itN=400;
  Double_t tolerance=0.01;

  for (Int_t i_acc=0; i_acc<N; i_acc++){
    for (Int_t i=0; i<itN; i++){
      x.setRange("window",sigE-stepSize*i,upB);
      RooAbsReal* isig = sigPdf_detResp.createIntegral(x,NormSet(x),Range("window"));
      if (TMath::Abs(isig->getVal()-(i_acc+1)*(0.96/N)) < tolerance) {
	Window_ELow=sigE-stepSize*i;
	//std::cout << isig->getVal() << std::endl;
	std::cout << "Signal Acceptance in Momentum Window: " << TMath::Abs(isig->getVal()) << std::endl;
	std::cout << "Lower Bound of E-Window: " << Window_ELow << std::endl;
	break;
      }
    }
    x.setRange("window",Window_ELow,Window_ELow+stepSize*itN);  
    x.setRange("Norm",  lowB,       Window_ELow+stepSize*itN);  
    OptWindowMin[i_acc]=Window_ELow;
    OptWindowMax[i_acc]=Window_ELow+stepSize*itN;

    RooAbsReal* isig = sigPdf_detResp.createIntegral(x,NormSet(x),Range("window"));
    RooAbsReal* irmc_detResp = poly_detResp.createIntegral(x,NormSet(x),Range("window"));
    RooAbsReal* irmc_detResp_Norm = poly_detResp.createIntegral(x,NormSet(x),Range("Norm"));

    SigAcceptance[i_acc] = TMath::Abs(isig->getVal());
    RMCAcceptance[i_acc] =  TMath::Abs(irmc_detResp->getVal()/irmc_detResp_Norm->getVal());
    std::cout << "Lower Bound of E-Window: " << Window_ELow << std::endl;
    std::cout << "Signal Acceptance in Momentum Window: " << SigAcceptance[i_acc] << std::endl;
    std::cout << "RMC Acceptance in Momentum Window: " << RMCAcceptance[i_acc] << std::endl;    
    b[i_acc] = NumOfStoppedMu*BR_rmc*N_mom/t->GetEntries()*Acceptance*RMCAcceptance[i_acc];
    if (ifUseInternalRMC==1) b[i_acc]*=2;
    std::pair<double,double> pair = FindClosestKey(BkgToSig,b[i_acc]);
    Nup[i_acc] = pair.second;
    Double_t Acceptance = 0.01 * TimeFOM/TimeFOM_Al * SigAcceptance[i_acc]; // * (i_acc+1)*(0.95/N);
    Sensitivity[i_acc] = Nup[i_acc]/(NumOfStoppedMu*fcap*ProbGND*Acceptance);
  }

  int opt_Index;
  double opt_Sens=1;
  double opt_Acceptance;
  double opt_Acceptance_RMC;
  double opt_winMin, opt_winMax;

  for (int i=0; i<N; i++){
    std::cout << "Signal Acceptance: " << left << std::setw(12) << SigAcceptance[i] << " "
	      << "b: " << left << std::setw(12) << b[i] << " "
	      << "N up: "  << left << std::setw(8)  << Nup[i] << " "
	      << "Sensitivity: " << left << std::setw(12) << Sensitivity[i] << std::endl;
    if (Sensitivity[i]<opt_Sens){
      opt_Index=i;
      opt_Sens=Sensitivity[opt_Index];
      opt_Acceptance=SigAcceptance[opt_Index];
      opt_Acceptance_RMC=RMCAcceptance[opt_Index];
      opt_winMin=OptWindowMin[opt_Index];
      opt_winMax=OptWindowMax[opt_Index];
    }
  }
  std::cout << "Optimized Sensitivity: " << opt_Sens << std::endl;
  std::cout << "Optimized SignalWidow: " << opt_winMin << "~" << opt_winMax << std::endl;




  ////////////// Test Convolved detection response ///////////////
  TCanvas* c5 = new TCanvas("test","test",600,600);

  RooRealVar x_test("x_test","x_test",-10,10);
  RooELossPdf eLoss_test ("eLoss_test","eLoss_test",x_test, landauMean,landauVar);    
  RooGaussian smearPdf_test("smear_test","smear_test",x_test,gausMean,gausVar);
  RooNumConvPdf LandauGaus_test("landauGaus_test","landauGaus_test",x_test,eLoss_test,smearPdf_test);
  
  RooPlot* testFrame = x_test.frame(Title("")); 
  LandauGaus_test.plotOn(testFrame);
  eLoss_test.plotOn(testFrame,LineColor(kRed));
  smearPdf_test.plotOn(testFrame,LineColor(kGreen));
  testFrame->Draw();
  ////////////////////////////////////////////////////////////////


  /////////////////////////////////////
  // Set Composite Model For Calcium //
  /////////////////////////////////////

  Int_t rmcNum_part=NumOfStoppedMu*fcap*BR_rmc*N_mom/t->GetEntries()*Acceptance; // RMC number from lowB to upB (range of x)

  if (element=="Al") opt_winMax=upB;

  x.setRange("part",lowB,opt_winMax);
  x.setRange("total",lowerBound.getValV(),opt_winMax) ;
  x.setRange("window",opt_winMin,opt_winMax);  


  if (element=="Al"){
    x.setRange("part",lowB,upB-15);
    x.setRange("total",lowerBound.getValV(),upB-15); 
  }

  RooAbsReal* irmc_part   = poly_detResp.createIntegral(x,NormSet(x),Range("part"));
  RooAbsReal* irmc_total  = poly_detResp.createIntegral(x,NormSet(x),Range("total"));
  RooAbsReal* irmc_window = poly_detResp.createIntegral(x,NormSet(x),Range("window"));

  Int_t rmcNum=rmcNum_part*irmc_total->getVal()/irmc_part->getVal(); // RMC number from lowerBound.getValV() to upB (range of x)


  if (element=="Al"){
    x.setRange("redundant",105,upB);  
    RooAbsReal* irmc_redundant = poly_detResp.createIntegral(x,NormSet(x),Range("redundant"));
    Int_t rmcNum_redundant = rmcNum*irmc_redundant->getVal()/irmc_total->getVal();
    rmcNum = rmcNum+rmcNum_redundant;
  }

  if (ifUseInternalRMC==1)rmcNum*=2; // Count Internal RMC
  //Double_t muepBr=opt_Sens; 
  RooRandom::randomGenerator()->SetSeed(1);  
  Double_t muepBr=TMath::Power(10,-14); 
  if (element=="Al") muepBr=1.7*TMath::Power(10,-13);
  Int_t sigNum=NumOfStoppedMu*fcap*muepBr*Acceptance; // SIG number

  RooRealVar sigFrac("sigFrac", "Fraction of RMC", Double_t(sigNum)/(sigNum+rmcNum),0.,1.);
  RooRealVar rmcFrac("rmcFrac", "Fraction of RMC", Double_t(rmcNum)/(sigNum+rmcNum));
  //RooBinning bins(BinNumber,lowB,upB);
  RooBinning bins(BinNumber,lowerBound.getValV(),upB);

  /* Make Composite Comp */
  RooAddPdf compPdf("comp", "comp", RooArgList(sigPdf_detResp, poly_detResp),sigFrac);
  RooDataSet *compData=compPdf.generate(x, (sigNum+rmcNum));

  /* Plot Composite Comp */
  RooPlot* compFrame = x.frame(Title(""));
  compFrame->SetTitle("");
  compData->plotOn(compFrame,Binning(bins));
  compPdf.plotOn(compFrame, Components(poly_detResp), LineStyle(9), LineColor(kBlue));
  compPdf.plotOn(compFrame, LineStyle(kDashed), LineColor(kRed));

  TCanvas* c3 = new TCanvas("Composite_Comp","Composite_Comp",600,600);  
  c3->cd() ; 
  gPad->SetLeftMargin(0.15) ; 
  compFrame->GetYaxis()->SetTitleOffset(1.23) ; 
  compFrame->GetXaxis()->SetTitle("E_{e^{+}} (MeV)");
  compFrame->GetXaxis()->SetTitleOffset(1.23);
  if (element=="Al") compFrame->GetYaxis()->SetTitleOffset(1.5);

  Double_t xAxisMin=101.0;
  Double_t xAxisMax=104.5;
  Double_t yAxisMin=0;
  Double_t yAxisMax=10.0;

  if (element=="Ca") {
    xAxisMin=100.5;
    xAxisMax=104.5;
    yAxisMax=40.0;
  }

  else if (element=="S") {
    xAxisMin=100.5;
    xAxisMax=103.0;
    yAxisMax=70.0;
  }

  else if (element=="Ti") {
    xAxisMin=97.5;
    xAxisMax=100.0;
    yAxisMax=40.0;
  }

  else if (element=="Fe" || element=="Cr") {
    xAxisMin=100.5;
    xAxisMax=104.5;
    yAxisMax=20.0;
  }

  else if (element=="Ni") {
    xAxisMin=100.5;
    xAxisMax=104.5;
    yAxisMax=10.0;
  }

  else if (element=="Zn") {
    xAxisMin=100.0;
    xAxisMax=103.5;
  }
  else if (element=="Ge") {
    xAxisMin=98.5;
    xAxisMax=101.0;
    yAxisMax=20.0;
  }
  else if (element=="Al") {
    xAxisMin=90.5;
    xAxisMax=96.5;
    //xAxisMax=upB;
    //yAxisMin=3*Power(10,4);
    yAxisMin=0*Power(10,5);
    yAxisMax=2.0*Power(10,5);
  }
  compFrame->GetXaxis()->SetLimits(xAxisMin,xAxisMax);
  compFrame->SetMaximum(yAxisMax);
  compFrame->SetMinimum(yAxisMin);

  TLine *Line = new TLine(opt_winMin, 0,opt_winMin,yAxisMax);
  Line->SetLineStyle(kDashed);
  TLatex *latexLine  = new TLatex(opt_winMin, yAxisMax*0.65,Form("  E_{e^{+}}>%.2f MeV",opt_winMin));
  latexLine->SetTextSize(0.044);
  TLatex *latexLine2  = new TLatex(xAxisMin+(xAxisMax-xAxisMin)*0.85, yAxisMin+(yAxisMax-yAxisMin)*0.91,Form("^{%d}"+TString(element),A));
  latexLine2->SetTextSize(0.047);



  if(element!="Al"){
    compFrame->addObject(Line);    
    compFrame->addObject(latexLine);    
  }
  compFrame->addObject(latexLine2);    
  compFrame->Draw();

  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << "RMC Num: (90.3~105) " << rmcNum_part << "  RMC Num: (90~105) "  << rmcNum << "  " << "Sig Num: " << sigNum << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << Double_t(N_pair)/t->GetEntries() << "  " << Double_t(N_VatT)/Double_t(N_pair) << "  " << Double_t(N_mom)/Double_t(N_VatT) << std::endl;
  std::cout << "///////////////////////////////////////" << std::endl;
  std::cout << Double_t(N_pair)/t->GetEntries() << "  " << Double_t(N_VatT)/Double_t(N_pair) << "  " << Double_t(N_mom)/Double_t(N_VatT) << std::endl; 
  std::cout << "///////////////////////////////////////" << std::endl;

  //////////////////////////
  // Draw Boxed Histogram //
  //////////////////////////

  if (element=="S") {
    xAxisMin=100.5;
    xAxisMax=103.0;
    yAxisMax=100.0;
  }

  else if (element=="Ti") {
    xAxisMin=97.5;
    xAxisMax=99.5;
    yAxisMax=50.0;
  }


  Double_t HistBinSize=0.1;
  Int_t HistBinNumber= (xAxisMax-xAxisMin)/HistBinSize;
  TCanvas* c6 = new TCanvas("hist","hist",600,600);
  TH1D *h_rmc = new TH1D("hist_rmc","hist_rmc", HistBinNumber ,xAxisMin,xAxisMax);
  TH1D *h_sig = new TH1D("hist_sig","hist_sig", HistBinNumber ,xAxisMin,xAxisMax);
  if (element=="Al") TGaxis::SetMaxDigits(3);

  Double_t bin_start;
  Double_t bin_end;
  Double_t histY_rmc;
  Double_t histY_sig;
  RooAbsReal* isig_part = sigPdf_detResp.createIntegral(x,NormSet(x),Range("total"));

  for (Int_t i=0; i< HistBinNumber; i++){
    bin_start=xAxisMin+i*HistBinSize;
    bin_end=xAxisMin+(i+1)*HistBinSize;
    x.setRange("singleBin",bin_start,bin_end);    

    RooAbsReal* irmc_singleBin = poly_detResp.createIntegral(x,NormSet(x),Range("singleBin"));
    histY_rmc = rmcNum*irmc_singleBin->getVal()/irmc_total->getVal();
    RooAbsReal* isig_singleBin = sigPdf_detResp.createIntegral(x,NormSet(x),Range("singleBin"));
    histY_sig = sigNum*isig_singleBin->getVal()/isig_part->getVal();
    h_rmc->SetBinContent(i+1,histY_rmc);
    h_sig->SetBinContent(i+1,histY_sig);   

    std::cout << i << "  " << bin_start << "  " << bin_end << "  " << histY_rmc << "  " << histY_sig << std::endl;
  }
  h_rmc->SetLineColor(kBlue);
  h_sig->SetLineColor(kRed);
  THStack *st = new THStack("st","Stacked 2D histograms");
  st->Add(h_rmc);
  st->Add(h_sig);
  st->SetTitle("");
  st->Draw();
  c6->Update();

  Double_t yHistMax = gPad->GetUymax();
  TLine *Line_hist = new TLine(opt_winMin, 0,opt_winMin,yHistMax);
  Line_hist->SetLineStyle(kDashed);
  TLatex *latexLine_hist  = new TLatex(opt_winMin, yHistMax*0.7,Form("  E_{e^{+}}>%.2f MeV",opt_winMin));
  latexLine_hist->SetTextSize(0.044);
  yAxisMin=0;
  TLatex *latexLine2_hist  = new TLatex(xAxisMin+(xAxisMax-xAxisMin)*0.85, (yHistMax)*0.91,Form("^{%d}"+TString(element),A));
  latexLine2_hist->SetTextSize(0.047);

  if(element!="Al"){
    Line_hist->Draw();    
    latexLine_hist->Draw();    
  }
  latexLine2_hist->Draw();    

  st->GetXaxis()->SetTitleOffset(1.23);
  st->GetYaxis()->SetTitleOffset(1.23) ; 
  st->GetXaxis()->SetTitle("E_{e^{+}} (MeV)");
  st->GetYaxis()->SetTitle("Counts per 0.1 MeV");

 
  ////////////////////////////////
  // RooStats Likelihood Method //
  ////////////////////////////////

  /* Set up Model with ProfileLikelihoodCalculator */

  if (element=="Al"){
      
    //Fitting PolyDetResp 
    gStyle->SetOptFit(111);

    TF1 *func = new TF1("fitf",fitf,lowerBound.getValV(),upB-15,3); // Same with the range of "total"
    func->SetParNames("par0","par1","par2");
    func->SetParLimits(0,0.,1000.);
    func->SetParLimits(1,0.,1000.);
    func->SetParLimits(2,95.,upB-15);
    h_rmc->Fit("fitf","r");

    Double_t p0=func->GetParameter(0);
    Double_t p1=func->GetParameter(1);
    Double_t p2=func->GetParameter(2);

    RooRealVar P0("Par0_Al","Par0_Al",p0);
    RooRealVar P1("Par1_Al","Par1_Al",p1);
    RooRealVar P2("Par2_Al","Par2_Al",p2);

    Double_t x_Al_start = 90.5;
    Double_t x_Al_end   = 94.0;

    RooRealVar x_Al("x_Al","x_Al",x_Al_start,x_Al_end);    
    RooPolyFitPdf polyFit_Al("polyFit_Al","polyFit_Al",x_Al,P0,P1,P2);

    x.setRange("total_Al",x_Al_start,x_Al_end);  
    RooAbsReal* irmc_total_Al = poly_detResp.createIntegral(x,NormSet(x),Range("total_Al"));
    Int_t rmcNum_total_Al = rmcNum_part*irmc_total_Al->getVal()/irmc_part->getVal();
    //rmcNum = rmcNum-rmcNum_redundant;
    rmcNum=rmcNum_part*irmc_total_Al->getVal()/irmc_part->getVal();
    RooRealVar sigFrac_Al("sigFrac_Al", "Fraction of RMC", Double_t(sigNum)/(sigNum+rmcNum),0.,0.01);
    
    RooGenericPdf sigPdf_g4_Al("sigG4","sigG4","Landau(-x_Al+sigEVal,landauMean,landauVar)",RooArgSet(x_Al,sigEVal,landauMean,landauVar));  
    RooGaussian smearPdf_Al("smear","smear",x_Al,gausMean,gausVar);
    RooFFTConvPdf sigPdf_detResp_Al("sig_DetResp","sig_DetResp_Al",x_Al,sigPdf_g4_Al,smearPdf_Al);
    //RooAddPdf compPdf2("comp2", "comp2", RooArgList(sigPdf_detResp_Al, polyFit_Al),sigFrac_Al);    
    RooAddPdf compPdf2("comp2", "comp2", RooArgList(sigPdf_g4_Al, polyFit_Al),sigFrac_Al);    
    
    RooDataSet *compData2=compPdf2.generate(x_Al, (sigNum+rmcNum));
    Double_t BinNumber_Al=7;
    RooBinning bins_Al(BinNumber_Al,x_Al_start,x_Al_end);

    TCanvas* c7 = new TCanvas("comp2","comp2", 600,600);
    RooPlot* compFrame2 = x_Al.frame(Title(""));
    compFrame2->SetTitle("");
    compData2->plotOn(compFrame2,Binning(bins_Al));
    compPdf2.plotOn(compFrame2, Components(polyFit_Al), LineStyle(9), LineColor(kBlue));
    compPdf2.plotOn(compFrame2, LineStyle(kDashed), LineColor(kRed));
    compFrame2->Draw();
    
    ModelConfig model;
    RooWorkspace* wks = new RooWorkspace("wks");
    wks->import(compPdf2);    
    wks->import(*compData2, Rename("data2"));
    model.SetWorkspace(*wks);
    model.SetPdf("comp2");
    
    ProfileLikelihoodCalculator plc;
    plc.SetData(*(wks->data("data2")));
    RooRealVar* fsig = wks->var("sigFrac_Al");
    RooArgSet poi(*fsig);
    RooArgSet *nullParams = (RooArgSet*)poi.snapshot();
    nullParams->setRealValue("sigFrac_Al",0);
    
    plc.SetModel(model);
    plc.SetNullParameters(*nullParams);
    //plc.SetNullParameters(*nullParams);                                            
    plc.SetModel(model);
    // NOTE: using snapshot will import nullparams                                   
    // in the WS and merge with existing "mu"                                        
    // model.SetSnapshot(*nullParams);                                               
    
    //use instead setNuisanceParameters                                              
    plc.SetNullParameters( *nullParams);
    
    // We get a HypoTestResult out of the calculator, and we can query it.           
    HypoTestResult* htr = plc.GetHypoTest();
    cout << "-------------------------------------------------" << endl;
    cout << "The p-value for the null is " << htr->NullPValue() << endl;
    cout << "Corresponding to a significance of " << htr->Significance() << endl;
    cout << "-------------------------------------------------\n\n" << endl;
    
    //MakePlots(wks);
    delete wks;  
    
  }

}
