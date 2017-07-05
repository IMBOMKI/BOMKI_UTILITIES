#include "TCanvas.h"
#include "TGraph.h"
#include "TPad.h"
#include "TAxis.h"
#include "TH1.h"
#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

using namespace TMath;

void FeldmanCousinTest(){

  TF1* pois=new TF1("poiss2","Pow((x+[0]),[1])*Exp(-(x+[0]))/Factorial([1])",0,100);
  
  // Create the function and wrap it
  TF1 fsin("Sin Function", "sin(x)", TMath::PiOver2(), TMath::TwoPi() );
  ROOT::Math::WrappedTF1 wf1(fsin);
 
  // Create the Integrator
  ROOT::Math::BrentRootFinder brf;
 
  // Set parameters of the method
  brf.SetFunction( wf1, TMath::PiOver2(), TMath::TwoPi() );
  brf.Solve();
 
  cout << brf.Root() << endl;

}
