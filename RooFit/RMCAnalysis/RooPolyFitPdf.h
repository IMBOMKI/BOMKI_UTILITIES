/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOPOLYFITPDF
#define ROOPOLYFITPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooPolyFitPdf : public RooAbsPdf {
public:
  RooPolyFitPdf() {} ; 
  RooPolyFitPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _par0,
	      RooAbsReal& _par1,
	      RooAbsReal& _par2);
  RooPolyFitPdf(const RooPolyFitPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooPolyFitPdf(*this,newname); }
  inline virtual ~RooPolyFitPdf() { }

protected:

  RooRealProxy x ;
  RooRealProxy par0 ;
  RooRealProxy par1 ;
  RooRealProxy par2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooPolyFitPdf,1) // Your description goes here...
};
 
#endif
