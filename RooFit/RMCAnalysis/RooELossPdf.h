/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOELOSSPDF
#define ROOELOSSPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooELossPdf : public RooAbsPdf {
public:
  RooELossPdf() {} ; 
  RooELossPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _landauMean,
	      RooAbsReal& _landauVar);
  RooELossPdf(const RooELossPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooELossPdf(*this,newname); }
  inline virtual ~RooELossPdf() { }

protected:

  RooRealProxy x ;
  RooRealProxy landauMean ;
  RooRealProxy landauVar ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooELossPdf,1) // Your description goes here...
};
 
#endif
