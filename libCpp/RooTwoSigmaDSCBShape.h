#ifndef ROO_TWOSIGMADSCB_SHAPE
#define ROO_TWOSIGMADSCB_SHAPE


#include "RooAbsPdf.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h"
#include "Riostream.h"

class RooTwoSigmaDSCBShape : public RooAbsPdf {
public:
  RooTwoSigmaDSCBShape() {}
  RooTwoSigmaDSCBShape(const char* name, const char* title,
                  RooAbsReal& _x,
                  RooAbsReal& _mean,
                  RooAbsReal& _sigmaL,
                  RooAbsReal& _sigmaR,
                  RooAbsReal& _alphaL,
                  RooAbsReal& _nL,
                  RooAbsReal& _alphaR,
                  RooAbsReal& _nR);

  RooTwoSigmaDSCBShape(const RooTwoSigmaDSCBShape& other, const char* name);
  inline virtual TObject* clone(const char* newname) const { return new RooTwoSigmaDSCBShape(*this,newname);}
  inline ~RooTwoSigmaDSCBShape(){}
  Double_t evaluate() const ;
  
  ClassDef(RooTwoSigmaDSCBShape, 2)

protected:
  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy nL;
  RooRealProxy alphaR;
  RooRealProxy nR;

};

#endif
