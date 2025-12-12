#include "RooTwoSigmaDSCBShape.h"

ClassImp(RooTwoSigmaDSCBShape)

RooTwoSigmaDSCBShape::RooTwoSigmaDSCBShape(const char *name, const char *title,
                                 RooAbsReal& _x,
                                 RooAbsReal& _mean,
                                 RooAbsReal& _sigmaL,
                                 RooAbsReal& _sigmaR,
                                 RooAbsReal& _alphaL,
                                 RooAbsReal& _nL,
                                 RooAbsReal& _alphaR,
                                 RooAbsReal& _nR
                                ):
RooAbsPdf(name, title),
    x("x", "x", this, _x),
    mean("mean", "mean", this, _mean),
    sigmaL("sigmaL", "sigmaL", this, _sigmaL),
    sigmaR("sigmaR", "sigmaR", this, _sigmaR),
    alphaL("alphaL", "alphaL", this, _alphaL),
    nL("nL", "nL", this, _nL),
    alphaR("alphaR", "alphaR", this, _alphaR),
    nR("nR", "nR", this, _nR)

{}

RooTwoSigmaDSCBShape::RooTwoSigmaDSCBShape(const RooTwoSigmaDSCBShape& other, const char* name): 
    RooAbsPdf(other, name),
    x("x", this, other.x),
    mean("mean", this, other.mean),
    sigmaL("sigmaL", this, other.sigmaL),
    sigmaR("sigmaR", this, other.sigmaR),
    alphaL("alphaL", this, other.alphaL),
    nL("nL", this, other.nL),
    alphaR("alphaR", this, other.alphaR),
    nR("nR", this, other.nR)
{}

Double_t RooTwoSigmaDSCBShape::evaluate() const 
{
  Double_t xx = x;
  Double_t m = mean;
  Double_t sL = sigmaL;
  Double_t sR = sigmaR;
  Double_t aL = alphaL;
  Double_t nL_val = nL;
  Double_t aR = alphaR;
  Double_t nR_val = nR;

  if (xx <= m) {
    Double_t t = (xx - m) / sL;
    if (t >= -aL) {
      return exp(-0.5 * t * t);
    } else {
      Double_t A = TMath::Power(nL_val / aL, nL_val) * exp(-0.5 * aL * aL);
      Double_t B = nL_val / aL - aL;
      return A * TMath::Power(B - t, -nL_val);
    }
  } else {
    Double_t t = (xx - m) / sR;
    if (t <= aR) {
      return exp(-0.5 * t * t);
    } else {
      Double_t A = TMath::Power(nR_val / aR, nR_val) * exp(-0.5 * aR * aR);
      Double_t B = nR_val / aR - aR;
      return A * TMath::Power(B + t, -nR_val);
    }
  }
}