#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TH1.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"

/// include pdfs
#include "RooCBExGaussShape.h"
#include "RooTwoSigmaDSCBShape.h"
#include "RooCMSShape.h"

#include <vector>
#include <string>
#ifdef __CINT__
#pragma link C++ class std::vector<std::string>+;
#endif

using namespace RooFit;
using namespace std;

class tnpFitter {
public:
  tnpFitter( TFile *file, std::string histname  );
  tnpFitter( TH1 *hPass, TH1 *hFail, std::string histname  );
  ~tnpFitter(void) {if( _work != 0 ) delete _work; }
  void setZLineShapes(TH1 *hZPass, TH1 *hZFail, bool zeroLowerThreshold=false, double threshold=60.0);
  void setWorkspace(std::vector<std::string>, bool isaddGaus=false, bool useBreitWigner=false);
  void setOutputFile(TFile *fOut ) {_fOut = fOut;}
  void fits(bool mcTruth,bool isMC,std::string title = "", bool isaddGaus=false);
  void useMinos(bool minos = true) {_useMinos = minos;}
  void textParForCanvas(RooFitResult *resP, RooFitResult *resF, TPad *p);
  
  void fixSigmaFtoSigmaP(bool fix=true) { _fixSigmaFtoSigmaP= fix;}
  void fixBkgPassToFail(bool fix=true) { _fixBkgPassToFail = fix; }  // new method

  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut;
  double _nTotP, _nTotF;
  bool _useMinos;
  bool _fixSigmaFtoSigmaP;
  bool _fixBkgPassToFail = false;  // new member variable
  double _xFitMin,_xFitMax;
  int _nBins = 10000;
  void initFromHistograms(TH1* hPassOrig, TH1* hFailOrig, const std::string& histname);
  void zeroBinsOutsideRange(TH1* hPass, TH1* hFail, double low, double high);
};

tnpFitter::tnpFitter(TFile *filein, std::string histname) {
    TH1* hPass = static_cast<TH1*>(filein->Get(TString::Format("%s_Pass", histname.c_str()).Data()));
    TH1* hFail = static_cast<TH1*>(filein->Get(TString::Format("%s_Fail", histname.c_str()).Data()));
    initFromHistograms(hPass, hFail, histname);
}

// Constructor from TH1* — unchanged signature
tnpFitter::tnpFitter(TH1 *hPass, TH1 *hFail, std::string histname) {
    initFromHistograms(hPass, hFail, histname);
}

// Private common initializer
void tnpFitter::initFromHistograms(TH1* hPassOrig, TH1* hFailOrig, const std::string& histname) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    _histname_base = histname;
    _xFitMin = 60.0;
    _xFitMax = 120.0;
    _useMinos = false;
    _fixSigmaFtoSigmaP = false;

    // Clone to avoid modifying input histograms
    TH1* hPass = static_cast<TH1*>(hPassOrig->Clone());
    TH1* hFail = static_cast<TH1*>(hFailOrig->Clone());

    _nTotP = hPass->Integral();
    _nTotF = hFail->Integral();

    _work = new RooWorkspace("w");
    _work->factory("x[50,130]");

    // Full-range datasets
    RooDataHist rooPass("hPassFull", "hPassFull", *_work->var("x"), hPass);
    RooDataHist rooFail("hFailFull", "hFailFull", *_work->var("x"), hFail);
    _work->import(rooPass);
    _work->import(rooFail);

    // Fit dataset: [60, 120]
    zeroBinsOutsideRange(hPass, hFail, 60.0, 120.0);
    RooDataHist rooPassFit("hPass", "hPass", *_work->var("x"), hPass);
    RooDataHist rooFailFit("hFail", "hFail", *_work->var("x"), hFail);
    _work->import(rooPassFit);
    _work->import(rooFailFit);

    delete hPass;
    delete hFail;
}

// Helper: zero bins outside [low, high]
void tnpFitter::zeroBinsOutsideRange(TH1* hPass, TH1* hFail, double low, double high) {
    int nBins = hPass->GetXaxis()->GetNbins();
    for (int ib = 0; ib <= nBins + 1; ++ib) {
        double center = hPass->GetXaxis()->GetBinCenter(ib);
        if (center <= low || center >= high) {
            hPass->SetBinContent(ib, 0.0);
            hFail->SetBinContent(ib, 0.0);
        }
    }
}


void tnpFitter::setZLineShapes(TH1 *hZPass, TH1 *hZFail, bool zeroLowerThreshold, double threshold) {
  if (zeroLowerThreshold) {
    zeroBinsOutsideRange(hZPass, hZFail, threshold, 120.0);
  }
  RooDataHist rooPass("hGenZPass","hGenZPass",*_work->var("x"),hZPass);
  RooDataHist rooFail("hGenZFail","hGenZFail",*_work->var("x"),hZFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;  
}

void tnpFitter::setWorkspace(std::vector<std::string> workspace, bool isaddGaus, bool useBreitWigner) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }

  _work->var("x")->setBins(_nBins, "cache");
  if (useBreitWigner) {
      cout << "Using Breit-Wigner for Z line shape" << endl;
      _work->factory("mZ[91.1876]");
      _work->factory("widthZ[2.4952]");
      _work->var("mZ")->setConstant(kTRUE);
      _work->var("widthZ")->setConstant(kTRUE);

      // using Breit-Wigner
      _work->factory("BreitWigner::sigPhysPass(x, mZ, widthZ)");
      _work->factory("BreitWigner::sigPhysFail(x, mZ, widthZ)");
  } 
  else {
      _work->factory("HistPdf::sigPhysPass(x,hGenZPass,3)");
      _work->factory("HistPdf::sigPhysFail(x,hGenZFail,3)");
  }
  _work->factory("FCONV::sigPass(x, sigPhysPass , sigResPass)");
  _work->factory("FCONV::sigFail(x, sigPhysFail , sigResFail)");
  _work->factory(TString::Format("nSigP[%f,1e-10,%f]",_nTotP*0.9,_nTotP*1.5)); // the default effi will be 1e-4
  _work->factory(TString::Format("nBkgP[%f,1e-6,%f]",_nTotP*0.1,_nTotP*1.5));
  _work->factory(TString::Format("nSigF[%f,1e-6,%f]",_nTotF*0.9,_nTotF*1.5));
  _work->factory(TString::Format("nBkgF[%f,1e-6,%f]",_nTotF*0.1,_nTotF*1.5));
  _work->factory("SUM::pdfPass(nSigP*sigPass,nBkgP*bkgPass)");
  
  if (isaddGaus) {
    _work->factory("SUM::pdfFail(expr('sigFracF*nSigF',{sigFracF,nSigF})*sigFail,nBkgF*bkgFail, expr('(1.-sigFracF)*nSigF',{sigFracF,nSigF})*sigGaussFail)");
  } 
  else {
    _work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
  }
  _work->Print();			         
}

void tnpFitter::fits(bool mcTruth,bool isMC,string title, bool isaddGaus) {

  cout << " title : " << title << endl;

  RooRealVar *x = _work->var("x");
  x->setRange("fitMassRange",_xFitMin,_xFitMax);
  RooAbsPdf *pdfPass = _work->pdf("pdfPass");
  RooAbsPdf *pdfFail = _work->pdf("pdfFail");
  RooAbsPdf *bkgPass = _work->pdf("bkgPass");
  RooAbsPdf *bkgFail = _work->pdf("bkgFail");
  RooFitResult* resPass;  
  RooFitResult* resFail;

  RooDataHist *dataPassFull = (RooDataHist*)_work->data("hPassFull");
  RooDataHist *dataFailFull = (RooDataHist*)_work->data("hFailFull");
  RooDataHist *dataPass = (RooDataHist*)_work->data("hPass");
  RooDataHist *dataFail = (RooDataHist*)_work->data("hFail");

  if( mcTruth ) {
    _work->var("nBkgP")->setVal(0); _work->var("nBkgP")->setConstant();
    _work->var("nBkgF")->setVal(0); _work->var("nBkgF")->setConstant();
    if( _work->var("sosP")   ) { 
      // the minimum value in config is 0.5, so it can not be setted to 0 directly
      _work->var("sosP")->setRange(0,0.001);
      _work->var("sosP")->setVal(0);
      _work->var("sosP")->setConstant(); }
    if( _work->var("sosF")   ) { 
      _work->var("sosF")->setRange(0,0.001);
      _work->var("sosF")->setVal(0);
      _work->var("sosF")->setConstant(); }
    if( _work->var("acmsP")  ) _work->var("acmsP")->setConstant();
    if( _work->var("acmsF")  ) _work->var("acmsF")->setConstant();
    if( _work->var("betaP")  ) _work->var("betaP")->setConstant();
    if( _work->var("betaF")  ) _work->var("betaF")->setConstant();
    if( _work->var("gammaP") ) _work->var("gammaP")->setConstant();
    if( _work->var("gammaF") ) _work->var("gammaF")->setConstant();
  }

    // first perform a simple fit without Minos to get close to minimum
  bool doPreFit = true;
  
  // Helper lambda to copy background parameters from Fail to Pass
  auto copyBkgParamsFailToPass = [&]() {
    // List of background parameter pairs (Fail -> Pass)
    std::vector<std::pair<std::string, std::string>> bkgParamPairs = {
      {"acmsF", "acmsP"},
      {"betaF", "betaP"},
      {"gammaF", "gammaP"},
      {"peakF", "peakP"},
    };
    for (const auto& pair : bkgParamPairs) {
      RooRealVar* varF = _work->var(pair.first.c_str());
      RooRealVar* varP = _work->var(pair.second.c_str());
      if (varF && varP) {
        double valF = varF->getVal();
        varP->setVal(valF);
        cout << "  Set " << pair.second << " = " << valF << " (from " << pair.first << ")" << endl;
        if (_fixBkgPassToFail) {
          varP->setConstant(kTRUE);
          cout << "  Fixed " << pair.second << " to Fail value" << endl;
        }
      }
    }
  };

  // --- Fit Fail first, then Pass ---
  cout << "--- Fitting Fail category first... ---" << endl;
  
  if( doPreFit ) {
    cout << "--- Performing pre fit (no Minos) to stabilize parameters... ---" << endl;
    if( isMC ) {
      pdfFail->fitTo(*dataFail, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("fitMassRange"));
    } else {
      pdfFail->fitTo(*dataFail, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("fitMassRange"));
    }
  }

  if( isMC ) {
      resFail = pdfFail->fitTo(*dataFail, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("fitMassRange"));
  } else {
      resFail = pdfFail->fitTo(*dataFail, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("fitMassRange"));
  }
  
  cout << "Global fit results (failing):" << endl;
  resFail->Print("v");

  // Copy background parameters from Fail to Pass (set initial values or fix)
  cout << "--- Copying background parameters from Fail to Pass... ---" << endl;
  copyBkgParamsFailToPass();

  // --- Now fit Pass category ---
  cout << "--- Fitting Pass category... ---" << endl;
  
  if( doPreFit ) {
    if( isMC ) {
      pdfPass->fitTo(*dataPass, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("fitMassRange"));
    } else {
      pdfPass->fitTo(*dataPass, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("fitMassRange"));
    }
  }

  if( isMC ) {
      resPass = pdfPass->fitTo(*dataPass, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("fitMassRange"));
  } else {
      resPass = pdfPass->fitTo(*dataPass, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("fitMassRange"));
  }

  // print fit results
  cout << "Global fit results (passing):" << endl;
  resPass->Print("v");
  cout << "--- Global fit finished. ---" << endl;

  // reset range for plotting
  _work->var("x")->setRange(_xFitMin,_xFitMax);


  RooPlot *pPass = _work->var("x")->frame(60,120);
  RooPlot *pFail = _work->var("x")->frame(60,120);
  pPass->SetTitle("passing probe");
  pFail->SetTitle("failing probe");
  
  _work->data("hPass") ->plotOn( pPass );
  _work->pdf("pdfPass")->plotOn( pPass, LineColor(kRed) );
  _work->pdf("pdfPass")->plotOn( pPass, Components("bkgPass"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hPass") ->plotOn( pPass );
  
  _work->data("hFail") ->plotOn( pFail );
  _work->pdf("pdfFail")->plotOn( pFail, LineColor(kRed) );
  _work->pdf("pdfFail")->plotOn( pFail, Components("bkgFail"),LineColor(kBlue),LineStyle(kDashed));
  _work->data("hFail") ->plotOn( pFail );

  TCanvas c("c","c",1100,450);
  c.Divide(3,1);
  TPad *padText = (TPad*)c.GetPad(1);
  textParForCanvas( resPass,resFail, padText );
  c.cd(2); pPass->Draw();
  c.cd(3); pFail->Draw();

  _fOut->cd();
  c.Write(TString::Format("%s_Canv",_histname_base.c_str()),TObject::kOverwrite);
  resPass->Write(TString::Format("%s_resP",_histname_base.c_str()),TObject::kOverwrite);
  resFail->Write(TString::Format("%s_resF",_histname_base.c_str()),TObject::kOverwrite);

  
}





/////// Stupid parameter dumper /////////
void tnpFitter::textParForCanvas(RooFitResult *resP, RooFitResult *resF,TPad *p) {

  double eff = -1;
  double e_eff = 0;

  RooRealVar *nSigP = _work->var("nSigP");
  RooRealVar *nSigF = _work->var("nSigF");
  
  double nP   = nSigP->getVal();
  double e_nP = nSigP->getError();
  double nF   = nSigF->getVal();
  double e_nF = nSigF->getError();
  double nTot = nP+nF;
  eff = nP / (nP+nF);
  e_eff = 1./(nTot*nTot) * sqrt( nP*nP* e_nF*e_nF + nF*nF * e_nP*e_nP );

  TPaveText *text1 = new TPaveText(0,0.8,1,1);
  text1->SetFillColor(0);
  text1->SetBorderSize(0);
  text1->SetTextAlign(12);

  text1->AddText(TString::Format("* fit status pass: %d, fail : %d",resP->status(),resF->status()));
  text1->AddText(TString::Format("* eff = %1.4f #pm %1.4f",eff,e_eff));

  //  text->SetTextSize(0.06);

//  text->AddText("* Passing parameters");
  TPaveText *text = new TPaveText(0,0,1,0.8);
  text->SetFillColor(0);
  text->SetBorderSize(0);
  text->SetTextAlign(12);
  text->AddText("    --- parmeters " );
  RooArgList listParFinalP = resP->floatParsFinal();
  for( int ip = 0; ip < listParFinalP.getSize(); ip++ ) {
    TString vName = listParFinalP[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }

//  text->AddText("* Failing parameters");
  RooArgList listParFinalF = resF->floatParsFinal();
  for( int ip = 0; ip < listParFinalF.getSize(); ip++ ) {
    TString vName = listParFinalF[ip].GetName();
    text->AddText(TString::Format("   - %s \t= %1.3f #pm %1.3f",
				  vName.Data(),
				  _work->var(vName)->getVal(),
				  _work->var(vName)->getError() ) );
  }

  p->cd();
  text1->Draw();
  text->Draw();
}
