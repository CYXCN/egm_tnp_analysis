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
  tnpFitter(TFile *file, std::string histname, double xMcMin = 60.0);
  tnpFitter(TH1 *hPass, TH1 *hFail, std::string histname, double xMcMin = 60.0);
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
  void setMCFitLow(double xMcMin) { _xMcMin = xMcMin; }
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut;
  double _nTotP, _nTotF;
  bool _useMinos;
  bool _fixSigmaFtoSigmaP;
  bool _fixBkgPassToFail = false;  // new member variable
  double _xFitMin,_xFitMax;
  double _xMcMin = 60;
  int _nBins = 10000;
  void initFromHistograms(TH1* hPassOrig, TH1* hFailOrig, const std::string& histname, double xMcMin = 60.0);
  void zeroBinsOutsideRange(TH1* hPass, TH1* hFail, double low, double high);
};

tnpFitter::tnpFitter(TFile *filein, std::string histname, double xMcMin) {
    TH1* hPass = static_cast<TH1*>(filein->Get(TString::Format("%s_Pass", histname.c_str()).Data()));
    TH1* hFail = static_cast<TH1*>(filein->Get(TString::Format("%s_Fail", histname.c_str()).Data()));
    initFromHistograms(hPass, hFail, histname, xMcMin); // use provided xMcMin
}

// Constructor from TH1* — unchanged signature
tnpFitter::tnpFitter(TH1 *hPass, TH1 *hFail, std::string histname, double xMcMin) {
    initFromHistograms(hPass, hFail, histname, xMcMin); // use provided xMcMin
}

// Private common initializer
void tnpFitter::initFromHistograms(TH1* hPassOrig, TH1* hFailOrig, const std::string& histname, double xMcMin) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    _histname_base = histname;
    _xMcMin = xMcMin;        // store for later use if needed
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
    RooDataHist rooPass("hPass", "hPass", *_work->var("x"), hPass);
    RooDataHist rooFail("hFail", "hFail", *_work->var("x"), hFail);
    _work->import(rooPass);
    _work->import(rooFail);

    // Fit dataset: [60, 120]
    zeroBinsOutsideRange(hPass, hFail, 60.0, 120.0);
    RooDataHist rooPassFit("hPassFit", "hPassFit", *_work->var("x"), hPass);
    RooDataHist rooFailFit("hFailFit", "hFailFit", *_work->var("x"), hFail);
    _work->import(rooPassFit);
    _work->import(rooFailFit);

    // MC fit dataset: [_xMcMin, 120]
    zeroBinsOutsideRange(hPass, hFail, _xMcMin, 120.0);
    RooDataHist rooPassFitMC("hPassFitMC", "hPassFitMC", *_work->var("x"), hPass);
    RooDataHist rooFailFitMC("hFailFitMC", "hFailFitMC", *_work->var("x"), hFail);
    _work->import(rooPassFitMC);
    _work->import(rooFailFitMC);

    // Core dataset: [85, 95]
    zeroBinsOutsideRange(hPass, hFail, 85.0, 95.0);
    RooDataHist rooPassCore("hPassCore", "hPassCore", *_work->var("x"), hPass);
    RooDataHist rooFailCore("hFailCore", "hFailCore", *_work->var("x"), hFail);
    _work->import(rooPassCore);
    _work->import(rooFailCore);

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
  
  // define（Sidebands）and（Full range）
  // sideband for bkg fitting
  x->setRange("lowSide", 50, 75);        
  x->setRange("highSide", 105, 130);    
  x->setRange("full", _xFitMin, _xFitMax); // global fit range : 60-120 GeV
  x->setRange("core", 85, 95);          // signal core region

  RooAbsPdf *pdfPass = _work->pdf("pdfPass");
  RooAbsPdf *pdfFail = _work->pdf("pdfFail");
  RooAbsPdf *bkgPass = _work->pdf("bkgPass");
  RooAbsPdf *bkgFail = _work->pdf("bkgFail");
  
  // using full dataset for sideband fit
  RooDataHist *dataPassFull = (RooDataHist*)_work->data("hPass");
  RooDataHist *dataFailFull = (RooDataHist*)_work->data("hFail");

  bool do_bkgPrefit = true;

  // fix peak as this parameter is the same as bkg number
  if( _work->var("peakP")  ) _work->var("peakP")->setConstant();
  if( _work->var("peakF")  ) _work->var("peakF")->setConstant();

  if( mcTruth ) {
    cout << "--- Performing MC truth fit (no background) ... ---" << endl;
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
  } else {
    
    if ( do_bkgPrefit ) {
      cout << "--- Resetting background yields for background pre-fit... ---" << endl;

      // bkg prefit only in sidebands
      double nEventsPassLow = dataPassFull->sumEntries("1", "lowSide");
      double nEventsFailLow = dataFailFull->sumEntries("1", "highSide");

      double nEventsPassHigh = dataPassFull->sumEntries("1", "highSide");
      double nEventsFailHigh = dataFailFull->sumEntries("1", "highSide");

      int thresholdEvents = 10; // minimum number of events required in each sideband for fitting
      bool doPassFit = nEventsPassLow > thresholdEvents && nEventsPassHigh > thresholdEvents;
      bool doFailFit = nEventsFailLow > thresholdEvents && nEventsFailHigh > thresholdEvents;
      cout << "--- Performing background-only fit on sidebands (50-60, 120-130)... ---" << endl;
      
      if( !doPassFit ) {
        cout << "!!! Not enough events in passing sidebands for background fit. Skipping background fit for passing category. !!!" << endl;
      }
      else {
        bkgPass->fitTo(*dataPassFull, Minimizer("Minuit2", "MIGRAD"), Strategy(2), SumW2Error(isMC), Range("lowSide,highSide"), PrintLevel(1));
        cout << "Fitted background parameters (passing):" << endl;
        RooArgSet *bkgPassParams = bkgPass->getParameters(*x);
        bkgPassParams->Print("v");
      }
      if( !doFailFit ) {
        cout << "!!! Not enough events in failing sidebands for background fit. Skipping background fit for failing category. !!!" << endl;
      }
      else {
        bkgFail->fitTo(*dataFailFull, Minimizer("Minuit2", "MIGRAD"), Strategy(2), SumW2Error(isMC), Range("lowSide,highSide"), PrintLevel(1));
        cout << "Fitted background parameters (failing):" << endl;
        RooArgSet *bkgFailParams = bkgFail->getParameters(*x);
        bkgFailParams->Print("v");
      }
      
      cout << "--- Background-only fit finished. ---" << endl;
    }
  }

  // --- using dataset for full range fit---
  cout << "--- Loading datasets for full range fit (60-120)... ---" << endl;
  
  RooDataHist *dataPassFit = (RooDataHist*)_work->data("hPassFit");
  RooDataHist *dataFailFit = (RooDataHist*)_work->data("hFailFit");
  RooDataHist *dataPassFitMC = nullptr;
  RooDataHist *dataFailFitMC = nullptr;

  _work->var("x")->setRange(_xFitMin,_xFitMax);

  if( isMC ) {
    cout << "--- Using MC fit range (70-120) ---" << endl;
    dataPassFitMC = (RooDataHist*)_work->data("hPassFitMC");
    dataFailFitMC = (RooDataHist*)_work->data("hFailFitMC");
    _work->var("x")->setRange(70,_xFitMax);
  }


  // setting initial values for signal shape parameters
  // if (_work->var("sigmaP")) {
  //     _work->var("sigmaP")->setVal(2.0);
  //     _work->var("sigmaP")->setRange(0.5, 5.0);
  // }
  // if (_work->var("sigmaF")) {
  //     _work->var("sigmaF")->setVal(2.0);
  //     _work->var("sigmaF")->setRange(0.5, 5.0);
  // }
  auto find_max_bin_center = [&](RooDataHist* data) {
    double max_weight = -1.0;
    double max_center = -999.0;
    for (int i = 0; i < data->numEntries(); ++i) {
      const RooArgSet* row = data->get(i);
      if (data->weight() > max_weight) {
        max_weight = data->weight();
        max_center = ((RooRealVar*)row->find(x->GetName()))->getVal();
      }
    }
    return max_center;
  };
  // double shiftP = find_max_bin_center(dataPassFull) - 91.1876;
  // double shiftF = find_max_bin_center(dataFailFull) - 91.1876;

  // if (_work->var("meanP")) _work->var("meanP")->setVal(0);
  // if (_work->var("meanF")) _work->var("meanF")->setVal(0);
  // if (_work->var("meanP")) _work->var("meanP")->setVal(mcTruth ? 0 : shiftP);
  // if (_work->var("meanF")) _work->var("meanF")->setVal(mcTruth ? 0 : shiftF);
  bool doMCInitFit = false; // do initial fit only for MC
  if (isMC && doMCInitFit) {
      cout << "--- Performing MC initial fit in core region (85, 95)... ---" << endl;

      // Helper lambda to perform the core fit for Pass or Fail
      auto fitMCCore = [&](RooAbsPdf* pdf, RooDataHist* data, string suffix) {
          // Identify variables
          RooRealVar* v_mean    = _work->var(("mean" + suffix).c_str());
          RooRealVar* v_sigma   = _work->var(("sigma" + suffix).c_str());
          RooRealVar* v_sigma_2 = _work->var(("sigma" + suffix + "_2").c_str());
          
          RooRealVar* v_n    = _work->var(("n" + suffix).c_str());
          RooRealVar* v_alpha   = _work->var(("alpha" + suffix).c_str());
          RooRealVar* v_sos     = _work->var(("sos" + suffix).c_str());

          _work->var("x")->setRange(85,95);

          // 1. Fix nSig, alpha, sos (Yields and Tails)
          if (v_n)  v_n->setConstant(kTRUE);
          if (v_alpha) v_alpha->setConstant(kTRUE);
          if (v_sos)   v_sos->setConstant(kTRUE);

          // 2. Float mean, sigma, sigma_2 (Core Shape)
          if (v_mean)    v_mean->setConstant(kFALSE);
          if (v_sigma)   v_sigma->setConstant(kFALSE);
          if (v_sigma_2) v_sigma_2->setConstant(kFALSE);

          // 3. Fit in Core Region
          // Use SumW2Error(kTRUE) because it is MC
          pdf->fitTo(*data, Minimizer("Minuit2", "MIGRAD"), Strategy(2), 
                     SumW2Error(kTRUE), Range("core"), PrintLevel(1));

          // 4. Release Fixed Parameters for Global Fit
          // We must release nSig so yield can be calculated.
          // We also release alpha/sos to allow global fit to adjust tails if needed.
          if (v_n)  v_n->setConstant(kFALSE);
          if (v_alpha) v_alpha->setConstant(kFALSE);
          if (v_sos)   v_sos->setConstant(kFALSE);

          // only allow 10% variation from core fit values
          if (v_mean) {
              double val = v_mean->getVal();
              if (val > 0)
                v_mean->setRange(val*0.8, val*1.2);
              else
                v_mean->setRange(val*1.2, val*0.8);
          }
          if (v_sigma) {
              double val = v_sigma->getVal();
              v_sigma->setRange(val*0.8, val*1.2);
          }
          if (v_sigma_2) {
              double val = v_sigma_2->getVal();
              v_sigma_2->setRange(val*0.8, val*1.2);
          }
          
          RooArgSet *params = pdf->getParameters(*x);
          cout << "MC core fit parameters (" << (suffix == "P" ? "passing" : "failing") << "):" << endl;
          params->Print("v");
          _work->var("x")->setRange(_xFitMin,_xFitMax);
      };
      RooDataHist *dataPassCore = (RooDataHist*)_work->data("hPassCore");
      RooDataHist *dataFailCore = (RooDataHist*)_work->data("hFailCore");
      // Perform fit for Pass
      fitMCCore(pdfPass, dataPassCore, "P");
      // Perform fit for Fail
      fitMCCore(pdfFail, dataFailCore, "F");
      
      cout << "--- MC initial fit finished. ---" << endl;
  }
  // =========================================================================

  // performing global signal + background fit on full range (60-120)
  cout << "--- Performing global signal + background fit on full range (60-120)... ---" << endl;
  RooFitResult* resPass;  
  RooFitResult* resFail;

  if( _fixSigmaFtoSigmaP ) {
    _work->var("sigmaF")->setVal( _work->var("sigmaP")->getVal() );
    _work->var("sigmaF")->setConstant();
  }

  // first perform a simple fit without Minos to get close to minimum
  bool do_simpleFit = true;
  
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
  
  if( do_simpleFit ) {
    cout << "--- Performing simple fit (no Minos) to stabilize parameters... ---" << endl;
    if( isMC ) {
      pdfFail->fitTo(*dataFailFitMC, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("full"));
    } else {
      pdfFail->fitTo(*dataFailFit, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("full"));
    }
  }

  if( isMC ) {
      resFail = pdfFail->fitTo(*dataFailFitMC, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("full"));
  } else {
      resFail = pdfFail->fitTo(*dataFailFit, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("full"));
  }
  
  cout << "Global fit results (failing):" << endl;
  resFail->Print("v");

  // Copy background parameters from Fail to Pass (set initial values or fix)
  cout << "--- Copying background parameters from Fail to Pass... ---" << endl;
  copyBkgParamsFailToPass();

  // --- Now fit Pass category ---
  cout << "--- Fitting Pass category... ---" << endl;
  
  if( do_simpleFit ) {
    if( isMC ) {
      pdfPass->fitTo(*dataPassFitMC, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("full"));
    } else {
      pdfPass->fitTo(*dataPassFit, Minimizer("Minuit2", "MIGRAD"), Minos(kFALSE), Strategy(0), SumW2Error(kTRUE), Range("full"));
    }
  }

  if( isMC ) {
      resPass = pdfPass->fitTo(*dataPassFitMC, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("full"));
  } else {
      resPass = pdfPass->fitTo(*dataPassFit, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("full"));
  }

  // print fit results
  cout << "Global fit results (passing):" << endl;
  resPass->Print("v");
  cout << "--- Global fit finished. ---" << endl;

  // reset range for plotting
  _work->var("x")->setRange(_xFitMin,_xFitMax);

  RooPlot *pPass = x->frame(Title("passing probe"), Range("full"));
  RooPlot *pFail = x->frame(Title("failing probe"), Range("full"));
  
  dataPassFit->plotOn( pPass );
  pdfPass->plotOn( pPass, LineColor(kRed) );
  pdfPass->plotOn( pPass, Components(*bkgPass), LineColor(kBlue), LineStyle(kDashed));
  dataPassFit->plotOn( pPass );
  
  dataFailFit->plotOn( pFail );
  pdfFail->plotOn( pFail, LineColor(kRed) );
  pdfFail->plotOn( pFail, Components(*bkgFail), LineColor(kBlue), LineStyle(kDashed));
  dataFailFit->plotOn( pFail );

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