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
  void setZLineShapes(TH1 *hZPass, TH1 *hZFail );
  void setWorkspace(std::vector<std::string>, bool isaddGaus=false);
  void setOutputFile(TFile *fOut ) {_fOut = fOut;}
  void fits(bool mcTruth,bool isMC,std::string title = "", bool isaddGaus=false);
  void useMinos(bool minos = true) {_useMinos = minos;}
  void textParForCanvas(RooFitResult *resP, RooFitResult *resF, TPad *p);
  
  void fixSigmaFtoSigmaP(bool fix=true) { _fixSigmaFtoSigmaP= fix;}

  void setFitRange(double xMin,double xMax) { _xFitMin = xMin; _xFitMax = xMax; }
private:
  RooWorkspace *_work;
  std::string _histname_base;
  TFile *_fOut;
  double _nTotP, _nTotF;
  bool _useMinos;
  bool _fixSigmaFtoSigmaP;
  double _xFitMin,_xFitMax;
  int _nBins = 10000;
};

tnpFitter::tnpFitter(TFile *filein, std::string histname   ) : _useMinos(false),_fixSigmaFtoSigmaP(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  _histname_base = histname;  

  TH1 *hPass = (TH1*) filein->Get(TString::Format("%s_Pass",histname.c_str()).Data());
  TH1 *hFail = (TH1*) filein->Get(TString::Format("%s_Fail",histname.c_str()).Data());
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  // for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
  //  if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
  //    hPass->SetBinContent(ib,0);
  //    hFail->SetBinContent(ib,0);
  //  }
  
  _work = new RooWorkspace("w") ;
  _work->factory("x[50,130]");

  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;

  // create datasets for fitting
  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
   if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
     hPass->SetBinContent(ib,0);
     hFail->SetBinContent(ib,0);
   }
  RooDataHist rooPassFit("hPassFit","hPassFit",*_work->var("x"),hPass);
  RooDataHist rooFailFit("hFailFit","hFailFit",*_work->var("x"),hFail);
  _work->import(rooPassFit) ;
  _work->import(rooFailFit) ;

  _xFitMin = 60;
  _xFitMax = 120;
}

tnpFitter::tnpFitter(TH1 *hPass, TH1 *hFail, std::string histname  ) : _useMinos(false),_fixSigmaFtoSigmaP(false) {
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  _histname_base = histname;
  
  _nTotP = hPass->Integral();
  _nTotF = hFail->Integral();
  /// MC histos are done between 50-130 to do the convolution properly
  /// but when doing MC fit in 60-120, need to zero bins outside the range
  // for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
  //   if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
  //     hPass->SetBinContent(ib,0);
  //     hFail->SetBinContent(ib,0);
  //   }
  
  _work = new RooWorkspace("w") ;
  _work->factory("x[50,130]");
  
  RooDataHist rooPass("hPass","hPass",*_work->var("x"),hPass);
  RooDataHist rooFail("hFail","hFail",*_work->var("x"),hFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;

  for( int ib = 0; ib <= hPass->GetXaxis()->GetNbins()+1; ib++ )
   if(  hPass->GetXaxis()->GetBinCenter(ib) <= 60 || hPass->GetXaxis()->GetBinCenter(ib) >= 120 ) {
     hPass->SetBinContent(ib,0);
     hFail->SetBinContent(ib,0);
   }
  RooDataHist rooPassFit("hPassFit","hPassFit",*_work->var("x"),hPass);
  RooDataHist rooFailFit("hFailFit","hFailFit",*_work->var("x"),hFail);
  _work->import(rooPassFit) ;
  _work->import(rooFailFit) ;
  
  _xFitMin = 60;
  _xFitMax = 120;
  
}


void tnpFitter::setZLineShapes(TH1 *hZPass, TH1 *hZFail ) {
  RooDataHist rooPass("hGenZPass","hGenZPass",*_work->var("x"),hZPass);
  RooDataHist rooFail("hGenZFail","hGenZFail",*_work->var("x"),hZFail);
  _work->import(rooPass) ;
  _work->import(rooFail) ;  
}

void tnpFitter::setWorkspace(std::vector<std::string> workspace, bool isaddGaus) {
  for( unsigned icom = 0 ; icom < workspace.size(); ++icom ) {
    _work->factory(workspace[icom].c_str());
  }

  _work->var("x")->setBins(_nBins, "cache");
  _work->factory("HistPdf::sigPhysPass(x,hGenZPass,3)");
  _work->factory("HistPdf::sigPhysFail(x,hGenZFail,3)");
  _work->factory("FCONV::sigPass(x, sigPhysPass , sigResPass)");
  _work->factory("FCONV::sigFail(x, sigPhysFail , sigResFail)");
  _work->factory(TString::Format("nSigP[%f,0.5,%f]",_nTotP*0.9,_nTotP*1.5));
  _work->factory(TString::Format("nBkgP[%f,0.5,%f]",_nTotP*0.1,_nTotP*1.5));
  _work->factory(TString::Format("nSigF[%f,0.5,%f]",_nTotF*0.9,_nTotF*1.5));
  _work->factory(TString::Format("nBkgF[%f,0.5,%f]",_nTotF*0.1,_nTotF*1.5));
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

  RooAbsPdf *pdfPass = _work->pdf("pdfPass");
  RooAbsPdf *pdfFail = _work->pdf("pdfFail");
  RooAbsPdf *bkgPass = _work->pdf("bkgPass");
  RooAbsPdf *bkgFail = _work->pdf("bkgFail");
  
  // using full dataset for sideband fit
  RooDataHist *dataPassFull = (RooDataHist*)_work->data("hPass");
  RooDataHist *dataFailFull = (RooDataHist*)_work->data("hFail");

  if( mcTruth ) {
    _work->var("nBkgP")->setVal(0); _work->var("nBkgP")->setConstant();
    _work->var("nBkgF")->setVal(0); _work->var("nBkgF")->setConstant();
    if( _work->var("sosP")   ) { _work->var("sosP")->setVal(0);
      _work->var("sosP")->setConstant(); }
    if( _work->var("sosF")   ) { _work->var("sosF")->setVal(0);
      _work->var("sosF")->setConstant(); }
    if( _work->var("acmsP")  ) _work->var("acmsP")->setConstant();
    if( _work->var("acmsF")  ) _work->var("acmsF")->setConstant();
    if( _work->var("betaP")  ) _work->var("betaP")->setConstant();
    if( _work->var("betaF")  ) _work->var("betaF")->setConstant();
    if( _work->var("gammaP") ) _work->var("gammaP")->setConstant();
    if( _work->var("gammaF") ) _work->var("gammaF")->setConstant();
  } else {

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

  // --- using dataset for full range fit---
  cout << "--- Creating datasets for full range fit (60-120)... ---" << endl;
  
  RooDataHist *dataPassFit = (RooDataHist*)_work->data("hPassFit");
  RooDataHist *dataFailFit = (RooDataHist*)_work->data("hFailFit");

  _work->var("x")->setRange(_xFitMin,_xFitMax);

  // setting initial values for signal shape parameters
  if (_work->var("sigmaP")) {
      _work->var("sigmaP")->setVal(2.0);
      _work->var("sigmaP")->setRange(0.5, 5.0);
  }
  if (_work->var("sigmaF")) {
      _work->var("sigmaF")->setVal(2.0);
      _work->var("sigmaF")->setRange(0.5, 5.0);
  }
  if (_work->var("meanP")) _work->var("meanP")->setVal(0);
  if (_work->var("meanF")) _work->var("meanF")->setVal(0);

  // performing global signal + background fit on full range (60-120)
  cout << "--- Performing global signal + background fit on full range (60-120)... ---" << endl;
  RooFitResult* resPass;  
  RooFitResult* resFail;

  if( _fixSigmaFtoSigmaP ) {
    _work->var("sigmaF")->setVal( _work->var("sigmaP")->getVal() );
    _work->var("sigmaF")->setConstant();
  }

  if( isMC ) {
      resPass = pdfPass->fitTo(*dataPassFit, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("full"));
      resFail = pdfFail->fitTo(*dataFailFit, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kTRUE), Save(), Range("full"));
  } else {
      resPass = pdfPass->fitTo(*dataPassFit, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kFALSE), Save(), Range("full"));
      resFail = pdfFail->fitTo(*dataFailFit, Minimizer("Minuit2", "MIGRAD"), Minos(_useMinos), Strategy(2), SumW2Error(kFALSE), Save(), Range("full"));
  }
  cout << "--- Global fit finished. ---" << endl;

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