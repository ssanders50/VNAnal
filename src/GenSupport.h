#ifndef GENSUPPORT
#define GENSUPPORT

#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include <iostream>
#include <complex>
#include <cmath>
#include "Qvec.h"
#include "PixelTrackEfficiency.h"

typedef complex<double> comp;
static const int ERRCNT = 10;
static const int NumPtBins = 28;
class GenSupport{
 public:
  explicit GenSupport(int eporder, GetEventInfo * info,  
int ep2a, int ep2b, int ep2c,int ep3a, int ep3b, int ep3c, int mintrk, int maxtrk,int subevt);
  void flip();
  void add(TH2D * qxtrk, TH2D * qytrk, TH2D * qcnt, double * qx, double* qy, double * qw, TH2D * avpt, int noff,TH2D * qxav, TH2D * qyav);
  TH2D * getRescor();
  TH1D * getRescorAv();
  TH2D * getEPVN();
  TH2D * getEPVNobs();
  TH2D * getSPVN();
  TH2D * getSpecRaw();
  TH2D * getSpecEffCorrected();
  TH2D * getSPEPRatio();
  TH2D * getQAB(){return qab;}
  TH2D * getQAC(){return qac;}
  TH2D * getQBC(){return qbc;}
  TH2D * getQnA(){return qna;}
  TH1D * getEPInt(double minpt, double maxpt);
  TH1D * getSPInt(double minpt, double maxpt);
  TH1D * getSPEPIntRatio(double minpt, double maxpt);

  TH2D * getQnFraction();

  int getCount(){return count;}
  int getMinTrk(){return mintrk_;}
  int getMaxTrk(){return maxtrk_;}
  TH2D * getAvPt();
  TH2D * getAvNoff();
  TH1D * getVnObs(int ptbin){return vnobs[ptbin];}
  TH2D * getVnObs2D(int ptbin){return vnobs2D[ptbin];}

 private:
  double avcent_;
  TString AName2;
  TString BName2;
  TString CName2;
  TString AName3;
  TString BName3;
  TString CName3;
  TH2D * heffa;
  int ep2a_;
  int ep2b_;
  int ep2c_;
  int ep2aSet_;
  int ep2bSet_;
  int ep2cSet_;
  int ep3a_;
  int ep3b_;
  int ep3c_;
  int ep3aSet_;
  int ep3bSet_;
  int ep3cSet_;
  TH1D * ptbins;
  TH1D * etabins;
  double epord_;
  std::vector < std::vector<Qvec *> > q;
  std::vector < std::vector<Qvec *> > qerr[ERRCNT];
  Qvec *qav;
  Qvec *qaverr[ERRCNT];

  TH2D * nfrac;
  TH2D * rescor;
  TH2D * epvn;
  TH2D * epvnobs;
  TH2D * spvn;
  TH2D * specRaw;
  TH2D * specEffCorrected;
  TH1D * epInt;
  TH1D * spInt;
  TH1D * spepintratio;
  TH2D * avpt_;
  TH2D * avnoff;
  TH2D * spepvnratio;
  TH1D * rescorav;
  TH1D * rescorErr;
  TH1D * epvnErr;
  TH1D * epvnobsErr;
  TH1D * spvnErr;
  TH1D * avptErr;
  TH1D * avnoffErr;
  TH1D * vnobs[NumPtBins];
  TH2D * vnobs2D[NumPtBins];
  TH2D * qab;
  TH2D * qac;
  TH2D * qbc;
  TH2D * qna;

  int etacheck;
  int ptchecklo;
  int ptcheckhi;
  int count;
  int mintrk_;
  int maxtrk_;
  bool rescorCalculated;
  bool epvnobsCalculated;
  bool specRawCreated;
  bool specEffCorrectedCreated;
  TrackEfficiency * trkeff;
  GetEventInfo * info_;
  bool ispPb2011;
  bool isPbPb2015pixel;
};
void GenSupport::flip(){
  ep2a_ = flipENUM(ep2aSet_);
  ep2b_ = flipENUM(ep2bSet_);
  ep2c_ = flipENUM(ep2cSet_);
  ep3a_ = flipENUM(ep3aSet_);
  ep3b_ = flipENUM(ep3bSet_);
  ep3c_ = flipENUM(ep3cSet_);
}
GenSupport::GenSupport(int eporder, GetEventInfo * info, int ep2a, int ep2b, int ep2c,
int ep3a, int ep3b, int ep3c, int mintrk, int maxtrk, int subevt) {
  info_ = info;
  ispPb2011 = false;
  isPbPb2015pixel = false;
  if(info_->ispPb2011()) ispPb2011 = true;
  if(info_->isPbPb2015pixel()) isPbPb2015pixel = true;

  count = 0;
  epord_ = eporder;
  avcent_ = (maxtrk+mintrk)/2.;
  rescorCalculated = false;
  epvnobsCalculated = false;
  specRawCreated = false;
  specEffCorrectedCreated = false;
  if(ispPb2011) {
    TFile * feff = new TFile("EfficiencyTables/TrackCorrections_HIJING_538_OFFICIAL_Mar24.root");
    //cout<<"Efficiencies from: "<<feff->GetName()<<endl;
    heffa = (TH2D *) ((TH2D *) feff->Get("rTotalEff3D"))->Clone("heffa");
    heffa->SetDirectory(0);
    feff->Close();
  }
  if(isPbPb2015pixel) {
    TFile * feff = new TFile("EfficiencyTables/EffCorrectionsPixel_TT_pt_0_10_v2.root");
    //cout<<"Efficiencies from: "<<feff->GetName()<<endl;
    trkeff = new TrackEfficiency(feff);
  }
  mintrk_ = mintrk;
  maxtrk_ = maxtrk;
  ep2a_ = ep2a;
  ep2b_ = ep2b;
  ep2c_ = ep2c;
  ep2aSet_ = ep2a;
  ep2bSet_ = ep2b;
  ep2cSet_ = ep2c;
  ep3a_ = ep3a;
  ep3b_ = ep3b;
  ep3c_ = ep3c;
  ep3aSet_ = ep3a;
  ep3bSet_ = ep3b;
  ep3cSet_ = ep3c;

  AName2 = hi::EPNames[ep2a].data();
  BName2 = hi::EPNames[ep2b].data();
  CName2 = hi::EPNames[ep2c].data();
  AName3 = hi::EPNames[ep3a].data();
  BName3 = hi::EPNames[ep3b].data();
  CName3 = hi::EPNames[ep3c].data();
  TString utag = Form("_%s_%s_%s_%s_%s_%s_%s_%d_%d_%d",info->getFileName().Data(),
		      AName2.Data(),BName2.Data(),CName2.Data(),
		      AName3.Data(),BName3.Data(),CName3.Data(),
		      mintrk,maxtrk,subevt);
  string rnge = to_string(mintrk)+"-"+to_string(maxtrk)+"\%";
  ptbins = (TH1D*) info->getPtHist()->Clone("ptbins");
  double ptb[100];
  for(int i = 1; i<=ptbins->GetNbinsX()+1; i++) ptb[i-1]=ptbins->GetXaxis()->GetBinLowEdge(i);
  etabins = (TH1D*) info->getEtaHist()->Clone("etabins");
  double etab[100];
  for(int i = 1; i<=etabins->GetNbinsX()+1; i++) etab[i-1]=etabins->GetXaxis()->GetBinLowEdge(i);
  for(int i = 0; i<NumPtBins; i++) {
    vnobs[i] = new TH1D(Form("vnobs_%d_%s",i,utag.Data()),"vnobs",100,0,1);
    vnobs2D[i] = new TH2D(Form("vnobs2D_%d_%s",i,utag.Data()),"vnob2D",50,-1,1,50,-1,1);
    vnobs[i]->SetXTitle(Form("|q_{n}| (%4.1f < p_{T} < %4.1f GeV/c)",ptb[i],ptb[i+1]));
    vnobs2D[i]->SetOption("colz");
    vnobs2D[i]->SetXTitle(Form("q_{n,x} (%4.1f < p_{T} < %4.1f GeV/c)",ptb[i],ptb[i+1]) );
    vnobs2D[i]->SetYTitle("q_{n,y}");
  }
  rescor = new TH2D(Form("rescor%s",utag.Data()),"rescor",ptbins->GetNbinsX(),ptb,
		    etabins->GetNbinsX(),etab);
  rescor->Sumw2();
  rescor->SetOption("colz");
  rescor->Reset();
  rescor->SetYTitle("#eta");
  string pttitle = "p_{t} ("+rnge+")";
  string etatitle = "#eta_{t} ("+rnge+")";
  epvn             = (TH2D *) rescor->Clone(Form("epvn%s",utag.Data()));
  epvn->SetXTitle(pttitle.data());
  epvnobs          = (TH2D *) rescor->Clone(Form("epvnobs%s",utag.Data()));
  epvnobs->SetXTitle(pttitle.data());
  spvn             = (TH2D *) rescor->Clone(Form("spvn%s",utag.Data()));
  spvn->SetXTitle(pttitle.data());
  qab             = (TH2D *) rescor->Clone(Form("qab%s",utag.Data()));
  qab->SetXTitle(pttitle.data());
  qac             = (TH2D *) rescor->Clone(Form("qac%s",utag.Data()));
  qac->SetXTitle(pttitle.data());
  qbc             = (TH2D *) rescor->Clone(Form("qbc%s",utag.Data())); 
  qbc->SetXTitle(pttitle.data());
  qna             = (TH2D *) rescor->Clone(Form("qna%s",utag.Data()));
  qna->SetXTitle(pttitle.data());
  specRaw          = (TH2D *) rescor->Clone(Form("specRaw%s",utag.Data()));
  specRaw->SetXTitle(pttitle.data());
  specEffCorrected = (TH2D *) rescor->Clone(Form("specEffCorrected%s",utag.Data()));
  specEffCorrected->SetXTitle(pttitle.data());
  avpt_            = (TH2D *) rescor->Clone(Form("avpt%s",utag.Data()));
  avpt_->SetXTitle(pttitle.data());
  avnoff           = (TH2D *) rescor->Clone(Form("avnoff%s",utag.Data()));
  avnoff->SetXTitle(pttitle.data());
  nfrac           = (TH2D *) rescor->Clone(Form("nfrac%s",utag.Data()));
  nfrac->SetXTitle(pttitle.data());
  epInt = new TH1D(Form("epInt%s",utag.Data()),"rescor",etabins->GetNbinsX(),etabins->GetXaxis()->GetXbins()->GetArray());
  epInt->SetXTitle(etatitle.data());
  spInt = new TH1D(Form("spInt%s",utag.Data()),"rescor",etabins->GetNbinsX(),etabins->GetXaxis()->GetXbins()->GetArray());
  spInt->SetXTitle(etatitle.data());
  spepintratio = new TH1D(Form("spepintratio%s",utag.Data()),"rescor",etabins->GetNbinsX(),etabins->GetXaxis()->GetXbins()->GetArray());
  spepintratio->SetXTitle(etatitle.data());
  q.resize(etabins->GetNbinsX(), vector<Qvec *>(ptbins->GetNbinsX()));
  qav = new Qvec(0,0,"qav");
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      string name = Form("q_%d_%d",i,j);
      q.at(j).at(i) = new Qvec(i,j,name);
    }
  }

  for(int k = 0; k<ERRCNT; k++) {
    string avname = "qaverr_"+to_string(k);
    qaverr[k] = new Qvec(0,0,avname);
    qerr[k].resize(etabins->GetNbinsX(), vector<Qvec *>(ptbins->GetNbinsX()));
    for(int i = 0; i<ptbins->GetNbinsX(); i++) {
      for(int j = 0; j<etabins->GetNbinsX(); j++) {
	string name = Form("qerr_%d_%d_%d",k,i,j);
	qerr[k].at(j).at(i) = new Qvec(i,j,name);
      }
    }
  }

  rescorav   = new TH1D(Form("rescorav%s",utag.Data()),"rescorav",4,1,4); 
  rescorav->Sumw2(); 
  rescorErr  = new TH1D(Form("%sErr",rescor->GetName()),"rescorErr",1000,0,1);
  epvnErr    = new TH1D(Form("%sErr",epvn->GetName()),"epvnErr",2000,-0.4,0.6);
  epvnobsErr = new TH1D(Form("%sErr",epvnobs->GetName()),"epvnobsErr",2000,-0.1,0.5);
  spvnErr    = new TH1D(Form("%sErr",spvn->GetName()), "spvnErr",2000,-0.4,0.6);
  avptErr    = new TH1D(Form("%sErr",avpt_->GetName()), "avptErr",2000,0,10.);
  avnoffErr  = new TH1D(Form("%sErr",avnoff->GetName()), "avnoffErr",2000,0,300);
  etacheck = 5;
  ptchecklo = ptbins->FindBin(1.0);
  ptcheckhi = ptbins->FindBin(7.0);
}

void GenSupport::add(TH2D * qxtrk, TH2D * qytrk, TH2D * qcnt, double * qx, double * qy, double * qw, TH2D *  avpt, int evtchar, TH2D * qxav, TH2D * qyav){
  ++count;
  ep2a_ = ep2aSet_;
  ep2b_ = ep2bSet_;
  ep2c_ = ep2cSet_;
  ep3a_ = ep3aSet_;
  ep3b_ = ep3bSet_;
  ep3c_ = ep3cSet_;

  if(ispPb ) flip();
  comp q2A(qx[ep2a_],qy[ep2a_]);
  comp q2B(qx[ep2b_],qy[ep2b_]);
  comp q2C(qx[ep2c_],qy[ep2c_]);
  comp q3A(qx[ep3a_],qy[ep3a_]);
  comp q3B(qx[ep3b_],qy[ep3b_]);
  comp q3C(qx[ep3c_],qy[ep3c_]);
  comp n2A;
  comp n3A; 
  comp Qvec (0.,0.);
  comp n2Anorm;
  comp n3Anorm; 
  comp qDummy(1.,0.);
  int iran = (int) ran->Uniform(0,ERRCNT-0.0001);
  qav->add(1,qDummy,q2A,qw[ep2a_],q2B,qw[ep2b_],q2C,qw[ep2c_],
  	   q3A,qw[ep3a_],q3B,qw[ep3b_],q3C,qw[ep3c_],1,1);
  qaverr[iran]->add(1,qDummy,q2A,qw[ep2a_],q2B,qw[ep2b_],q2C,qw[ep2c_],
  		    q3A,qw[ep3a_],q3B,qw[ep3b_],q3C,qw[ep3c_],1,1);

  int jeta=0;
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    double ptcnt = 0;
    Qvec = zero;
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      jeta = j+1;
      double eta = etabins->GetXaxis()->GetBinCenter(jeta);
      if(ispPb) jeta = etabins->GetNbinsX()-j;
      comp qn(qxtrk->GetBinContent(i+1,jeta)-qxav->GetBinContent(i+1,jeta), qytrk->GetBinContent(i+1,jeta)-qyav->GetBinContent(i+1,jeta));
      int mult = qcnt->GetBinContent(i+1,jeta);
      double pt= avpt->GetBinContent(i+1,jeta);
      if(eta>-0.8&&eta<0.8&&mult>1){
      	Qvec+=qn;
      	ptcnt+=mult;
      }
      int iran = (int) ran->Uniform(0,ERRCNT-0.0001);
      q.at(j).at(i)->add(mult,qn,q2A,qw[ep2a_], q2B,qw[ep2b_], q2C,qw[ep2c_],
      			 q3A,qw[ep3a_], q3B,qw[ep3b_], q3C,qw[ep3c_],pt,evtchar);
      qerr[iran].at(j).at(i)->add(mult,qn,q2A,qw[ep2a_],q2B,qw[ep2b_],q2C,qw[ep2c_],
      				  q3A,qw[ep3a_], q3B,qw[ep3b_], q3C,qw[ep3c_],pt,evtchar);
    }
    if(ptcnt>1){
      Qvec/=ptcnt;
      if(i<NumPtBins){
    	vnobs[i]->Fill(std::abs(Qvec));
    	vnobs2D[i]->Fill(Qvec.real(),Qvec.imag());
      }
    }
  }
}

TH2D * GenSupport::getQnFraction(){
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      nfrac->SetBinContent(i+1,j+1, q.at(j).at(i)->getnfrac());
      double err = 0;
      for(int k = 0; k<ERRCNT; k++) {
	err+=pow(qerr[k].at(j).at(i)->getnfrac()-q.at(j).at(i)->getnfrac(),2);
      }
      err=sqrt(err/(ERRCNT-1.));
      err/=sqrt((double)ERRCNT);
      if(err<1e-5) err=0;
      nfrac->SetBinError(i+1,j+1,err);
    }
  }
  return nfrac;
}

TH2D * GenSupport::getRescor(){
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      rescor->SetBinContent(i+1,j+1, q.at(j).at(i)->getRescor());
      double err = 0;
      double ecnt = 0;
      if(ERRCNT == 1) {
      	err = q.at(j).at(i)->getRescorErr();
      	ecnt = 1.;
      } else {
      	for(int k = 0; k<ERRCNT; k++) {
      	  if(qerr[k].at(j).at(i)==0) continue;
      	  err+=pow(qerr[k].at(j).at(i)->getRescor()-q.at(j).at(i)->getRescor(),2);
      	  ecnt+=1.;
      	}
      	if(ecnt>5) {
      	  err=sqrt(err/(ecnt-1.));
      	  err/=sqrt((double)ecnt);
      	}
      }
      if(err<1e-5) err=0;
      rescor->SetBinError(i+1,j+1,err);
    }
  }
  
  rescorCalculated = true;
  return rescor;
}

TH1D * GenSupport::getRescorAv(){
  if(std::isnormal(qav->getRescor())) {
    rescorav->SetBinContent(1, qav->getRescor());
    rescorav->SetBinContent(3, qav->getRescor());
  } else {
    rescorav->SetBinContent(1, 0);
    rescorav->SetBinContent(3, 0);
    rescorav->SetBinError(1, 0);
    rescorav->SetBinError(3, 0);
  }
  if(std::isnormal(qav->getSPVNdenom())) {
    rescorav->SetBinContent(2, qav->getSPVNdenom());
    rescorav->SetBinContent(4, qav->getSPVNdenom());
  } else {
    rescorav->SetBinContent(2, 0);
    rescorav->SetBinContent(4, 0);
    rescorav->SetBinError(2, 0);
    rescorav->SetBinError(4, 0);
  }
  double err = 0;
  double errsp = 0;

  if(std::isnormal(qav->getRescorErr())) {
    rescorav->SetBinError(3, qav->getRescorErr());
  } else {
    rescorav->SetBinError(3, 0);
  }
  if(std::isnormal(qav->getSPVNdenomerr())) {
    rescorav->SetBinError(4, qav->getSPVNdenomerr());
  } else {
    rescorav->SetBinError(4, 0.);
  }
  if(ERRCNT == 1) {
    err = qav->getRescorErr();
    errsp = qav->getSPVNdenomerr();
  } else {
    double npt = 0;
    double nptsp = 0;
    for(int k = 0; k<ERRCNT; k++) {
      double rn = qaverr[k]->getRescor();
      double rav = qav->getRescor();
      if(std::isnormal(rn) && std::isnormal(rav)) {
	++npt;
	err+=pow(qaverr[k]->getRescor()-qav->getRescor(),2);
      }

      double rnsp = qaverr[k]->getSPVNdenom();
      double ravsp = qav->getSPVNdenom();
      if(std::isnormal(rnsp) && std::isnormal(ravsp)) {
	++nptsp;
	errsp+=pow(qaverr[k]->getSPVNdenom()-qav->getSPVNdenom(),2);
      }
    }
    err=sqrt(err/(npt-1.));
    err/=sqrt(npt);
    errsp=sqrt(errsp/(nptsp-1.));
    errsp/=sqrt(nptsp);
  }
  if(err<1e-6) err=0;
  if(errsp<1e-6) errsp=0;
  rescorav->SetBinError(1,err);
  rescorav->SetBinError(2,errsp);
  return rescorav;
}

TH2D * GenSupport::getEPVNobs(){
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      epvnobs->SetBinContent(i+1,j+1, q.at(j).at(i)->getEPVNobs());
      double err = 0;
      if(ERRCNT == 1) {
	err = q.at(j).at(i)->getEPVNobsErr();
	if(err<0.01) err=0;
      } else {
	for(int k = 0; k<ERRCNT; k++) {
	  err+=pow(qerr[k].at(j).at(i)->getEPVNobs()-q.at(j).at(i)->getEPVNobs(),2);
	}
	err=sqrt(err/(ERRCNT-1.));
	err/=sqrt((double)ERRCNT);
      }
      epvnobs->SetBinError(i+1,j+1,err);
    }
  }
  epvnobsCalculated = true;
  return epvnobs;
}


TH2D * GenSupport::getSPEPRatio(){
  spepvnratio = (TH2D *) getSPVN()->Clone("spepvnratio");
  spepvnratio->Divide(getEPVN());
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      double sp =  q.at(j).at(i)->getSPVN();
      double ep = q.at(j).at(i)->getEPVN();
      double err = 0;
      if(ERRCNT==1) {
	err = q.at(j).at(i)->getSPVNErr();
      } else {
	for(int k = 0; k<ERRCNT; k++) {
	  double spsub = qerr[k].at(j).at(i)->getSPVN();
	  double epsub = qerr[k].at(j).at(i)->getEPVN();

	  err+=pow(spsub/epsub-sp/ep,2);
	}
	err=sqrt(err/(ERRCNT-1.));
	err/=sqrt((double)ERRCNT);
      }
      spepvnratio->SetBinError(i+1,j+1,err);
    }
  }

  return spepvnratio;
}



TH2D * GenSupport::getSPVN(){
  TH1D * res = getRescorAv();
  double denom = res->GetBinContent(2);
  double denerr = res->GetBinError(2);
  if(!std::isnormal(denom)) {
    spvn->Reset();
    return spvn;
  }
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      if(!std::isnormal(q.at(j).at(i)->getSPVNnum())) {
	spvn->SetBinContent(i+1,j+1,0);
	spvn->SetBinError(i+1,j+1,0);
	qab->SetBinContent(i+1,j+1,0);
	qab->SetBinError(i+1,j+1,0);
	qac->SetBinContent(i+1,j+1,0);
	qac->SetBinError(i+1,j+1,0);
	qbc->SetBinContent(i+1,j+1,0);
	qbc->SetBinError(i+1,j+1,0);
	qna->SetBinContent(i+1,j+1,0);
	qna->SetBinError(i+1,j+1,0);
	continue;
      };
      spvn->SetBinContent(i+1,j+1, q.at(j).at(i)->getSPVNnum()/denom);
      qab->SetBinContent(i+1,j+1,q.at(j).at(i)->getQAB());
      qac->SetBinContent(i+1,j+1,q.at(j).at(i)->getQAC());
      qbc->SetBinContent(i+1,j+1,q.at(j).at(i)->getQBC());
      qna->SetBinContent(i+1,j+1,q.at(j).at(i)->getQnA());

      double err = 0;
      if(ERRCNT==1) {
	if(std::isnormal(q.at(j).at(i)->getSPVNnumerr())) err = q.at(j).at(i)->getSPVNnumerr();
      } else {
	double ncnt = 0;
	for(int k = 0; k<ERRCNT; k++) {
	  double spn = qerr[k].at(j).at(i)->getSPVNnum();
	  double spav = q.at(j).at(i)->getSPVNnum();
	  if(std::isnormal(spn) && std::isnormal(spav)) {
	    err+=pow(spn-spav,2);
	    ++ncnt;
	  }
	}
	if(ncnt > 1) {
	  err=sqrt(err/(ncnt-1.));
	  err/=sqrt(ncnt);
	}
      }
      double qret = q.at(j).at(i)->getSPVNnum();


      if(std::isnormal(denom) && std::isnormal(err) && std::isnormal(qret) && std::isnormal(denerr)) {
	double arg = pow(err/qret,2.)+pow(denerr/denom,2);
	if(std::isnormal(arg)) err=spvn->GetBinContent(i+1,j+1)*sqrt( arg);
      } else if(std::isnormal(err) && std::isnormal(qret)) {
	double arg = pow(err/qret,2.);
	if(std::isnormal(arg))err=spvn->GetBinContent(i+1,j+1)*sqrt( arg );
      } else {
	err = 0;
      }
      if(std::isnormal(err)) {
	spvn->SetBinError(i+1,j+1,err);
      } else {
	spvn->SetBinError(i+1,j+1,0.);
      }
    }
  }
  return spvn;
}


TH2D * GenSupport::getEPVN(){
  TH1D * res = getRescorAv();
  double resc = res->GetBinContent(1);
  double reserr = res->GetBinError(1);
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      epvn->SetBinContent(i+1,j+1, q.at(j).at(i)->getEPVNobs()/resc);
      double err = 0;
      if(ERRCNT==1) {
	err = q.at(j).at(i)->getEPVNobsErr()/resc;
      } else {
	for(int k = 0; k<ERRCNT; k++) {
	  err+=pow(qerr[k].at(j).at(i)->getEPVNobs()-q.at(j).at(i)->getEPVNobs(),2);
	}
	err=sqrt(err/(ERRCNT-1.0));
	err/=sqrt((double)ERRCNT);
      }
      err=epvn->GetBinContent(i+1,j+1)*sqrt( pow(err/q.at(j).at(i)->getEPVNobs(),2.)+pow(reserr/resc,2.) );
      epvn->SetBinError(i+1,j+1,err);
    }
  }
  return epvn;
}



TH2D * GenSupport::getAvPt(){
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      avpt_->SetBinContent(i+1,j+1, q.at(j).at(i)->getAvPt());
      double err = 0;
      if(ERRCNT==1) {
	err = q.at(j).at(i)->getAvPtErr();
      } else {
	for(int k = 0; k<ERRCNT; k++) {
	  err+=pow(qerr[k].at(j).at(i)->getAvPt()-q.at(j).at(i)->getAvPt(),2);
	}
	err=sqrt(err/(ERRCNT-1.));
	err/=sqrt((double)ERRCNT);
      }
      avpt_->SetBinError(i+1,j+1,err);
    }
  }
  return avpt_;
}


TH2D * GenSupport::getAvNoff(){
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      avnoff->SetBinContent(i+1,j+1, q.at(j).at(i)->getAvNoff());
      double err = 0;
      if(ERRCNT==1) {
	err = q.at(j).at(i)->getAvNoffErr();
      } else {
	for(int k = 0; k<ERRCNT; k++) {
	  err+=pow(qerr[k].at(j).at(i)->getAvNoff()-q.at(j).at(i)->getAvNoff(),2);
	}
	err=sqrt(err/(ERRCNT-1.));
	err/=sqrt((double)ERRCNT);
      }
      avnoff->SetBinError(i+1,j+1,err);
    }
  }
  return avnoff;
}

TH2D * GenSupport::getSpecRaw(){
  
  for(int i = 0; i<ptbins->GetNbinsX(); i++) {
    for(int j = 0; j<etabins->GetNbinsX(); j++) {
      specRaw->SetBinContent(i+1,j+1, q.at(j).at(i)->getCount());
      specRaw->SetBinError(i+1,j+1, sqrt(q.at(j).at(i)->getCount()));
    }
  }
  specRawCreated = true;
  return specRaw;
}
TH2D * GenSupport::getSpecEffCorrected(){
  if(!specRawCreated) {getSpecRaw();}
  for(int i = 1; i<=specRaw->GetNbinsX(); i++) {
    for(int j = 1; j<=specRaw->GetNbinsY(); j++) {
      double pt = specRaw->GetXaxis()->GetBinCenter(i);
      double eta = specRaw->GetYaxis()->GetBinCenter(j);
      double cent = avcent_;
      double eff = 0 ;
      if(ispPb2011) {
	eff= heffa->GetBinContent( heffa->GetXaxis()->FindBin(eta), heffa->GetYaxis()->FindBin(pt));
      }
      if(isPbPb2015pixel) {
	eff = trkeff->GetEff(cent,pt,eta);
      }
      if(eff<1e-4) eff = 1;
      double correctedCount = specRaw->GetBinContent(i,j)/eff;
      double correctedCountErr = specRaw->GetBinError(i,j)/eff;
      specEffCorrected->SetBinContent(i,j,correctedCount);
      specEffCorrected->SetBinError(i,j,correctedCountErr);
    }
  }
  specEffCorrectedCreated = true;
  return specEffCorrected;
}
TH1D * GenSupport::getEPInt(double minpt, double maxpt){
  if(!specEffCorrectedCreated) getSpecEffCorrected();
  int imin = specRaw->GetXaxis()->FindBin( minpt );
  int imax = specRaw->GetXaxis()->FindBin( maxpt - 0.0001);

  for(int j = 1; j<= epInt->GetNbinsX(); j++) {
    double sum = 0;
    double sumw = 0;
    double sumerr = 0;

    for(int i =imin; i<=imax; i++) {
      double vn = epvn->GetBinContent(i,j);
      double vnerr = epvn->GetBinError(i,j);
      double w =specEffCorrected->GetBinContent(i,j);
      double werr = specEffCorrected->GetBinError(i,j);
      sum+=w*vn;
      sumw+=w;
      sumerr+=pow(w*vn,2)*(pow(vnerr/vn,2)+pow(werr/w,2));
    }
    epInt->SetBinContent(j, sum/sumw);
    epInt->SetBinError(j,sqrt(sumerr)/sumw);
  }
  return epInt;
};
TH1D * GenSupport::getSPInt(double minpt, double maxpt){
  TH1D * res = getRescorAv();
  double denom = res->GetBinContent(2);
  if(!std::isnormal(denom)) {
    spInt->Reset();
    return spInt;
  }
  if(!specEffCorrectedCreated) getSpecEffCorrected();
  int imin = specRaw->GetXaxis()->FindBin( minpt );
  int imax = specRaw->GetXaxis()->FindBin( maxpt - 0.0001);

  for(int j = 1; j<= epInt->GetNbinsX(); j++) {
    double sum = 0;
    double sumw = 0;
    double sumerr = 0;
    for(int i =imin; i<=imax; i++) {
      double vn = spvn->GetBinContent(i,j);
      double vnerr = spvn->GetBinError(i,j);
      double w =specEffCorrected->GetBinContent(i,j);
      double werr = specEffCorrected->GetBinError(i,j);
      if(std::isnormal(vn)&&std::isnormal(vnerr)&&std::isnormal(w)&&std::isnormal(werr)) {
	sum+=w*vn;
	sumw+=w;
	sumerr+=pow(w*vn,2)*(pow(vnerr/vn,2)+pow(werr/w,2));
      }
    }
    spInt->SetBinContent(j, sum/sumw);
    spInt->SetBinError(j,sqrt(sumerr)/sumw);
  }
  return spInt;
};
TH1D * GenSupport::getSPEPIntRatio(double minpt, double maxpt){
  TH1D * hres = getRescorAv();
  double res = hres->GetBinContent(1);
  double den = hres->GetBinContent(2);

  int imin = specRaw->GetXaxis()->FindBin( minpt );
  int imax = specRaw->GetXaxis()->FindBin( maxpt - 0.0001);
  TH1D * spint = getSPInt(minpt, maxpt);
  TH1D * epint = getEPInt(minpt, maxpt);
  for(int j = 1; j<=spint->GetNbinsX(); j++) {

    double r = spint->GetBinContent(j)/epint->GetBinContent(j);
    spepintratio->SetBinContent(j,r);
    double err = 0;
    if(ERRCNT==1) {
      double sp = spint->GetBinContent(j);
      err = spint->GetBinError(j);
      err = r*err/sp;
    } else {
      double sumr = 0;
      double sumw = 0;

      for(int i = imin; i<=imax; i++) {
	double sp = q.at(j-1).at(i-1)->getSPVNnum()/den;
	double ep = q.at(j-1).at(i-1)->getEPVNobs()/res;
	r = sp/ep;
	double w = specEffCorrected->GetBinContent(i,j);
	sumr+=w*r;
	sumw+=w;
      }
      double avr = sumr/sumw;
      
      for(int k = 0; k<ERRCNT; k++) {	
	sumr = 0;
	sumw = 0;
	for(int i = imin; i<=imax; i++) {
	  double sp = qerr[k].at(j-1).at(i-1)->getSPVNnum()/den;
	  double ep = qerr[k].at(j-1).at(i-1)->getEPVNobs()/res;
	  r = sp/ep;
	  double w = specEffCorrected->GetBinContent(i,j);
	  sumr+=w*r;
	  sumw+=w;
	}
	double avrk = sumr/sumw;
	err+=pow(avrk - avr,2);
	}
	err=sqrt(err/(ERRCNT-1.));
	err/=sqrt((double)ERRCNT);
    }
    spepintratio->SetBinError(j,err);
  }
  return spepintratio;
};

#endif
