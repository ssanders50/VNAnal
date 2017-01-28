#include "MCEvent.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TDirectory.h"
#include "HiEvtPlaneList.h"
#include "TRandom2.h"
using namespace hi;
static const int npt = 18;
static const double ptbins[]={0.3,0.4,0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
				       2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 
				       10.0, 12.};

static const int netabinsDefault = 12;
static const double etabinsDefault[]={-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};

Double_t bounds(int ord, double ang) {
  while(ang>TMath::Pi()/ord)  ang-=TMath::TwoPi()/ord;
  while(ang< -TMath::Pi()/ord) ang+=TMath::TwoPi()/ord;
  return ang;
}
void RunCase(Int_t nevents=20, double dndeta=100, Double_t setv1=0.01, Double_t setv2=0.06,
	     Double_t setv3 = 0., Double_t setv4 = 0., Double_t setv5 = 0., Double_t setv6 = 0., Double_t setv7 = 0.);

void toyMC() {
  RunCase(5000000, 600, 0.0, 0.1, 0.03, 0.02, 0.01, 0.01,0.01);
  //RunCase(50000,400,0.005,0.06,kFALSE);
}

void RunCase(Int_t nevents, double  dndeta, Double_t setv1, Double_t setv2, Double_t setv3, Double_t setv4, Double_t setv5,
	     Double_t setv6, Double_t setv7) {
  string fname = "/rfs/sanders/toyMC/MC_vn_"+to_string((int)(100*setv2))+"_"
    +to_string((int)(100*setv3))+"_"+to_string((int)(100*setv4))+
    "_"+to_string((int)(100*setv5))+"_"+to_string((int)(100*setv6))+"_"+to_string((int)(100*setv7))+"_E.root";
  TFile * tf = new TFile(fname.data(),"recreate");
  TDirectory * dir = tf->mkdir("MC");
  dir->cd();
  double centval;
  int Noff;
  double vtx;
  Double_t epang[NumEPNames];
  Double_t eporig[NumEPNames];
  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t sumw[NumEPNames];
  Double_t qcnt[NumEPNames];
  Double_t vn[NumEPNames];
  Double_t epmult[NumEPNames];
  unsigned int runno_;
  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  TH2D * qxtrk2;
  TH2D * qytrk2;
  TH2D * qxtrk3;
  TH2D * qytrk3;
  TH2D * qxtrk4;
  TH2D * qytrk4;
  TH2D * qxtrk5;
  TH2D * qytrk5;
  TH2D * qxtrk6;
  TH2D * qytrk6;
  TH2D * qxtrk7;
  TH2D * qytrk7;

  TH2D * qxycnt;
  TH2D * avpt;
  TH1D * hpt = new TH1D("hpt","hpt",1000,0,12);
  TH1D * hphi = new TH1D("hphi","hphi",1000,-4,4);
  TH1D * hv2 = new TH1D("hv2","hv2",1000,-1.1,1.1);
  TH1D * hv3 = new TH1D("hv3","hv3",1000,-1.1,1.1);
  TH1D * hv4 = new TH1D("hv4","hv4",1000,-1.1,1.1);
  TH1D * hv5 = new TH1D("hv5","hv5",1000,-1.1,1.1);
  TH1D * hv6 = new TH1D("hv6","hv6",1000,-1.1,1.1);
  TH1D * hv7 = new TH1D("hv7","hv7",1000,-1.1,1.1);
  qxtrk2 = new TH2D("qxtrk_v2","qxtrk_v2",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk2 = new TH2D("qytrk_v2","qytrk_v2",npt,ptbins, netabinsDefault, etabinsDefault);
  qxycnt = new TH2D("qxycnt"  ,"qxycnt"  ,npt,ptbins, netabinsDefault, etabinsDefault);
  avpt =  new TH2D("avpt", "avpt",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk2->Sumw2();
  qxtrk2->SetOption("colz");
  qytrk2->Sumw2();
  qytrk2->SetOption("colz");
  qxycnt->Sumw2();
  qxycnt->SetOption("colz");
  qxtrk3 = new TH2D("qxtrk_v3","qxtrk_v3",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk3 = new TH2D("qytrk_v3","qytrk_v3",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk3->Sumw2();
  qxtrk3->SetOption("colz");
  qytrk3->Sumw2();
  qytrk3->SetOption("colz");
  qxtrk4 = new TH2D("qxtrk_v4","qxtrk_v4",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk4 = new TH2D("qytrk_v4","qytrk_v4",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk4->Sumw2();
  qxtrk4->SetOption("colz");
  qytrk4->Sumw2();
  qytrk4->SetOption("colz");
  qxtrk5 = new TH2D("qxtrk_v5","qxtrk_v5",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk5 = new TH2D("qytrk_v5","qytrk_v5",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk5->Sumw2();
  qxtrk5->SetOption("colz");
  qytrk5->Sumw2();
  qytrk5->SetOption("colz");
  qxtrk6 = new TH2D("qxtrk_v6","qxtrk_v6",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk6 = new TH2D("qytrk_v6","qytrk_v6",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk6->Sumw2();
  qxtrk6->SetOption("colz");
  qytrk6->Sumw2();
  qytrk6->SetOption("colz");
  qxtrk7 = new TH2D("qxtrk_v7","qxtrk_v7",npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk7 = new TH2D("qytrk_v7","qytrk_v7",npt,ptbins, netabinsDefault, etabinsDefault);
  qxtrk7->Sumw2();
  qxtrk7->SetOption("colz");
  qytrk7->Sumw2();
  qytrk7->SetOption("colz");

  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";

  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
  }
  MCEvent * event = new MCEvent(setv1,setv2,setv3,setv4,setv5,setv6,setv7);
  TRandom2 * ran = new TRandom2(0);
  TTree * tree = new TTree("tree","EP tree");
  tree->Branch("Cent",&centval,"cent/D");
  tree->Branch("NtrkOff",&Noff,"Noff/I");
  tree->Branch("Vtx",&vtx,"vtx/D");
  tree->Branch("epang",&epang, epnames.Data());
  tree->Branch("eporig",&eporig, epnames.Data());
  tree->Branch("sumw",&sumw, epnames.Data());
  tree->Branch("qx",      &qx,       epnames.Data());
  tree->Branch("qy",      &qy,       epnames.Data());
  tree->Branch("q",       &q,       epnames.Data());
  tree->Branch("vn", &vn, epnames.Data());
  tree->Branch("mult",    &epmult,  epnames.Data());
  tree->Branch("Run",     &runno_,   "run/i");
  tree->Branch("Rescor",  &rescor,   epnames.Data());
  tree->Branch("RescorErr",  &rescorErr,   epnames.Data());
  tree->Branch("qxtrk_v2",   "TH2D",  &qxtrk2, 128000, 0);
  tree->Branch("qytrk_v2",   "TH2D",  &qytrk2, 128000, 0);
  tree->Branch("qxtrk_v3", "TH2D",  &qxtrk3, 128000, 0);
  tree->Branch("qytrk_v3", "TH2D",  &qytrk3, 128000, 0);
  tree->Branch("qxtrk_v4", "TH2D",  &qxtrk4, 128000, 0);
  tree->Branch("qytrk_v4", "TH2D",  &qytrk4, 128000, 0);
  tree->Branch("qxtrk_v5", "TH2D",  &qxtrk5, 128000, 0);
  tree->Branch("qytrk_v5", "TH2D",  &qytrk5, 128000, 0);
  tree->Branch("qxtrk_v6", "TH2D",  &qxtrk6, 128000, 0);
  tree->Branch("qytrk_v6", "TH2D",  &qytrk6, 128000, 0);
  tree->Branch("qxtrk_v7", "TH2D",  &qxtrk7, 128000, 0);
  tree->Branch("qytrk_v7", "TH2D",  &qytrk7, 128000, 0);
  tree->Branch("qcnt",  "TH2D",  &qxycnt, 128000, 0);
  tree->Branch("avpt",    "TH2D",  &avpt, 128000, 0);
  for(int ievent = 0; ievent < nevents; ievent++){
    if(fmod( ievent, nevents/20)==0) cout<<(double)ievent/(double)nevents<<endl;
    qxtrk2->Reset();
    qytrk2->Reset();
    qxycnt->Reset();
    qxtrk3->Reset();
    qytrk3->Reset();
    qxtrk4->Reset();
    qytrk4->Reset();
    qxtrk5->Reset();
    qytrk5->Reset();
    qxtrk6->Reset();
    qytrk6->Reset();
    qxtrk7->Reset();
    qytrk7->Reset();
    avpt->Reset();
    for(int i = 0; i< NumEPNames; i++) {
      qx[i] = 0;
      qy[i] = 0;
      q[i] = 0;
      vn[i] = 0;
      qcnt[i] = 0;
      epmult[i] = 0;
      rescor[i] = 0;
      sumw[i] = 0;
      epang[i] = -10;
      eporig[i] = -10;
      
      
    }
    event->SetPsiRandom();
    double Psi = event->GetPsi();
    double phiarr[20000];
    double ptarr[20000];
    int mult = 10.*dndeta;
    event->SetMult(mult);
    event->GetThrowPhi(phiarr);
    event->GetPtRandom(ptarr);
    centval = ran->Uniform(0,70);
    vtx = 0;
    for(int imult = 0; imult<mult; imult++) {
      double eta = ran->Uniform(-5,5);
      double phi = phiarr[imult];
      hphi->Fill(phi);
      hv2->Fill(cos(2*phi));
      hv3->Fill(cos(3*phi));
      hv4->Fill(cos(4*phi));
      hv5->Fill(cos(5*phi));
      hv6->Fill(cos(6*phi));
      hv7->Fill(cos(7*phi));
      phi+=Psi;
      double pt = ptarr[imult];
      hpt->Fill(pt);
      qxtrk2->Fill(pt,eta, TMath::Cos(2.*phi));
      qytrk2->Fill(pt,eta, TMath::Sin(2.*phi));
      qxycnt->Fill(pt,eta);
      avpt->Fill(pt,eta,pt);
      qxtrk3->Fill(pt,eta, TMath::Cos(3.*phi));
      qytrk3->Fill(pt,eta, TMath::Sin(3.*phi));
      qxtrk4->Fill(pt,eta, TMath::Cos(4.*phi));
      qytrk4->Fill(pt,eta, TMath::Sin(4.*phi));
      qxtrk5->Fill(pt,eta, TMath::Cos(5.*phi));
      qytrk5->Fill(pt,eta, TMath::Sin(5.*phi));
      qxtrk6->Fill(pt,eta, TMath::Cos(6.*phi));
      qytrk6->Fill(pt,eta, TMath::Sin(6.*phi));
      qxtrk7->Fill(pt,eta, TMath::Cos(7.*phi));
      qytrk7->Fill(pt,eta, TMath::Sin(7.*phi));
      for(int epindx = 0; epindx<NumEPNames; epindx++){
	if(eta>=EPEtaMin1[epindx] && eta<EPEtaMax1[epindx]) {
	  qx[epindx] = TMath::Cos(EPOrder[epindx]*phi);
	  qy[epindx] = TMath::Sin(EPOrder[epindx]*phi);
	  ++epmult[epindx];
	  ++sumw[epindx];
	  ++qcnt[epindx];
	}
      }

    }
    for(int epindx = 0; epindx<NumEPNames; epindx++){
      if(qcnt[epindx]>5) {
	epang[epindx] = bounds(EPOrder[epindx],atan2(qy[epindx],qx[epindx])/EPOrder[epindx]);
	eporig[epindx] = epang[epindx];
      }
    }
    tree->Fill();
    //hpt->Draw();
  }
  hphi->Write();
  hv2->Write();
  hv3->Write();
  hv4->Write();
  hv5->Write();
  hv6->Write();
  hv7->Write();
  qxtrk2->Write("qxtrk_v2");
  qytrk2->Write("qytrk_v2");
  qxycnt->Write("qcnt_v2");
  //tf->Close();
}
