#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include <iostream>
#include <unistd.h>
#include <complex>
#include <cmath>

TRandom3 * ran;

typedef complex<double> comp;

int epord_ = 2.;
bool trkoff = true;
bool RECENTERTRACKS = false;

static const int MaxFiles = 1000;
static const int MaxEvents = -1;
static const int ntrkbins = 15;
static const double trkBins[]={0,20,30,40,50,60,80,100,120,150,185,220,260,300,350,500};
static const int ncentbins = 12;
static const double centBins[]={0,5,10,15,20,25,30,35,40,45,50,60,70};
static const int nanals = 28;
enum AnalType {
  N2SUB2,   N2SUB3,    N3SUB2,   N3SUB3,   N4SUB2,  N4SUB3,  
  N42SUB2,  N42SUB3,   N5SUB2,   N5SUB3,   N6SUB2,  N6SUB3,   
  N7SUB2,   N7SUB3,    N523SUB2, N523SUB3, N723SUB2,N723SUB3, 
  N63SUB2,  N63SUB3,   D24SUB2,  D24SUB3,  D34SUB2, D34SUB3,   
  D2232SUB2,D2232SUB3, D2432SUB2, D2432SUB3
};
string AnalNames[]={
  "N2SUB2",   "N2SUB3",    "N3SUB2",   "N3SUB3",   "N4SUB2",  "N4SUB3",  
  "N42SUB2",  "N42SUB3",   "N5SUB2",   "N5SUB3",   "N6SUB2",  "N6SUB3",   
  "N7SUB2",   "N7SUB3",    "N523SUB2", "N523SUB3", "N723SUB2","N723SUB3", 
  "N63SUB2",  "N63SUB3",   "D24SUB2",  "D24SUB3",  "D34SUB2", "D34SUB3",   
  "D2232SUB2","D2232SUB3", "D2432SUB2","D2432SUB3"
};
int ANAL;
int epa;
int epb;
int epc;
int ep3a;
int ep3b;
int ep3c;

#include "HiEvtPlaneList.h"
using namespace hi;
Bool_t ispPb ;

int GetMidIndx(Int_t epord, TString midn) {
  int mid = -1;
  string match = midn.Data();
  for(int i = 0; i<hi::NumEPNames; i++) {
    string ct = hi::EPNames[i];
    if(epord==3) ct=hi::EPNames[i];
    if(match==ct) mid = i;
  }
  if(mid<0) {cout<<midn.Data()<< "not found"<<endl; return -1;}
  return mid;
}


Int_t flipENUM(int ein){
  if(ein<0) {
    cout<<"flipENUM called with ein = "<<ein<<"  EXPECT TO CRASH!"<<endl;
  }
  if(!ispPb) return ein;
  string ename = hi::EPNames[ein];
  if(ename.find("mid") != std::string::npos) return ein;
  if(ename.find("m") != std::string::npos) {
    ename.replace(ename.find("m"),1,"p");
  } else if (ename.find("p") != std::string::npos){
    ename.replace(ename.find("p"),1,"m");
  }
  //int epord = 2;
  //if(ename.find("3")!=std::string::npos) epord = 3;
  return GetMidIndx(epord_,ename);
}

#include "src/GetEventInfo.h"

string rpnames[hi::NumEPNames];
//----------------------------------
// Tree Variables:
//
double centval;
int noff;
double vtx;
double epang[hi::NumEPNames];
Double_t qx[hi::NumEPNames];
Double_t qy[hi::NumEPNames];
Double_t q[hi::NumEPNames];
Double_t epmult[hi::NumEPNames];
unsigned int  runno_;
Double_t rescor[hi::NumEPNames];
Double_t rescorErr[hi::NumEPNames];
Double_t sumw[hi::NumEPNames];
TH2D * qxtrk_;
TH2D * qytrk_;
TH2D * qxtrk3_;
TH2D * qytrk3_;
TH2D * qcnt_;
TH2D * avpt_;

Int_t NumEvents[40];
Int_t TotNumEvents; 
TString KeyNames[40];
int NumKeys;
string reac_;
//----------------------------------
#include "src/GenSupport.h"
#include "src/Qvec.h"
#include "src/Setup.h"

void ReadTree(GetEventInfo * info, string prename, TString trig, TString trig2, TString trig3);
int GetMidIndx(Int_t epord, TString midn);
void  GetNumEvents(string prename, TString trig, TString trig2, TString trig3);
//void GenerateVN(string anal, TString reac = "pPb", TString prename = "rpflat160", TString trig = "crab_pPb_v2_HM160A", TString trig2 = "crab_pPb_v2_HM160B", TString trig3 = "", 
//TString mid2n =  "trackmid2", TString mid3n = "trackmid3"){
void GenerateVN(string anal="", TString reac="", int NTRK = 160) {

  string prename;
  TString trig;
  TString trig2;
  TString trig3;
  TString mid2n;
  TString mid3n;
#include "src/PbPbSetup.h"
#include "src/pPbSetup.h"

  cout<<"mid2n: "<<trig2<<endl;
  cout<<"mid3n: "<<trig3<<endl;
  ran = new TRandom3();
  reac_ = reac.Data();
  
  //
  //Locate information about data structure
  //
  GetEventInfo * info=0;
  for(int iset = 0; iset < 500; iset++) {
    TString inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig.Data(),prename.data(), iset);
    FILE *ftest = fopen(inFile.Data(),"r");
    if(ftest==NULL) continue;
    fclose(ftest);
    TFile * tf    = new TFile(inFile.Data(),"read");
    if(tf->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
    if(tf->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
    cout<<"Found: "<<inFile.Data()<<endl;
    info = new GetEventInfo(inFile);
    if(info->status == 0) {cout<<inFile.Data()<<" not found or has error"<<endl; return;}
    //    Setup(anal, info, mid2n, mid3n);
    Setup(anal, info);
    tf->Close();
    break;
  }
  //TString v3 = "";
  //if(epord_==3) v3 = "_v3";
  cout<<"Setup is complete epord_ set "<<epord_<<endl;
  GetNumEvents(prename, trig,trig2,trig3);
  for(int i = 0; i<NumKeys; i++) cout<<i<<"\t"<<KeyNames[i].Data()<<"\t"<<NumEvents[i]<<endl;
  ReadTree(info, prename, trig, trig2,trig3);
  OutputResults(reac, mid2n, trig);
}

void ReadTree(GetEventInfo * info, string prename, TString trig, TString trig2, TString trig3){ 
  trkbins->Reset();
  centbins->Reset();
  int ntrig = 1;
  if(trig2.Length()>0) ntrig = 2;
  if(trig3.Length()>0) ntrig = 3;
  TFile * tfin;
  int NumEvnts = 0;
  int nbins = ncentbins;
  if(trkoff) nbins=ntrkbins;
  int filecnt=0;
  if(RECENTERTRACKS) {
    for(int iset = 0; iset < 1000 && filecnt<MaxFiles; iset++) {
      for(int nt = 0; nt<ntrig; nt++){
      if(filecnt>=MaxFiles) continue;
	TString inFile="";
	if( nt==0) {
	  inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig.Data(),prename.data(),iset);
	} else if (nt==1) { 
	  inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig2.Data(),prename.data(),iset);
	} else if (nt==2) {
	  inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig3.Data(),prename.data(),iset);
	}
	FILE *ftest = fopen(inFile.Data(),"r");
	
	if(ftest==NULL) continue;
	fclose(ftest);
	tfin    = new TFile(inFile.Data(),"read");
	if(tfin->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
	if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
	tfin->ResetErrno();
	++filecnt;
	for(int i = 0; i<info->getNumKeys(); i++) {      

	  TTree * tree = (TTree * ) tfin->Get(Form("%s/tree",info->getKeyName(i).Data())); 
	  tree->SetBranchAddress("NtrkOff",    &noff);
	  tree->SetBranchAddress("Cent",  &centval);
	  tree->SetBranchAddress("Vtx",        &vtx);
	  tree->SetBranchAddress(Form("qxtrk_v%d",epord_),      &qxtrk_);
	  tree->SetBranchAddress(Form("qytrk_v%d",epord_),      &qytrk_);
	  if(ANAL==D2232SUB2 || ANAL==D2432SUB2||ANAL==D2232SUB3 || ANAL==D2432SUB3) {
	    tree->SetBranchAddress("qxtrk_v3",      &qxtrk3_);
	    tree->SetBranchAddress("qytrk_v3",      &qytrk3_);
	  }
	  for(int ievent = 0; ievent<tree->GetEntries(); ievent++) {
	    if(MaxEvents>0&&NumEvnts>=MaxEvents) {
	      cout<<"MaxEvents: "<<MaxEvents<<endl;
	      cout<<"NumEvnts: "<<NumEvnts<<endl;
	      continue;
	    }
	    tree->GetEntry(ievent);
	    if(fabs(vtx)>15.) continue;
	    int bin = -1;      
	    if(trkoff) {
	      bin =  trkbins->FindBin(noff)-1;
	      if(bin>=ntrkbins) continue;
	      if(bin<0) continue;
	    } else {
	      bin =  centbins->FindBin(centval)-1;
	      if(bin>=ncentbins) continue;
	      if(bin<0) continue;
	    }	 
	    qxav[bin]->Add(qxtrk_);
	    qyav[bin]->Add(qxtrk_);
	    ++NumEvnts;
	  }
	  tfin->Close();
	}
      }
    }
    for(int i = 0; i<nbins; i++) {
      if(qxav[i]) qxav[i]->Scale(1./(double) NumEvnts);
      if(qyav[i]) qyav[i]->Scale(1./(double) NumEvnts);
    }
  }
  int NEvt = TotNumEvents;
  cout<<"Nevt: "<<NEvt<<endl;
  NumEvnts = 0;
  filecnt = 0;
  for(int iset = 0; iset < 1000 && filecnt<MaxFiles; iset++) {
    for(int nt = 0; nt<ntrig; nt++){
      if(filecnt>=MaxFiles) continue;
      TString inFile="";
      if( nt==0) {
	inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig.Data(),prename.data(),iset);
      } else if (nt==1) { 
	inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig2.Data(),prename.data(),iset);
      } else if (nt==2) {
	inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig3.Data(),prename.data(),iset);
      }
      ispPb = kFALSE;
      if(inFile.Contains("pPb")) ispPb = kTRUE;
      string pPbTag = "Pbp";
      if(ispPb) pPbTag = "pPb";
      if(reac_.find("PbPb")!=std::string::npos) pPbTag="PbPb";
      FILE *ftest = fopen(inFile.Data(),"r");
      if(ftest==NULL) continue;
      fclose(ftest);
      tfin    = new TFile(inFile.Data(),"read");
      if(tfin->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
      if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
      tfin->ResetErrno();
      ++filecnt;
      for(int i = 0; i<info->getNumKeys(); i++) {      
	TTree * tree = (TTree * ) tfin->Get(Form("%s/tree",info->getKeyName(i).Data()));
	tree->SetBranchAddress("Cent",       &centval);
	tree->SetBranchAddress("NtrkOff",    &noff);
	tree->SetBranchAddress("Vtx",        &vtx);
	tree->SetBranchAddress("epang",      &epang);
	tree->SetBranchAddress("sumw", &sumw);
	tree->SetBranchAddress("qx",         &qx);
	tree->SetBranchAddress("qy",         &qy);
	tree->SetBranchAddress("q",          &q);
	tree->SetBranchAddress("mult",       &epmult);
	tree->SetBranchAddress("Run",        &runno_);
	tree->SetBranchAddress("Rescor",     &rescor);
	tree->SetBranchAddress("RescorErr",  &rescorErr);
	tree->SetBranchAddress(Form("qxtrk_v%d",epord_),      &qxtrk_);
	tree->SetBranchAddress(Form("qytrk_v%d",epord_),      &qytrk_);
	if(ANAL==D2232SUB2 || ANAL==D2432SUB2||ANAL==D2232SUB3 || ANAL==D2432SUB3) {
	  tree->SetBranchAddress("qxtrk_v3",      &qxtrk3_);
	  tree->SetBranchAddress("qytrk_v3",      &qytrk3_);
	}
	tree->SetBranchAddress(Form("qcnt_v%d",epord_),       &qcnt_);
	tree->SetBranchAddress("avpt",        &avpt_);
	for(int ievent = 0; ievent<tree->GetEntries(); ievent++) {
	  if(MaxEvents>0&&NumEvnts>=MaxEvents) {
	    cout<<"MaxEvents: "<<MaxEvents<<endl;
	    cout<<"NumEvnts: "<<NumEvnts<<endl;
	    continue;
	  }
	  tree->GetEntry(ievent);
	  if(fabs(vtx)>15.) continue;
	  int bin = -1; 
	  if(trkoff) {
	    bin =  trkbins->FindBin(noff)-1;
	    if(bin>=ntrkbins) continue;
	    if(bin<0) continue;
	  } else {
	    bin =  centbins->FindBin(centval)-1;
	    if(bin>=ncentbins) continue;
	    if(bin<0) continue;
	  }
	  if((int)fmod( NumEvnts, NEvt/20) == 0 ) cout<<(int) (100*(NumEvnts/(double)NEvt)+0.5)<<endl;
	  trkbins->Fill(noff);
	  centbins->Fill(centval);
	  int evtchar = centval;
	  if(trkoff) evtchar = noff;
	  int j=0;

	  if(ANAL==D2232SUB2 ||ANAL==D2232SUB3 ) {
	    TH2D * qxt = (TH2D *) qxtrk_->Clone("qxt");
	    TH2D * qyt = (TH2D *) qytrk_->Clone("qyt");
	    qxt->Reset();
	    qyt->Reset();
	    for(int ix = 1; ix<=qxt->GetNbinsX(); ix++) {
	      for(int jy = 1; jy<=qyt->GetNbinsY(); jy++) {
		double a = qxtrk_->GetBinContent(ix,jy);
		double b = qytrk_->GetBinContent(ix,jy);
		double c = qxtrk3_->GetBinContent(ix,jy);
		double d = qytrk3_->GetBinContent(ix,jy);
		qxt->SetBinContent(ix,jy,a*c-b*d);
		qyt->SetBinContent(ix,jy,b*c+a*d);
		qcnt_->SetBinContent(ix,jy,pow(qcnt_->GetBinContent(ix,jy),2));
	      }
	    }
	    genp[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    genm[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    j=(int)(ran->Uniform(0,9.999))+1;
	    genp[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    genm[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    qxt->Delete();
	    qyt->Delete();
	  } else if (ANAL==D2432SUB2||ANAL==D2432SUB3){
	    TH2D * qxt = (TH2D *) qxtrk_->Clone("qxt");
	    TH2D * qyt = (TH2D *) qytrk_->Clone("qyt");
	    qxt->Reset();
	    qyt->Reset();
	    for(int ix = 1; ix<=qxt->GetNbinsX(); ix++) {
	      for(int jy = 1; jy<=qyt->GetNbinsY(); jy++) {
		double a = qxtrk_->GetBinContent(ix,jy);
		double b = qytrk_->GetBinContent(ix,jy);
		double c = qxtrk3_->GetBinContent(ix,jy);
		double d = qytrk3_->GetBinContent(ix,jy);

		double cc = a*c - b*d;
		double dd = b*c + a*d;
		qxt->SetBinContent(ix,jy,a*cc-b*dd);
		qyt->SetBinContent(ix,jy,b*cc+a*dd);
		qcnt_->SetBinContent(ix,jy,pow(qcnt_->GetBinContent(ix,jy),3));
	      }
	    }
	    genp[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    genm[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    j=(int)(ran->Uniform(0,9.999))+1;
	    genp[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    genm[bin][j]->add(qxt,qyt,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    qxt->Delete();
	    qyt->Delete();

	  }else {
	    genp[bin][j]->add(qxtrk_,qytrk_,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    genm[bin][j]->add(qxtrk_,qytrk_,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    j=(int)(ran->Uniform(0,9.999))+1;
	    genp[bin][j]->add(qxtrk_,qytrk_,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	    genm[bin][j]->add(qxtrk_,qytrk_,qcnt_,qx,qy,sumw,avpt_,evtchar,qxav[bin],qyav[bin]);
	  }
	  ++NumEvnts;
	}
	tfin->Close();
      }
    }
  }
  //cout<<"Leaving readtree"<<endl;
  //  for(int i = 0; i<nbins; i++) {
  //  cout<<genp[i][0]->getCount()<<endl;
  //  cout<<genp[i][0]->getRescor()<<endl;
  //}
 
}
void GetNumEvents(string prename, TString trig, TString trig2, TString trig3){ 
  //Loop over datasets
  int ntrig = 1;
  TotNumEvents = 0;
  for(int i = 0; i<40; i++) NumEvents[i] = 0;
  if(trig2.Length()>0) ntrig = 2;
  if(trig3.Length()>0) ntrig = 3;
  TString tr = "";
  TFile * tfin;
  int filecnt = 0;
  for(int iset = 0; iset < 1000 && filecnt <= MaxFiles; iset++) {
    for(int nt = 0; nt<ntrig; nt++) {
      if(filecnt>=MaxFiles) continue;
      TString inFile;
      if(nt==0) {
	tr = trig;
	inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig.Data(),prename.data(),iset);
      } else if(nt==1){
 	tr = trig2;
	inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig2.Data(),prename.data(),iset);
      } else if (nt==2) {
 	tr = trig3;
	inFile = Form("/rfs/sanders/crab_projects/%s/%s_%d.root",trig3.Data(),prename.data(),iset);
      }
     FILE *ftest = fopen(inFile.Data(),"r");
      if(ftest==NULL) continue;
      fclose(ftest);
      TList * flist;
      int nkeys = 0;
      tfin    = new TFile(inFile.Data(),"read");
      if(tfin->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
      if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
      tfin->ResetErrno();
      flist = tfin->GetListOfKeys();
      ++filecnt;
      while(flist->At(nkeys) != flist->Last()) ++nkeys;
      ++nkeys;
      for(int i = 0; i<nkeys; i++) {
	TTree * tree = (TTree *) tfin->Get(Form("%s/tree",flist->At(i)->GetName()));
	if(NumKeys==0) KeyNames[i] = flist->At(i)->GetName();
	NumEvents[i] += tree->GetEntries();
	TotNumEvents += tree->GetEntries();
	cout<<KeyNames[i].Data()<<"\t"<<tree->GetEntries()<<"\t"<<NumEvents[i]<<"\t"<<inFile.Data()<<endl;
      }
      if(NumKeys==0) NumKeys = nkeys;
      tfin->Close();
    }
    //for(int i = 0; i<NumKeys; i++) {
    //  cout<<tr.Data()<<"\t"<<KeyNames[i].Data()<<"\t NumEvents: "<<NumEvents[i]<<endl;
    // }
  }
  cout<<"Total number of events: "<<TotNumEvents<<endl;
  return;
}
