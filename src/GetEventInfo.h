#ifndef GETEVENTINFO
#define GETEVENTINFO
#include "TString.h"
#include "TList.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TArrayD.h"
#include <string>

class GetEventInfo{
 public:
  GetEventInfo(TString inFile);
  ~GetEventInfo(){tfin->Delete();}
  int getNumKeys(){return NumKeys;}
  int getKeyCount(int key) {return NumEvents[key];}
  TString getKeyName(int key) {return KeyNames[key];}
  TTree * getTree(int key){ return tree[key]; }
  int status ;
  int getNumEtaBins(){return netabins;}
  int getNumPtBins(){return nptbins;}
  TH1D * getEtaHist(){return heta;}
  TH1D * getPtHist(){return hpt;}
  TH2D * getTemplate(){return hTemplate;}
  TString getFileName(){return filename;}
  bool ispPb2011(){return pPb2011;}
  bool isPbPb2015pixel(){return PbPb2015pixel;} 
 private:
  bool pPb2011;
  bool PbPb2015pixel;
  int NumKeys;
  int NumEvents[20];
  TString KeyNames[20];
  TTree * tree[20];
  int netabins;
  int nptbins;
  TH1D * heta;
  TH1D * hpt;
  TH2D * hTemplate;
  TFile * tfin;
  TString filename;
};

GetEventInfo::GetEventInfo(TString inFile){
  pPb2011 = false;
  PbPb2015pixel=false;
  if(inFile.Contains("2011")) pPb2011 = true;
  if(inFile.Contains("pPb"))  pPb2011 = true;
  if(inFile.Contains("2015")) PbPb2015pixel = true; 
  if(pPb2011) cout<<"pPb2011 set true"<<endl;
  if(PbPb2015pixel) cout<<"PbPb2015pixel set true"<<endl;
  status = 0;
  NumKeys = 0;
  filename = inFile;
  hTemplate = 0;
  FILE *ftest = fopen(inFile.Data(),"r");
  if(ftest==NULL) return;
  fclose(ftest);
  TList * flist;
  int nkeys = 0;
  tfin    = new TFile(inFile.Data(),"read");
  if(tfin->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; return;}
  if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; return;}
  tfin->ResetErrno();
  flist = tfin->GetListOfKeys();

  while(flist->At(nkeys) != flist->Last()) ++nkeys;
  ++nkeys;
  cout<<"nkeys: "<<nkeys<<endl;
  const TArrayD * etabins;
  const TArrayD * ptbins;
  for(int i = 0; i<nkeys; i++) {
    tree[i] = (TTree *) tfin->Get(Form("%s/tree",flist->At(i)->GetName()));
    cout<<"tree: "<<Form("%s/tree",flist->At(i)->GetName())<<endl;
    if(NumKeys==0) KeyNames[i] = flist->At(i)->GetName();
    NumEvents[i] = tree[i]->GetEntries();
    cout<<"Key "<<i<<"\t"<<KeyNames[i].Data()<<endl;

    TH2D * h2 = (TH2D *) tfin->Get(Form("%s/qxtrk_v2",flist->At(i)->GetName()));
    if(h2 == 0) {
      h2 = (TH2D *) tfin->Get(Form("%s/qxtrk_v3",flist->At(i)->GetName()));
    }
    if(h2 == 0) {
      h2 = (TH2D *) tfin->Get(Form("%s/qxtrk_v4",flist->At(i)->GetName()));
    }
    if(h2 == 0) {
      h2 = (TH2D *) tfin->Get(Form("%s/qxtrk_v5",flist->At(i)->GetName()));
    }
    if(h2 == 0) {
      h2 = (TH2D *) tfin->Get(Form("%s/qxtrk_v6",flist->At(i)->GetName()));
    }
    if(h2 == 0) {
      h2 = (TH2D *) tfin->Get(Form("%s/qxtrk_v7",flist->At(i)->GetName()));
    }
    if(h2==0) return;
    netabins = h2->GetNbinsY();
    nptbins = h2->GetNbinsX();
    etabins = h2->GetYaxis()->GetXbins();
    ptbins =  h2->GetXaxis()->GetXbins();
    if(!hTemplate) hTemplate = (TH2D *) h2->Clone("hTemplate");
  }
  if(hTemplate) hTemplate->Reset();
  if(NumKeys==0) NumKeys = nkeys;
  double pttmp[100];
  double etatmp[100];

  for(int i = 0; i<=nptbins; i++) {
    pttmp[i] = ptbins->At(i);
  }
  for(int i = 0; i<=netabins; i++) {
    etatmp[i] = etabins->At(i);
  }
  hpt = new TH1D(Form("hpt_%s",inFile.Data()),"hpt",nptbins,pttmp);
  heta = new TH1D(Form("heta_%s",inFile.Data()),"heta",netabins,etatmp);
  std::cout<<"Replay Setup: "<<std::endl;
  std::cout<<"   pt bins: "<<nptbins<<std::endl;
  std::cout<<"  eta bins: "<<netabins<<std::endl; 
  status = 1;
  return;
}

#endif
