#ifndef QREF
#define QREF
#include <complex>
#include <vector>
#include <string>
#include <iostream>
#include "TH2D.h"
#include "TH1D.h"

typedef complex<double> comp;
static comp zero (0,0);

class QRef {
 public:
  QRef(void){;}
  QRef(GetEventInfo * info, int epord, double minpt, double maxpt, double etamin1, 
       double etamax1, double etamin2, double etamax2, TH2D * eff, string name);
  QRef(GetEventInfo * info, int epord, double minpt, double maxpt, double etamin1, 
       double etamax1, double etamin2, double etamax2, TrackEfficiency * trkeff, string name);
  void Add(TH2D * qxtrk, TH2D * qytrk, TH2D * qcnt);
 private:
  int epord;
  int minetabin1;
  int maxetabin1;
  int minetabin2;
  int maxetabin2;
  int minptbin;
  int maxptbin;
  string name;
  comp Q12;
  int nevents;
  TH2D * eff;
  TrackEfficiency * trkeff;
  bool isPbPb;
};

QRef::QRef(GetEventInfo * info, int epn, double minpt, double maxpt, double etamin1, double etamax1, double etamin2, double etamax2, TH2D * heff, string sname){
  TH2D * temp = (TH2D *) info->getTemplate();
  minetabin1 = temp->GetYaxis()->FindBin(etamin1);
  maxetabin1 = temp->GetYaxis()->FindBin(etamax1-0.001);
  minetabin2 = temp->GetYaxis()->FindBin(etamin2);
  maxetabin2 = temp->GetYaxis()->FindBin(etamax2-0.001);
  minptbin = temp->GetXaxis()->FindBin(minpt);
  maxptbin = temp->GetXaxis()->FindBin(maxpt-0.001);
  std::cout<<minetabin1<<"\t"<<maxetabin1<<"\t"<<minetabin2<<"\t"<<maxetabin2<<"\t"<<minptbin<<"\t"<<maxptbin<<std::endl;
  name = sname;
  epord = epn;
  Q12 = zero;
  nevents = 0;
  eff = heff;
  isPbPb= false;
}
QRef::QRef(GetEventInfo * info, int epn, double minpt, double maxpt, double etamin1, double etamax1, double etamin2, double etamax2, TrackEfficiency * teff, string sname){
  TH2D * temp = (TH2D *) info->getTemplate();
  minetabin1 = temp->GetYaxis()->FindBin(etamin1);
  maxetabin1 = temp->GetYaxis()->FindBin(etamax1-0.001);
  minetabin2 = temp->GetYaxis()->FindBin(etamin2);
  maxetabin2 = temp->GetYaxis()->FindBin(etamax2-0.001);
  minptbin = temp->GetXaxis()->FindBin(minpt);
  maxptbin = temp->GetXaxis()->FindBin(maxpt-0.001);
  std::cout<<minetabin1<<"\t"<<maxetabin1<<"\t"<<minetabin2<<"\t"<<maxetabin2<<"\t"<<minptbin<<"\t"<<maxptbin<<std::endl;
  name = sname;
  epord = epn;
  Q12 = zero;
  nevents = 0;
  trkeff = teff;
  isPbPb = true;
}

void QRef::Add(TH2D * qxtrk, TH2D * qytrk, TH2D * qcnt){
  comp q1 = zero;
  comp q2 = zero;
  int summult1 = 0;
  int summult2 = 0;
  for(int i = minptbin; i<=maxptbin; i++) {
    for(int j = minetabin1; j<=maxetabin1; j++){
      comp val = (qxtrk->GetBinContent(i,j),qytrk->GetBinContent(i,j));
      val*=qcnt->GetBinContent(i,j);
      q1+=val;
      summult1+=qcnt->GetBinContent(i,j);
    }
    for(int j = minetabin2; j<=maxetabin2; j++){
      comp val = (qxtrk->GetBinContent(i,j),qytrk->GetBinContent(i,j));
      val*=qcnt->GetBinContent(i,j);
      q2+=val;
      summult2+=qcnt->GetBinContent(i,j);
    }
  }
}
#endif

