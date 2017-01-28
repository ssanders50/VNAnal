
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

using namespace std;

class MCEvent{

  public:
        
  MCEvent(Double_t tv1,Double_t tv2,Double_t tv3,Double_t tv4,Double_t tv5,Double_t tv6,Double_t tv7);
    void SetEventParms();
    void SetMult(Int_t value);
    void SetPsiRandom();
    void SetPsi(Double_t value){Psi = value;}
    void GetPtRandom(Double_t ptTrackArray[]);
    void Setv1(Double_t value);
    void Setv2(Double_t value);
    void Setv3(Double_t value);
    void Setv4(Double_t value);
    void Setv5(Double_t value);
    void Setv6(Double_t value);
    Int_t GetMult();
    void GetThrowPhi(Double_t phiTrackArray[]);
    Double_t GetPsi();
    void SetSeed(Int_t seed);

  private:
    
    Double_t v1;
    Double_t v2;
    Double_t v3;
    Double_t v4;
    Double_t v5;
    Double_t v6;
    Double_t v7;
    Double_t Psi;
    Int_t mult;
    TRandom * ran;
    TF1 * phidist;
    TF1 * ptdist;

};

//-----------------------------------------------------------------------------
//                                  METHODS
//-----------------------------------------------------------------------------

MCEvent::MCEvent(Double_t tv1,Double_t tv2,Double_t tv3,Double_t tv4,Double_t tv5,Double_t tv6, Double_t tv7){
  ran = new TRandom3(0);
  v1 = tv1;
  v2 = tv2;
  v3 = tv3;
  v4 = tv4;
  v5 = tv5;
  v6 = tv6;
  v7 = tv7;
  phidist = new TF1("phidist","1+2*[0]*cos(x)+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x)+2*[5]*cos(6*x)+2*[6]*cos(7*x)",-TMath::Pi(),TMath::Pi());
  phidist->SetParameters(v1,v2,v3,v4,v5,v6);
  ptdist = new TF1("ptdist","x*((([0]-1)*([0]-2))/([1]*[1]))*[2]*TMath::Power(1+x/[1],-[0])",0,12);
  ptdist->SetParameters(6.06271, 1.08372, 586.474);

    
}

//-----------------------------------------------------------------------------
void MCEvent::SetMult(Int_t value) {mult=value;}
void MCEvent::SetPsiRandom(){ 
  Psi = ran->Uniform(-TMath::Pi(),TMath::Pi());    
}

void MCEvent::SetEventParms() {    
  phidist->SetParameters(v1,v2,v3,v4,v5,v6,v7);   
 
}
//-----------------------------------------------------------------------------

Int_t MCEvent::GetMult() {return mult;}

//-----------------------------------------------------------------------------

void MCEvent::GetPtRandom(Double_t ptTrackArray[]) {
    
  for(Int_t i=0; i<mult;i++){
    ptTrackArray[i] = ptdist->GetRandom();
  }
    
}

//-----------------------------------------------------------------------------

Double_t MCEvent::GetPsi() {
    
  return Psi;
    
}

//-----------------------------------------------------------------------------

void MCEvent::GetThrowPhi(Double_t phiTrackArray[]){
  for(Int_t i=0;i<mult;i++){
    //    Double_t phi = phidist->GetRandom()+Psi;
    Double_t phi = phidist->GetRandom();
    phiTrackArray[i] = phi;
  }
}

//-----------------------------------------------------------------------------

void MCEvent::Setv1(Double_t value) { v1 = value; }
void MCEvent::Setv2(Double_t value) { v2 = value; }
void MCEvent::Setv3(Double_t value) { v3 = value; }
void MCEvent::Setv4(Double_t value) { v4 = value; }
void MCEvent::Setv5(Double_t value) { v5 = value; }
void MCEvent::Setv6(Double_t value) { v6 = value; }
void MCEvent::SetSeed(Int_t seed) {
  ran->SetSeed(seed);
  gRandom->SetSeed(seed);
}
