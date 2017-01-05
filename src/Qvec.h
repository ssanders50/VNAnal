#ifndef QVEC
#define QVEC
#include <complex>
#include <vector>
typedef complex<double> comp;
comp zero (0,0);

class Qvec {
 public:
  Qvec(void){;}
  Qvec(int ptbin, int etabin, string name);
  void add(int mult, comp Qn, 
	   comp Q2A, double w2A, comp Q2B, double w2B,  comp Q2C, double w2C, 
	   comp Q3A, double w3A, comp Q3B, double w3B,  comp Q3C, double w3C, 
	   double pt, int noff); 
  double getRescor();
  double getRescorErr();
  double getEPVNobs();
  double getEPVNobsErr();
  double getSPVN();
  double getSPVNErr();
  double getSPVNnum();
  double getSPVNnumerr();
  double getSPVNdenom();
  double getSPVNdenomerr();
  double getEPVN();
  double getEPVNErr();
  double getAvPt();
  double getAvPtErr();
  double getAvNoff();
  double getAvNoffErr();
  double getnfrac(){return (double)(totevents-multlost)/(double)totevents;}
  double getafrac(){return (double)(qatot-qalost)/(double)qatot;}
  double getbfrac(){return (double)(qbtot-qblost)/(double)qbtot;}
  double getcfrac(){return (double)(qctot-qclost)/(double)qctot;}
  double getQAB() {return QAB.real()/wABsum;}
  double getQAC() {return QAC.real()/wACsum;}
  double getQBC() {return QBC.real()/wBCsum;}
  double getQnA() {return QnA.real()/wnAsum;}
  int getCount(){return ncount;}
  int getNevts(){return nevts;}
  string getName(){return name_;}
 private:
  comp QAB;
  comp QAC;
  comp QBC;
  comp QnA;
  comp QABnorm;
  comp QACnorm;
  comp QBCnorm;
  comp QnAnorm;

  double wnAsum;
  double wABsum;
  double wACsum;
  double wBCsum;
  double wnA2sum;
  double wAB2sum;
  double wAC2sum;
  double wBC2sum;


  double rQAB;
  double rQAC;
  double rQBC;
  double rQnA;
  double rQABnorm;
  double rQACnorm;
  double rQBCnorm;
  double rQnAnorm;
  double rQABnorm2;
  double rQACnorm2;
  double rQBCnorm2;
  double rQnAnorm2;

  double rQAB2;
  double rQAC2;
  double rQBC2;
  double rQnA2;

  string name_;
  int ptbin_;
  int etabin_;
  int nevts;
  int ncount;
  double rescor;
  double avnoff;
  double avpt;
  double avnoff2;
  double avpt2;
  int multlost;
  int qalost;
  int qblost;
  int qclost;
  int totevents;
  int qatot;
  int qbtot;
  int qctot;
  bool twosub;
};

Qvec::Qvec(int ptbin, int etabin, string name){
  ptbin_ = ptbin;
  etabin_ = etabin;
  name_ = name;
  nevts = 0;
  ncount = 0;
  wnAsum=0;
  wABsum=0;
  wACsum=0;
  wBCsum=0;
  wnA2sum=0;
  wAB2sum=0;
  wAC2sum=0;
  wBC2sum=0;

  QAB = zero;
  QAC = zero;
  QBC = zero;
  QABnorm = zero;
  QACnorm = zero;
  QBCnorm = zero;
  QnAnorm = zero;

  rQAB = 0;
  rQAC = 0;
  rQBC = 0;
  rQABnorm = 0;
  rQACnorm = 0;
  rQBCnorm = 0;
  rQnAnorm = 0;
  rQABnorm2 = 0;
  rQACnorm2 = 0;
  rQBCnorm2 = 0;
  rQnAnorm2 = 0;

  rQAB2 = 0;
  rQAC2 = 0;
  rQBC2 = 0;
  totevents = 0;
  qatot = 0;
  qbtot = 0;
  qctot = 0;
  multlost = 0;
  qalost = 0;
  qblost = 0;
  qclost = 0;
  rescor = -1;
  avpt = 0;
  twosub = false;
  if(ANAL == N2SUB2 || ANAL == N3SUB2 || ANAL == N4SUB2 || ANAL==N42SUB2 ||
     ANAL == N5SUB2 || ANAL == N6SUB2 || ANAL == N7SUB2 || ANAL==N523SUB2 ||
     ANAL == N63SUB2|| ANAL == N723SUB2 || ANAL == D24SUB2 || ANAL == D34SUB2 ||
     ANAL == D2232SUB2 || ANAL == D2432SUB2) {
    twosub = true;
  }
}

void Qvec::add(int imult, comp Qn, 
	       comp Q2A, double w2A, comp Q2B, double w2B,  comp Q2C, double w2C, 
	       comp Q3A, double w3A, comp Q3B, double w3B,  comp Q3C, double w3C, 
	       double pt, int noff) {
  ++totevents;
  if(imult==0) {
    ++multlost;
    return;
  }
  ++qatot;
  if(std::abs(Q2A) < 1e-5) {
    ++qalost;
    return;
  }
  ++qbtot;
  if(std::abs(Q2B) < 1e-5) {
    ++qblost;
    return;
  }
  ++qctot;
  if(std::abs(Q2C) < 1e-5) {
    ++qclost;
    return;
  }
  //
  // Commented is the correct code with weights, the
  // actual code simplifies
  //
  /* AB = QA * std::conj(QB)/(wA*wB); */
  /* AC = QA * std::conj(QC)/(wA*wC); */
  /* BC = QB * std::conj(QC)/(wB*wC);   */
  /* nA = Qn * std::conj(QA)/(wA*(double)mult) ; */
  /* QAB += wA*wB*AB; */
  /* QAC += wA*wC*AC; */
  /* QBC += wB*wC*BC; */
  /* QnA += ((double)mult)*wA*nA; */
  /* rQAB += wA*wB*AB.real(); */
  /* rQAC += wA*wC*AC.real(); */
  /* rQBC += wB*wC*BC.real(); */
  /* rQnA += ((double)mult)*wA*nA.real(); */
  /* rQAB2 += wA*wB*pow(AB.real(),2); */
  /* rQAC2 += wA*wC*pow(AC.real(),2); */
  /* rQBC2 += wB*wC*pow(BC.real(),2); */
  /* rQnA2 += ((double)mult)*wA*pow(nA.real(),2); */

  comp unit(1.,0);
  comp AB;
  comp AC;
  comp BC;
  comp nA;
  comp ABnorm(0.,0.);
  comp ACnorm(0.,0.);
  comp BCnorm(0.,0.);
  comp nAnorm(0.,0.);
  //  wACsum = 1;
  // wBCsum = 1;
  // wAC2sum = 1;
  // wBC2sum = 1;
  // AC = 1;
  double mult = (double) imult;
  if(ANAL == N2SUB2 || ANAL == N3SUB2 || ANAL == N4SUB2 ||
     ANAL == N5SUB2 || ANAL == N6SUB2 || ANAL == N7SUB2) {
    wnAsum+=w2A*mult;
    wABsum+=w2A*w2B;
    wnA2sum+=pow(w2A*mult,2);
    wAB2sum+=pow(w2A*w2B,2);
    AB = Q2A * std::conj(Q2B);
    nA = Qn * std::conj(Q2A) ;
    QAB += AB;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2B);
    rQnA2 += pow(nA.real(),2)/(mult*w2A);
    ABnorm = (Q2A/std::abs(Q2A)) * (std::conj(Q2B)/std::abs(Q2B));
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A));

  } else if(ANAL == N42SUB2 || ANAL == N63SUB2) {
    wnAsum+=w2A*w2A*mult;
    wABsum+=w2A*w2A*w2B*w2B;
    wnA2sum+=pow(w2A*w2A*mult,2);
    wAB2sum+=pow(w2A*w2A*w2B*w2B,2);
    AB = Q2B * Q2B* std::conj(Q2A) * std::conj(Q2A);
    nA = Qn * std::conj(Q2A)*std::conj(Q2A) ;
    QAB += AB;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2A*w2B*w2B);
    rQnA2 += pow(nA.real(),2)/(mult*w2A*w2A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) * (std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q2A)/std::abs(Q2A));
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A));

  } else if(ANAL == D24SUB2 || ANAL == D34SUB2) {
    wnAsum+=w2A*w2A*mult*mult;
    wABsum+=w2A*w2A*w2B*w2B;
    wnA2sum+=pow(w2A*w2A*mult*mult,2);
    wAB2sum+=pow(w2A*w2A*w2B*w2B,2);
    AB = Q2B * Q2B* std::conj(Q2A) * std::conj(Q2A);
    nA = Qn * Qn*std::conj(Q2A)*std::conj(Q2A) ;
    QAB += AB;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2A*w2B*w2B);
    rQnA2 += pow(nA.real(),2)/(mult* mult* w2A*w2A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) * (std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q2A)/std::abs(Q2A));
    nAnorm = Qn* Qn *(std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A));

  } else if(ANAL == N523SUB2 || ANAL == D2232SUB2 || ANAL == D2432SUB2) {
    wnAsum+=w2A*w3A*mult;
    wABsum+=w2A*w3A*w2B*w3B;
    wnA2sum+=pow(w2A*w3A*mult,2);
    wAB2sum+=pow(w2A*w3A*w2B*w3B,2);
    AB = Q2B * Q3B* std::conj(Q2A) * std::conj(Q3A);
    nA = Qn * std::conj(Q2A)*std::conj(Q3A) ;
    QAB += AB;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w3A*w2B*w3B);
    rQnA2 += pow(nA.real(),2)/( mult*w2A*w3A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q3B/std::abs(Q3B)) * (std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q3A)/std::abs(Q3A));
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q3A)/std::abs(Q3A));

  } else if(ANAL == D2432SUB2) {
    wnAsum+=w2A*w2A*w3A*mult;
    wABsum+=w2A*w2A*w3A*w2B*w2B*w3B;
    wnA2sum+=pow(w2A*w2A*w3A*mult,2);
    wAB2sum+=pow(w2A*w2A*w3A*w2B*w2B*w3B,2);
    AB = Q2B * Q2B * Q3B* std::conj(Q2A) * std::conj(Q2A)*std::conj(Q3A);
    nA = Qn * std::conj(Q2A)*std::conj(Q2A)* std::conj(Q3A) ;
    QAB += AB;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2A*w3A*w2B*w2B*w3B);
    rQnA2 += pow(nA.real(),2)/( mult*w2A*w2A*w3A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) *(Q3B/std::abs(Q3B)) * (std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q3A)/std::abs(Q3A));
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q3A)/std::abs(Q3A));

  } else if(ANAL == N723SUB2) {
    wnAsum+=w2A*w2A*w3A*mult;
    wABsum+=w2A*w2A*w3A*w2B*w2B*w3B;
    wnA2sum+=pow(w2A*w2A*w3A*mult,2);
    wAB2sum+=pow(w2A*w2A*w3A*w2B*w2B*w3B,2);
    AB = Q2B * Q2B * Q3B* std::conj(Q2A) * std::conj(Q2A)* std::conj(Q3A);
    nA = Qn * std::conj(Q2A)*std::conj(Q2A)*std::conj(Q3A) ;
    QAB += AB;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2A*w3A*w2B*w2B*w3B);
    rQnA2 += pow(nA.real(),2)/(mult*w2A*w2A*w3A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) *(Q3B/std::abs(Q3B)) * (std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q2A)/std::abs(Q2A)) *(std::conj(Q3A)/std::abs(Q3A));
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q3A)/std::abs(Q3A));

  } else if(ANAL == N2SUB3 || ANAL == N3SUB3 || ANAL == N4SUB3 || 
	    ANAL == N5SUB3 || ANAL == N6SUB3 || ANAL == N7SUB3) {
    wnAsum+=w2A*mult;
    wABsum+=w2A*w2B;
    wACsum+=w2A*w2C;
    wBCsum+=w2B*w2C;
    wnA2sum+=pow(w2A*mult,2);
    wAB2sum+=pow(w2A*w2B,2);
    wAC2sum+=pow(w2A*w2C,2);
    wBC2sum+=pow(w2B*w2C,2);
    AB = Q2A * std::conj(Q2B);
    AC = Q2A * std::conj(Q2C);
    BC = Q2B * std::conj(Q2C);  
    nA = Qn * std::conj(Q2A) ;
    QAB += AB;
    QAC += AC;
    QBC += BC;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2B);
    rQAC2 += pow(AC.real(),2)/(w2A*w2C);
    rQBC2 += pow(BC.real(),2)/(w2B*w2C);
    rQnA2 += pow(nA.real(),2)/(mult*w2A);
    ABnorm = (Q2B/std::abs(Q2B)) * (std::conj(Q2A)/std::abs(Q2A));
    ACnorm = (Q2C/std::abs(Q2C)) * (std::conj(Q2A)/std::abs(Q2A));
    BCnorm = (Q2C/std::abs(Q2C)) * (std::conj(Q2B)/std::abs(Q2B));  
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A));

  } else if(ANAL == N42SUB3 || ANAL == N63SUB3) {
    wnAsum+=w2A*w2A*mult;
    wABsum+=w2A*w2A*w2B*w2B;
    wACsum+=w2A*w2A*w2C*w2C;
    wBCsum+=w2B*w2B*w2C*w2C;
    wnA2sum+=pow(w2A*w2A*mult,2);
    wAB2sum+=pow(w2A*w2B,4);
    wAC2sum+=pow(w2A*w2C,4);
    wBC2sum+=pow(w2B*w2C,4);
    AB = Q2B * Q2B * std::conj(Q2A) * std::conj(Q2A);
    AC = Q2C * Q2C * std::conj(Q2A) * std::conj(Q2A);
    BC = Q2C * Q2C * std::conj(Q2B) * std::conj(Q2B);  
    nA = Qn * std::conj(Q2A) * std::conj(Q2A) ;
    QAB += AB;
    QAC += AC;
    QBC += BC;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/pow(w2A*w2B,2);
    rQAC2 += pow(AC.real(),2)/pow(w2A*w2C,2);
    rQBC2 += pow(BC.real(),2)/pow(w2B*w2C,2);
    rQnA2 += pow(nA.real(),2)/(mult*w2A*w2A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q2A)/std::abs(Q2A));
    ACnorm = (Q2C/std::abs(Q2C)) *(Q2C/std::abs(Q2C)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q2A)/std::abs(Q2A));
    BCnorm = (Q2C/std::abs(Q2C)) *(Q2C/std::abs(Q2C)) * (std::conj(Q2B)/std::abs(Q2B))* (std::conj(Q2B)/std::abs(Q2B));  
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q2A)/std::abs(Q2A));

  } else if(ANAL == D24SUB3 || ANAL == D34SUB3) {
    wnAsum+=w2A*w2A*mult*mult;
    wABsum+=w2A*w2A*w2B*w2B;
    wACsum+=w2A*w2A*w2C*w2C;
    wBCsum+=w2B*w2B*w2C*w2C;
    wnA2sum+=pow(w2A*w2A*mult*mult,2);
    wAB2sum+=pow(w2A*w2B,4);
    wAC2sum+=pow(w2A*w2C,4);
    wBC2sum+=pow(w2B*w2C,4);
    AB = Q2B * Q2B * std::conj(Q2A) * std::conj(Q2A);
    AC = Q2C * Q2C * std::conj(Q2A) * std::conj(Q2A);
    BC = Q2C * Q2C * std::conj(Q2B) * std::conj(Q2B);  
    nA = Qn * Qn * std::conj(Q2A) * std::conj(Q2A) ;
    QAB += AB;
    QAC += AC;
    QBC += BC;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/pow(w2A*w2B,2);
    rQAC2 += pow(AC.real(),2)/pow(w2A*w2C,2);
    rQBC2 += pow(BC.real(),2)/pow(w2B*w2C,2);
    rQnA2 += pow(nA.real(),2)/(mult*mult*w2A*w2A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q2A)/std::abs(Q2A));
    ACnorm = (Q2C/std::abs(Q2C)) *(Q2C/std::abs(Q2C)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q2A)/std::abs(Q2A));
    BCnorm = (Q2C/std::abs(Q2C)) *(Q2C/std::abs(Q2C)) * (std::conj(Q2B)/std::abs(Q2B))* (std::conj(Q2B)/std::abs(Q2B));  
    nAnorm = Qn * Qn* (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q2A)/std::abs(Q2A));

  } else if(ANAL == N523SUB3 || ANAL == D2232SUB3 ) {
    wnAsum+=w2A*w3A*mult;
    wABsum+=w2A*w3A*w2B*w3B;
    wACsum+=w2A*w3A*w2C*w3C;
    wBCsum+=w2B*w3B*w2C*w3C;
    wnA2sum+=pow(w2A*w3A*mult,2);
    wAB2sum+=pow(w2A*w3A*w2B*w3B,2);
    wAC2sum+=pow(w2A*w3A*w2C*w3C,2);
    wBC2sum+=pow(w2B*w3B*w2C*w3C,2);
    AB = Q2B * Q3B * std::conj(Q2A) * std::conj(Q3A);
    AC = Q2C * Q3C * std::conj(Q2A) * std::conj(Q3A);
    BC = Q2C * Q3C * std::conj(Q2B) * std::conj(Q3B);  
    nA = Qn * std::conj(Q2A) * std::conj(Q3A) ;
    QAB += AB;
    QAC += AC;
    QBC += BC;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w3A*w2B*w3B);
    rQAC2 += pow(AC.real(),2)/(w2A*w3A*w2C*w3C);
    rQBC2 += pow(BC.real(),2)/(w2B*w3B*w2C*w3C);
    rQnA2 += pow(nA.real(),2)/(mult*w2A*w3A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q3B/std::abs(Q3B)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
    ACnorm = (Q2C/std::abs(Q2C)) *(Q3C/std::abs(Q3C)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
    BCnorm = (Q2C/std::abs(Q2C)) *(Q3C/std::abs(Q3C)) * (std::conj(Q2B)/std::abs(Q2B))* (std::conj(Q3B)/std::abs(Q3B));  
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));

  } else if(ANAL==D2432SUB3) {
    wnAsum+=w2A*w2A*w3A*mult;
    wABsum+=w2A*w2A*w3A*w2B*w2B*w3B;
    wACsum+=w2A*w2A*w3A*w2C*w2C*w3C;
    wBCsum+=w2B*w2B*w3B*w2C*w2C*w3C;
    wnA2sum+=pow(w2A*w2A*w3A*mult,2);
    wAB2sum+=pow(w2A*w2A*w3A*w2B*w2B*w3B,2);
    wAC2sum+=pow(w2A*w2A*w3A*w2C*w2C*w3C,2);
    wBC2sum+=pow(w2B*w2B*w3B*w2C*w2C*w3C,2);
    AB = Q2B * Q2B * Q3B * std::conj(Q2A) * std::conj(Q2A) * std::conj(Q3A);
    AC = Q2C * Q2C * Q3C * std::conj(Q2A) * std::conj(Q2A) * std::conj(Q3A);
    BC = Q2C * Q2C * Q3C * std::conj(Q2B) *std::conj(Q2B) * std::conj(Q3B);  
    nA = Qn * std::conj(Q2A) *std::conj(Q2A) * std::conj(Q3A) ;
    QAB += AB;
    QAC += AC;
    QBC += BC;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2A*w3A*w2B*w2B*w3B);
    rQAC2 += pow(AC.real(),2)/(w2A*w2A*w3A*w2C*w2C*w3C);
    rQBC2 += pow(BC.real(),2)/(w2B*w2B*w3B*w2C*w2C*w3C);
    rQnA2 += pow(nA.real(),2)/(mult*w2A*w2A*w3A);
    ABnorm = (Q2B/std::abs(Q2B)) * (Q2B/std::abs(Q2B)) * (Q3B/std::abs(Q3B)) * (std::conj(Q2A)/std::abs(Q2A)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
    ACnorm = (Q2C/std::abs(Q2C)) * (Q2C/std::abs(Q2C)) * (Q3C/std::abs(Q3C)) * (std::conj(Q2A)/std::abs(Q2A)) * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
    BCnorm = (Q2C/std::abs(Q2C)) * (Q2C/std::abs(Q2C)) * (Q3C/std::abs(Q3C)) * (std::conj(Q2B)/std::abs(Q2B)) * (std::conj(Q2B)/std::abs(Q2B))* (std::conj(Q3B)/std::abs(Q3B));  
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));

  } else if(ANAL == N723SUB3 ) {
    wnAsum+=w2A*w2A*w3A*mult;
    wABsum+=w2A*w2A*w3A*w2B*w2B*w3B;
    wACsum+=w2A*w2A*w3A*w2C*w2C*w3C;
    wBCsum+=w2B*w2B*w3B*w2C*w2C*w3C;
    wnA2sum+=pow(w2A*w2A*w3A*mult,2);
    wAB2sum+=pow(w2A*w2A*w3A*w2B*w2B*w3B,2);
    wAC2sum+=pow(w2A*w2A*w3A*w2C*w2C*w3C,2);
    wBC2sum+=pow(w2B*w2B*w3B*w2C*w2C*w3C,2);
    AB = Q2B * Q2B *Q3B * std::conj(Q2A) * std::conj(Q2A) * std::conj(Q3A);
    AC = Q2C * Q2C *Q3C * std::conj(Q2A) * std::conj(Q2A) * std::conj(Q3A);
    BC = Q2C * Q2C *Q3C * std::conj(Q2B) * std::conj(Q2B) * std::conj(Q3B);  
    nA = Qn * std::conj(Q2A) *std::conj(Q2A) * std::conj(Q3A) ;
    QAB += AB;
    QAC += AC;
    QBC += BC;
    QnA += nA;
    rQAB2 += pow(AB.real(),2)/(w2A*w2A*w3A*w2B*w2B*w3B);
    rQAC2 += pow(AC.real(),2)/(w2A*w2A*w3A*w2C*w2C*w3C);
    rQBC2 += pow(BC.real(),2)/(w2B*w2B*w3B*w2C*w2C*w3C);
    rQnA2 += pow(nA.real(),2)/(mult*w2A*w2A*w3A);
    ABnorm = (Q2B/std::abs(Q2B)) *(Q2B/std::abs(Q2B)) *(Q3B/std::abs(Q3B)) * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
    ACnorm = (Q2C/std::abs(Q2C)) *(Q2C/std::abs(Q2C)) *(Q3C/std::abs(Q3C)) * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
    BCnorm = (Q2C/std::abs(Q2C)) *(Q2C/std::abs(Q2C)) *(Q3C/std::abs(Q3C)) * (std::conj(Q2B)/std::abs(Q2B))*(std::conj(Q2B)/std::abs(Q2B))* (std::conj(Q3B)/std::abs(Q3B));  
    nAnorm = Qn * (std::conj(Q2A)/std::abs(Q2A))*(std::conj(Q2A)/std::abs(Q2A))* (std::conj(Q3A)/std::abs(Q3A));
  }

  QABnorm += ABnorm;
  QACnorm += ACnorm;
  QBCnorm += BCnorm;
  QnAnorm += nAnorm;

  rQABnorm += ABnorm.real();
  rQACnorm += ACnorm.real();
  rQBCnorm += BCnorm.real();
  rQnAnorm += nAnorm.real();
  rQABnorm2 += pow(ABnorm.real(),2);
  rQACnorm2 += pow(ACnorm.real(),2);
  rQBCnorm2 += pow(BCnorm.real(),2);
  rQnAnorm2 += pow(nAnorm.real(),2);
  
  ++nevts;
  ncount+=mult;
  avpt+=pt;
  avpt2+=pow(pt,2.);
  avnoff+=noff;
  avnoff2+=pow(noff,2.);
  return;
}

double Qvec::getRescor(){
  rescor = 0;
  double ab = 0;
  if( twosub ) {
    if(nevts>0){ 
      double ab = QABnorm.real()/(double) nevts;
      if(ab>0) {
	rescor = sqrt( ab );
      }  
    }
  } else {
    if(nevts>0) {
      double ab = QABnorm.real()/(double) nevts;
      double ac = QACnorm.real()/(double) nevts;
      double bc = QBCnorm.real()/(double) nevts;
      if(bc>0 && ab>0 && ac>0) {
	rescor = sqrt( ab*ac/bc );
      } 
    }
  }
    return rescor;
}


double Qvec::getRescorErr(){
  if( twosub ) {
    double muab = rQABnorm/(double) nevts;
    double varab = rQABnorm2/(double)nevts - pow(muab,2.);
    double sigab = sqrt(varab/(double) nevts);
    double err = 0;
    try{ err = 0.5*getRescor()*sqrt( pow(sigab/muab,2.));}
    catch(const std::exception&) {err = 0;};
    return err;
  } else {
    double muab = rQABnorm/(double) nevts;
    double muac = rQACnorm/(double) nevts;
    double mubc = rQBCnorm/(double) nevts;
    double varab = rQABnorm2/(double)nevts - pow(muab,2.);
    double varac = rQACnorm2/(double)nevts - pow(muac,2.);
    double varbc = rQBCnorm2/(double)nevts - pow(mubc,2.);
    double sigab = sqrt(varab/(double) nevts);
    double sigac = sqrt(varac/(double) nevts);
    double sigbc = sqrt(varbc/(double) nevts);
    double err = 0;
    if( pow(sigab/muab,2.) + pow(sigac/muac,2.) + pow(sigbc/mubc,2.)>0) err = 0.5*getRescor()*sqrt( pow(sigab/muab,2.) + pow(sigac/muac,2.) + pow(sigbc/mubc,2.));
    return err;
  }
}

double Qvec::getEPVNobs(){
  double na = QnAnorm.real()/(double) ncount;
  return na;
}

double Qvec::getEPVNobsErr(){
  double muna = rQnAnorm/(double) ncount;
  double varna = rQnAnorm2/(double) ncount - pow(muna,2.);
  double signa = sqrt(varna/(double) ncount);
  double err = 0;
  try{ err = 0.5*getEPVNobs()*signa/muna;}
  catch(const std::exception&) {err = 0;};
  return err;
}

double Qvec::getSPVN(){
  if( twosub ) {
    double ab = QAB.real()/wABsum;
    double na = QnA.real()/wnAsum;
    double spres = 0;
    if( ab <=0 ) {
      if(name_.find("qav") !=std::string::npos) std::cout<<"Qvec::getSPVN ("<<name_<<") error ab = "<<ab<<std::endl;
      return 0;
    } else {
      spres = sqrt( ab );
      return na/spres;
    }
  } else {
    double ab = QAB.real()/wABsum;
    double ac = QAC.real()/wACsum;
    double bc = QBC.real()/wBCsum;
    double na = QnA.real()/wnAsum;
    double spres = 0;

    if( ab*ac/bc <=0 ) {
      if(name_.find("qav") !=std::string::npos) std::cout<<"Qvec::getSPVN ("<<name_<<") error spres = "<<ab<<"\t"<<ac<<"\t"<<bc<<std::endl;
      return 0;
    } else {
      spres = sqrt( ab*ac/bc );
      return na/spres;
    }
  }
}
double Qvec::getSPVNErr(){
  if( twosub ) {
    double muab = QAB.real()/wABsum;
    double muna = QnA.real()/wnAsum;
    double varab = rQAB2/wABsum - pow(muab,2.);
    double varna = rQnA2/wnAsum - pow(muna,2.);
    
    double sigab = sqrt(varab/(double)nevts);
    double signa = sqrt(varna/(double)nevts);
    double err = 0;
    try{ err = 0.5*getSPVN()*sqrt( pow(sigab/muab,2.));}
    catch(const std::exception&) {err = 0;};
    return err;
  } else {
    double muab = QAB.real()/wABsum;
    double muac = QAC.real()/wACsum;
    double mubc = QBC.real()/wBCsum;
    double muna = QnA.real()/wnAsum;
    double varab = rQAB2/wABsum - pow(muab,2.);
    double varac = rQAC2/wACsum - pow(muac,2.);
    double varbc = rQBC2/wBCsum - pow(mubc,2.);
    double varna = rQnA2/wnAsum - pow(muna,2.);
    
    double sigab = sqrt(varab/(double)nevts);
    double sigac = sqrt(varac/(double)nevts);
    double sigbc = sqrt(varbc/(double)nevts);
    double signa = sqrt(varna/(double)nevts);
    double err = 0;
    if(pow(sigab/muab,2.) + pow(sigac/muac,2.) + pow(sigbc/mubc,2.) + pow(signa/muna,2.)>0) err = 0.5*getSPVN()*sqrt( pow(sigab/muab,2.) + pow(sigac/muac,2.) + pow(sigbc/mubc,2.) + pow(signa/muna,2.));
    return err;
  }
}

double Qvec::getSPVNnum(){
  double na = QnA.real()/wnAsum;
  return na;
}

double Qvec::getSPVNnumerr(){
  double muna = rQnA/wnAsum;
  double varna = rQnA2/wnAsum - pow(muna,2.);
  double signa = sqrt(varna/(double) ncount);
  double err = 0;
  if(muna>0) { err = getSPVNnum()*signa/muna;}
  return err;
}
double Qvec::getSPVNdenom(){
  double spres = 0;
  if( twosub ) {
    if(wABsum>0) {
      double ab = QAB.real()/wABsum;
      if(ab>0) { spres = sqrt( ab ); 
      } else {
	std::cout<<"getSPVNdenom(2) ("<<name_<<") "<<ab<<std::endl;
      } 
    }
    return spres;
  } else {
    if(wABsum>0 && wACsum>0 && wBCsum>0 ) {
      double ab = QAB.real()/wABsum;
      double ac = QAC.real()/wACsum;
      double bc = QBC.real()/wBCsum;
      if(ab*ac/bc > 0) {
	spres = sqrt( ab*ac/bc ); 
      } else {
	std::cout<<"getSPVNdenom(3) ("<<name_<<") "<<ab*ac/bc<<std::endl;
      }
    }
    return spres;
  }
}
double Qvec::getSPVNdenomerr(){
  //NOTE:  THERE MAY BE AN ERROR IN THIS CALCULATION.  THE WEIGHT 
  //       USED FOR RQXY2 MAY BE INCORRECT.  IF SO, THE FIX IS NEEDED IN
  //       THE BASE REPLAY.  THE WORKAROUND IS TO CALCULATE THE UNCERTAINTIES
  //       BY SUBDIVIDING THE DATA.
  if( twosub ) {
    double muab = rQAB/wABsum;
    double varab = rQAB2/wABsum - pow(muab,2.);
    double sigab = sqrt(varab/(double) nevts);
    double err = 0;
    try{ err = 0.5*getSPVN()*sqrt( pow(sigab/muab,2.));}
    catch(const std::exception&) {err = 0;};
    return err;
  } else {
    double muab = rQAB/wABsum;
    double muac = rQAC/wACsum;
    double mubc = rQBC/wBCsum;
    double varab = rQAB2/wABsum - pow(muab,2.);
    double varac = rQAC2/wACsum - pow(muac,2.);
    double varbc = rQBC2/wBCsum - pow(mubc,2.);
    double sigab = sqrt(varab/(double) nevts);
    double sigac = sqrt(varac/(double) nevts);
    double sigbc = sqrt(varbc/(double) nevts);
    double err = 0;
    try{ err = 0.5*getSPVN()*sqrt( pow(sigab/muab,2.) + pow(sigac/muac,2.) + pow(sigbc/mubc,2.));}
    catch(const std::exception&) {err = 0;};
    return err;
  }
}


double Qvec::getEPVN(){
  double epvn = 0;
  if( twosub ) {
    if(nevts>0 && ncount > 0) {
      double ab = QABnorm.real()/(double) nevts;
      double na = QnAnorm.real()/(double) ncount;
      double res = 0;
      if( ab>0) { 
	res = sqrt( ab ); 
	epvn = na/res;
      } 
    }
    return epvn;
  } else {
    if(nevts>0&&ncount>0) { 
      double ab = QABnorm.real()/(double) nevts;
      double ac = QACnorm.real()/(double) nevts;
      double bc = QBCnorm.real()/(double) nevts;
      double na = QnAnorm.real()/(double) ncount;
      double res = 0;
      if(bc>0 && ab*ac>0) { 
	res = sqrt( ab*ac/bc ); 
	epvn = na/res;
      }
    }
    return epvn;
  }
}

double Qvec::getEPVNErr(){
  if( twosub ) {
    double muab = rQABnorm/(double) nevts;
    double muna = rQnAnorm/(double) ncount;
    double varab = rQABnorm2/(double)nevts - pow(muab,2.);
    double varna = rQnAnorm2/(double)ncount - pow(muna,2.);
    double sigab = sqrt(varab/(double) nevts);
    double signa = sqrt(varna/(double) ncount);
    double err = 0;
    try{ err = 0.5*getEPVN()*sqrt( pow(sigab/muab,2.) + pow(signa/muna,2.));}
    catch(const std::exception&) {err = 0;};
    return err;
  } else {
    double muab = rQABnorm/(double) nevts;
    double muac = rQACnorm/(double) nevts;
    double mubc = rQBCnorm/(double) nevts;
    double muna = rQnAnorm/(double) ncount;
    double varab = rQABnorm2/(double)nevts - pow(muab,2.);
    double varac = rQACnorm2/(double)nevts - pow(muac,2.);
    double varbc = rQBCnorm2/(double)nevts - pow(mubc,2.);
    double varna = rQnAnorm2/(double)ncount - pow(muna,2.);
    double sigab = sqrt(varab/(double) nevts);
    double sigac = sqrt(varac/(double) nevts);
    double sigbc = sqrt(varbc/(double) nevts);
    double signa = sqrt(varna/(double) ncount);
    double err = 0;
    try{ err = 0.5*getEPVN()*sqrt( pow(sigab/muab,2.) + pow(sigac/muac,2.) + pow(sigbc/mubc,2.) + pow(signa/muna,2.));}
    catch(const std::exception&) {err = 0;};
    return err;
  }
}

double Qvec::getAvPt(){
  double ravpt = 0;
  if(ncount> 0) {
    ravpt =  avpt/(double) ncount;
  }
  return ravpt;
}

double Qvec::getAvPtErr(){
  double mupt = avpt/(double)ncount;
  double varpt = avpt2/(double)ncount - pow(mupt,2.);
  double sigpt = sqrt(varpt/(double)ncount);
  double err = 0;
  try{ err = 0.5*getAvPt()*sigpt/mupt;}
  catch(const std::exception&) {err = 0;};

  return err;
}

double Qvec::getAvNoff(){
  return avnoff/nevts;
}
double Qvec::getAvNoffErr(){
  double munoff = avnoff/(double)nevts;
  double varnoff = avnoff2/(double)nevts - pow(munoff,2.);
  double signoff = sqrt(varnoff/(double)nevts);
  double err = 0;
  try{ err = 0.5*getAvNoff()*signoff/munoff;}
  catch(const std::exception&) {err = 0;};

  return err;
  
}
#endif
