#include "stdio.h"

void HIN_11_005_Compare_doit(string pre="v4", string dir = "v4anal_2015", string AnalName="N4SUB2", string EP = "HFm4", int mincent=35, int maxcent=40){
  string fpre = pre;
  if(pre=="v42") fpre = "v4_2";
  std::string f005 = "data/HIN-11-005_EP/PtDists/"+fpre+"_"+to_string(mincent)+"_"+to_string(maxcent)+".txt";
  FILE * fin = fopen(f005.data(),"r");
  char buf[80];
  double x005[40];
  double y005[40];
  double stat005[40];
  double sys005[40];
  int n005 = 0;
  while(fgets(buf,80,fin)!=NULL) {
    sscanf(buf,"%lf\t%lf\t%lf\t%lf",&x005[n005],&y005[n005],&stat005[n005],&sys005[n005]);
    ++n005;
  }
  TGraphErrors * g005 = new TGraphErrors(n005,x005,y005,0,stat005);
  g005->SetMarkerStyle(20);
  g005->SetMarkerColor(kGreen);
  g005->SetLineColor(kGreen);
  TCanvas * c = new TCanvas(Form("c_%s_%s_%d_%d",pre.data(),EP.data(),mincent,maxcent),"c",800,600);
  TH1D * h = new TH1D("h","h",100,0,8);
  h->SetMaximum(0.12);
  if(pre=="v42") h->SetMaximum(0.08);
  double scale = 1.;
  if(pre=="v42") scale = 0.08/0.12;
  h->SetMinimum(0.0);
  h->SetXTitle("p_{T}");
  h->SetYTitle(pre.data());
  h->Draw();
  g005->Draw("p");
  string mid = "results_trackmid4";
  if(pre=="v42") mid = "results_trackmid2";
  std::string rfile = "/home/sanders/VNAnal/"+dir+"/v4results/"+mid+"/2015/PbPb_2015_"+AnalName+"_v4_"+EP+".root";
  cout<<rfile<<endl;
  TFile * tf = new TFile(rfile.data(),"read");
  std::string cspdir = to_string(mincent)+"_"+to_string(maxcent)+"/spVn";
  TH2D * spVn = (TH2D *) tf->Get(cspdir.data());
  TH1D * sp = spVn->ProjectionX("sp",5,8);
  sp->Scale(0.25);
  sp->SetMarkerColor(kBlue);
  sp->SetLineColor(kBlue);
  sp->Draw("same");

  std::string cepdir = to_string(mincent)+"_"+to_string(maxcent)+"/epVn";
  TH2D * epVn = (TH2D *) tf->Get(cepdir.data());
  TH1D * ep = epVn->ProjectionX("ep",5,8);
  ep->Scale(0.25);
  ep->SetMarkerColor(kCyan);
  ep->SetLineColor(kCyan);
  ep->Draw("same");
  TLegend * leg = new TLegend(0.4,0.2,0.5,0.4);
  leg->SetTextFont(43);
  leg->SetTextSize(22);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->AddEntry(g005,"EP 2.76 TeV HIN-11-005","lp");
  leg->AddEntry(sp,"SP 5.02 TeV Re-reco PbPb","lp");
  leg->AddEntry(ep,"EP 5.02 TeV Re-reco PbPb","lp");
  leg->Draw();
  TLatex * crange = new TLatex(1,scale*0.1,Form("%d - %d%c",mincent,maxcent,'%'));
  crange->SetTextFont(43);
  crange->SetTextSize(32);
  crange->Draw();
  c->Print(Form("figures/%s_%s_%d_%d.pdf",pre.data(),EP.data(),mincent,maxcent));
  TLatex * epn = new TLatex(1,scale* 0.09,EP.data());
  epn->SetTextFont(43);
  epn->SetTextSize(28);
  epn->Draw();
}

void HIN_11_005_Compare(){
  //HIN_11_005_Compare_doit("v4","v4anal_2015","N4SUB2","HFm4",5,10);
  //HIN_11_005_Compare_doit("v4","v4anal_2015","N4SUB2","HFm4",20,25);
  //HIN_11_005_Compare_doit("v4","v4anal_2015","N4SUB2","HFm4",35,40);
  //HIN_11_005_Compare_doit("v4","v4anal_2015","N4SUB2","HFm4",50,60);
  HIN_11_005_Compare_doit("v42","v42anal","N42SUB2","HFm4",5,10);
  HIN_11_005_Compare_doit("v42","v42anal","N42SUB2","HFm4",20,25);
  HIN_11_005_Compare_doit("v42","v42anal","N42SUB2","HFm4",35,40);
  HIN_11_005_Compare_doit("v42","v42anal","N42SUB2","HFm4",50,60);
}
