TGraphErrors * getIntGraph(string infile="../PbPb_2015/vnanal_2015/results/results_trackmid2/2015_N42SUB2/PbPb_2015_N42SUB2_v4_HFm4.root", string rnge = "25-30"){
  std::string urnge = rnge;
  urnge.replace(urnge.find("-"),1, "_");
  string sp = urnge+"/spInt";
  TFile * tf[11];
  TH1D * thm[11];
  TH1D * thp[11];
  TFile * out = new TFile("tmp.root","recreate");
  TDirectory * dir = gDirectory;
  tf[0] = new TFile(infile.data(),"r");
  string savpt = urnge+"/avpt";
  TH2D * avpt = (TH2D *) tf[0]->Get(savpt.data());
  TH1D * avpt1 = avpt->ProjectionX("avpt1",5,8);
  avpt1->Scale(0.25);
  avpt1->SetDirectory(0);
  for(int i = 1; i<=10; i++) {
    string ctmp = infile;
    ctmp.replace(ctmp.find(".root"),5,"_"+to_string(i)+".root");
    tf[i] = new TFile(ctmp.data(),"r");
  }
  for(int i = 0; i<=10; i++) {
    string name1d = rnge+to_string(i);
    thm[i] = (TH1D *) tf[i]->Get(sp.data());
    thm[i]->SetDirectory(0);
    tf[i]->Close();
  }
  infile.replace(infile.find("HFm"),3,"HFp");
  tf[0] = new TFile(infile.data(),"r");
  for(int i = 1; i<=10; i++) {
    string ctmp = infile;
    ctmp.replace(ctmp.find(".root"),5,"_"+to_string(i)+".root");
    tf[i] = new TFile(ctmp.data(),"r");
  }
  for(int i = 0; i<=10; i++) {
    TH2D * tmp = (TH2D *) tf[i]->Get(sp.data());
    string name1d = rnge+"p"+to_string(i);
    thp[i] = (TH1D *) tf[i]->Get(sp.data());
    thp[i]->SetDirectory(0);
    tf[i]->Close();
  }

  int npts = 0;
  double x[40];
  double y[40];
  double ex[40];
  double ey[40];
  double val;
  double err;
  
  out->cd();
  TF1 * fg = new TF1("fg","gaus",-2,2);
  fg->FixParameter(0,0.);
  fg->SetParameter(1,1.);
  TH1D * hsys = (TH1D *) thm[0]->Clone("hsys");
  for(int i = 0; i<hsys->GetNbinsX(); i++) {
    x[i] = hsys->GetBinCenter(i+1);
    TH1D * hg = new TH1D(Form("hg_%d",i),"hg",20,-5,5);
    double avm = 0;
    double avm2 = 0;
    double avp = 0;
    double avp2 = 0;
    for(int j = 1; j<=10; j++) {
      avm+=thm[j]->GetBinContent(i+1);
      avm2+=pow(thm[j]->GetBinContent(i+1),2.);
    }
    for(int j = 1; j<=10; j++) {
      avp+=thp[j]->GetBinContent(i+1);
      avp2+=pow(thp[j]->GetBinContent(i+1),2.);
    }
    double sigm = sqrt( (10*avm2-avm*avm)/(10.*9.));
    avm/=10.;
    double sigp = sqrt( (10*avp2-avp*avp)/(10.*9.));
    avp/=10.;
    for(int j = 1; j<=10; j++) {
      double err = (thm[j]->GetBinContent(i+1)-avm)/sigm;
      hg->Fill(err);
    }
    for(int j = 1; j<=10; j++) {
      double err = (thp[j]->GetBinContent(i+1)-avp)/sigp;
      hg->Fill(err);
    }
    hg->Fit(fg,"NQ");
    ++npts;
    y[i]=(avm*sigm+avp*sigp)/(sigm+sigp);
    ey[i]=sqrt(sigm*sigm+sigp*sigp);
    hsys->SetBinContent(i+1, fabs(avp-avm));
    hsys->SetBinError(i+1,ey[i]);
    //cout<<"..."<<hsys->GetBinCenter(i+1)<<"\t"<<hsys->GetBinContent(i+1)<<"\t"<<hsys->GetBinError(i+1)<<endl;
  }
  
  TF1 * fitsys = new TF1("fitsys","pol2",-2.5,2.5);
  hsys->Fit(fitsys,"QR");
  TGraphErrors * g = new TGraphErrors(hsys->GetNbinsX(),x,y,0,ey);
  for(int i = 0; i<hsys->GetNbinsX(); i++) {
    g->GetEX()[i] = fitsys->Eval(g->GetX()[i]);
  }
  
  return g;
}
