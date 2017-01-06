void DoIt(string reac = "PbPb", string dir="../PbPb_2015/vnanal_2015/results/", string mid = "trackmid2",string scase="2015_N2SUB",string rnge = "25-30", double vmin = 0, double vmax=0.5){
  string redstr = scase.substr(scase.find("_")+1,scase.length()-scase.find("_")-1);
  string scalc;
  int vn;
  int epord;
  if(redstr == "N2SUB" ) {
    scalc = "v_{2} {#Psi_{2}}";
    vn = 2;
    epord = 2; 
  }
  if(redstr == "N3SUB" ) {
    scalc = "v_{3} {#Psi_{3}}";
    vn = 3;
    epord = 3; 
  }
  if(redstr == "N4SUB" ) {
    scalc = "v_{4} {#Psi_{4}}"; 
    vn = 4;
    epord = 4; 
  }
  if(redstr == "N42SUB") {
    scalc = "v_{4} {#Psi_{2}}"; 
    vn = 4;
    epord = 4; 
  }
  if(redstr == "D24SUB") {
    scalc = "v_{2}^{2} {#Psi_{2}}"; 
    vn = 2;
    epord = 2; 
  }
  if(redstr == "N4SUB" ) {
    scalc = "v_{4} {#Psi_{4}}"; 
    vn = 4;
    epord = 4; 
  }
  if(redstr == "N5SUB" ) {
    scalc = "v_{5} {#Psi_{5}}";
    vn = 5;
    epord = 5; 
  }
  if(redstr == "N6SUB" ) {
    scalc = "v_{6} {#Psi_{6}}";
    vn = 6;
    epord = 6; 
  }
  if(redstr == "N63SUB" ) {
    scalc = "v_{6} {#Psi_{3}}";
    vn = 6;
    epord = 6; 
  }
  if(redstr == "D34SUB" ) {
    scalc = "v_{3}^{2} {#Psi_{3}}";
    vn = 3;
    epord = 3; 
  }
  if(redstr == "N7SUB" ) {
    scalc = "v_{7} {#Psi_{7}}";
    vn = 7;
    epord = 7; 
  }
  if(redstr == "N523SUB" ) {
    scalc = "v_{5} {#Psi_{23}}"; 
    vn = 5;
    epord = 2;
  } 
  if(redstr == "N723SUB" ) {
    scalc = "v_{7} {#Psi_{23}}"; 
    vn = 7;
    epord = 2;
  } 
  string mname2 = dir+"results_"+mid+"/"+scase+"2/"+reac+"_"+scase+"2_v"+to_string(vn)+"_HFm"+to_string(vn)+".root";
  string pname2 = dir+"results_"+mid+"/"+scase+"2/"+reac+"_"+scase+"2_v"+to_string(vn)+"_HFp"+to_string(vn)+".root";
  cout<<mname2<<endl;
  TFile * tm2 = new TFile(mname2.data(),"read");
  TFile * tp2 = new TFile(pname2.data(),"read");
  string mname3 = dir+"results_"+mid+"/"+scase+"3/"+reac+"_"+scase+"3_v"+to_string(vn)+"_HFm"+to_string(vn)+".root";
  string pname3 = dir+"results_"+mid+"/"+scase+"3/"+reac+"_"+scase+"3_v"+to_string(vn)+"_HFp"+to_string(vn)+".root";
  TFile * tm3 = new TFile(mname3.data(),"read");
  TFile * tp3 = new TFile(pname3.data(),"read");
  std::string urnge = rnge;
  urnge.replace(urnge.find("-"),1, "_");
  string sp = urnge+"/spVn";
  string sint = urnge+"/spInt";
  string ytitle = scalc+" ("+rnge+"\%)";
  string iytitle = scalc+" (0.3 < p_{t}^{"+rnge+"} < 3.0 GeV/c)";
  string savpt = urnge+"/avpt";
  TH2D * avpt = (TH2D *) tm2->Get(savpt.data());
  TH1D * avpt1 = avpt->ProjectionX("avpt1",5,8);
  avpt1->Scale(0.25);
  TH2D * spm2 = (TH2D *) tm2->Get(sp.data());
  TH2D * spp2 = (TH2D *) tp2->Get(sp.data());
  TH1D * ptm2 = (TH1D *) spm2->ProjectionX("ptm2",5,8);
  ptm2->Scale(0.25);
  ptm2->SetMarkerColor(kRed);
  TH1D * ptp2 = (TH1D *) spp2->ProjectionX("ptp2",5,8);
  ptp2->Scale(0.25);

  TH1D * pt2 = (TH1D *) spm2->ProjectionX("pt2a",7,8);
  pt2->Add( (TH1D *) spp2->ProjectionX("pt2b",5,6));
  pt2->Scale(0.25);

  TH1D * hpt = new TH1D("hpt","hpt",100,0,8);
  hpt->SetMinimum(vmin);
  hpt->SetMaximum(vmax);
  hpt->SetXTitle("p_{t}");
  hpt->SetYTitle(ytitle.data());
  TH1D * im2 = (TH1D *) tm2->Get(sint.data());
  im2->SetMarkerColor(kRed);
  im2->SetMinimum(vmin/2.);
  im2->SetMaximum(vmax/4.);
  //im2->SetYTitle(Form("v_{%d} (0.3 < p_{t} < 3.0 GeV/c)",epord));
  im2->SetYTitle(iytitle.data());
  im2->SetXTitle("#eta");
  TH1D * ip2 = (TH1D *) tp2->Get(sint.data());
  ip2->SetMarkerColor(kBlue);
    
  TH2D * spm3 = (TH2D *) tm3->Get(sp.data());
  TH2D * spp3 = (TH2D *) tp3->Get(sp.data());
  TH1D * ptm3 = (TH1D *) spm3->ProjectionX("ptm3",5,8);
  ptm3->SetMarkerColor(kRed);
  TH1D * ptp3 = (TH1D *) spp3->ProjectionX("ptp3",5,8);
  ptp3->SetMarkerColor(kBlue);
  ptm3->Scale(0.25);
  ptm3->SetMinimum(0.0);
  ptm3->SetMaximum(vmax);
  ptm3->SetXTitle("p_{t} (GeV/c)");
  ptm3->SetYTitle(ytitle.data());
  ptp3->Scale(0.25);


  TH1D * pt3 = (TH1D *) spm3->ProjectionX("pt3a",7,8);
  pt3->Add( (TH1D *) spp2->ProjectionX("pt3b",5,6));
  pt3->Scale(0.25);

  TH1D * im3 = (TH1D *) tm3->Get(sint.data());
  im3->SetMarkerColor(kRed);
  im3->SetMinimum(vmin/2.);
  im3->SetMaximum(vmax/4.);
  //im3->SetYTitle(Form("v_{%d} (0.3 < p_{t} < 3.0 GeV/c)",epord));
  im3->SetYTitle(iytitle.data());
  im3->SetXTitle("#eta");
  TH1D * ip3 = (TH1D *) tp3->Get(sint.data());
  ip3->SetMarkerColor(kBlue);

  int npts = 0;
  double x[40];
  double y[40];
  double ey[40];
  for(int i = 0; i<ptm2->GetNbinsX(); i++) {
    x[i] = avpt1->GetBinContent(i+1);
    y[i] = ptm2->GetBinContent(i+1);
    ey[i] = ptm2->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gptm2 = new TGraphErrors(npts,x,y,0,ey);
  gptm2->SetMarkerStyle(20);
  gptm2->SetMarkerColor(kRed);
  gptm2->SetLineColor(kRed);
  
  npts = 0;
  for(int i = 0; i<ptp2->GetNbinsX(); i++) {
    x[i] = avpt1->GetBinContent(i+1);
    y[i] = ptp2->GetBinContent(i+1);
    ey[i] = ptp2->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gptp2 = new TGraphErrors(npts,x,y,0,ey);
  gptp2->SetMarkerStyle(24);
  gptp2->SetMarkerColor(kBlue);
  gptp2->SetLineColor(kBlue);
  
  npts = 0;
  for(int i = 0; i<pt2->GetNbinsX(); i++) {
    x[i] = avpt1->GetBinContent(i+1);
    y[i] = pt2->GetBinContent(i+1);
    ey[i] = pt2->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gpt2 = new TGraphErrors(npts,x,y,0,ey);
  gpt2->SetMarkerStyle(24);
  gpt2->SetMarkerColor(kGreen);
  gpt2->SetLineColor(kGreen);
  
  npts = 0;
  for(int i = 0; i<ptm3->GetNbinsX(); i++) {
    x[i] = avpt1->GetBinContent(i+1);
    y[i] = ptm3->GetBinContent(i+1);
    ey[i] = ptm3->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gptm3 = new TGraphErrors(npts,x,y,0,ey);
  gptm3->SetMarkerStyle(20);
  gptm3->SetMarkerColor(kRed);
  gptm3->SetLineColor(kRed);
  
  npts = 0;
  for(int i = 0; i<ptp3->GetNbinsX(); i++) {
    x[i] = avpt1->GetBinContent(i+1);
    y[i] = ptp3->GetBinContent(i+1);
    ey[i] = ptp3->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gptp3 = new TGraphErrors(npts,x,y,0,ey);
  gptp3->SetMarkerStyle(24);
  gptp3->SetMarkerColor(kBlue);
  gptp3->SetLineColor(kBlue);
  
  npts = 0;
  for(int i = 0; i<pt3->GetNbinsX(); i++) {
    x[i] = avpt1->GetBinContent(i+1);
    y[i] = pt3->GetBinContent(i+1);
    ey[i] = pt3->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gpt3 = new TGraphErrors(npts,x,y,0,ey);
  gpt3->SetMarkerStyle(24);
  gpt3->SetMarkerColor(kGreen);
  gpt3->SetLineColor(kGreen);

  npts = 0;
  for(int i = 0; i<im2->GetNbinsX(); i++) {
    x[i] = im2->GetBinCenter(i+1);
    y[i] = im2->GetBinContent(i+1);
    ey[i] = im2->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gim2 = new TGraphErrors(npts,x,y,0,ey);
  gim2->SetMarkerStyle(20);
  gim2->SetMarkerColor(kRed);
  gim2->SetLineColor(kRed);

  npts = 0;
  for(int i = 0; i<ip2->GetNbinsX(); i++) {
    x[i] = ip2->GetBinCenter(i+1);
    y[i] = ip2->GetBinContent(i+1);
    ey[i] = ip2->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gip2 = new TGraphErrors(npts,x,y,0,ey);
  gip2->SetMarkerStyle(24);
  gip2->SetMarkerColor(kBlue);
  gip2->SetLineColor(kBlue);

  npts = 0;
  for(int i = 0; i<im3->GetNbinsX(); i++) {
    x[i] = im3->GetBinCenter(i+1);
    y[i] = im3->GetBinContent(i+1);
    ey[i] = im3->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gim3 = new TGraphErrors(npts,x,y,0,ey);
  gim3->SetMarkerStyle(20);
  gim3->SetMarkerColor(kRed);
  gim3->SetLineColor(kRed);

  npts = 0;
  for(int i = 0; i<ip3->GetNbinsX(); i++) {
    x[i] = ip3->GetBinCenter(i+1);
    y[i] = ip3->GetBinContent(i+1);
    ey[i] = ip3->GetBinError(i+1); 
    ++npts;
  }
  TGraphErrors * gip3 = new TGraphErrors(npts,x,y,0,ey);
  gip3->SetMarkerStyle(24);
  gip3->SetMarkerColor(kBlue);
  gip3->SetLineColor(kBlue);

  string cname = scase+rnge;
  TCanvas * c = new TCanvas(cname.data(),"c",1200,1000);
  c->Divide(2,2);
  c->cd(1);
  gPad->SetGrid(1,1);
  hpt->Draw();
  gptm2->Draw("p");
  gptp2->Draw("p");
  c->cd(2);
  gPad->SetGrid(1,1);
  hpt->Draw();
  gptm3->Draw("p");
  gptp3->Draw("p");
  c->cd(3);
  gPad->SetGrid(1,1);
  im2->Reset();
  im2->Draw();
  gim2->Draw("p");
  gip2->Draw("p");
  c->cd(4);
  gPad->SetGrid(1,1);
  im3->Reset();
  im3->Draw();
  gim3->Draw("p");
  gip3->Draw("p");

  c->cd(1);
  TPaveText * tt1 = new TPaveText(0.2,0.8,0.7,0.9,"NDC");
  tt1->SetFillColor(kWhite);
  tt1->SetBorderSize(0);
  tt1->SetTextFont(43);
  tt1->SetTextSize(28);
  tt1->AddText(Form("%s2",redstr.data()));
  tt1->AddText("2 subevent");
  tt1->Draw(); 
  c->cd(2);
  TPaveText * tt2 = new TPaveText(0.2,0.8,0.7,0.9,"NDC");
  tt2->SetFillColor(kWhite);
  tt2->SetBorderSize(0);
  tt2->SetTextFont(43);
  tt2->SetTextSize(28);
  tt2->AddText(Form("%s3",redstr.data()));
  tt2->AddText("3 subevent");
  tt2->Draw(); 
  c->cd(3);
  TLegend * leg = new TLegend(0.2,0.2,0.5,0.3);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(24);
  leg->AddEntry(gptm2,"HF-","p");
  leg->AddEntry(gptp2,"HF+","p");
  leg->Draw();
  FILE * figout;
  string figdir  = dir+"figures";
  if((figout = fopen(figdir.data(),"r")) != NULL) {
    fclose(figout);
  } else {
    system(Form("mkdir %s",figdir.data()));
  }

  string datadir = dir+"data";
  if((figout = fopen(datadir.data(),"r")) != NULL) {
    fclose(figout);
  } else {
    system(Form("mkdir %s",datadir.data()));
  }

  string ptdir = dir+"data/pt";
  if((figout = fopen(ptdir.data(),"r")) != NULL) {
    fclose(figout);
  } else {
    system(Form("mkdir %s",datadir.data()));
  }

  string etadir = dir+"data/eta";
  if((figout = fopen(etadir.data(),"r")) != NULL) {
    fclose(figout);
  } else {
    system(Form("mkdir %s",etadir.data()));
  }

  c->Print(Form("%s/%s_%s.pdf",figdir.data(),redstr.data(),rnge.data()),"pdf");
  FILE * out ;

  out = fopen(Form("%s/pt/%s2_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gptm2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptm2->GetX()[i],gptm2->GetY()[i],gptm2->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/pt/%s2_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gptp2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptp2->GetX()[i],gptp2->GetY()[i],gptp2->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/pt/%s2_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gpt2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gpt2->GetX()[i],gpt2->GetY()[i],gpt2->GetEY()[i]);
  fclose(out);

  out = fopen(Form("%s/pt/%s3_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gptm3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptm3->GetX()[i],gptm3->GetY()[i],gptm3->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/pt/%s3_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gptp3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptp3->GetX()[i],gptp3->GetY()[i],gptp3->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/pt/%s3_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gpt3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gpt3->GetX()[i],gpt3->GetY()[i],gpt3->GetEY()[i]);
  fclose(out);


  out = fopen(Form("%s/eta/%s2_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gim2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gim2->GetX()[i],gim2->GetY()[i],gim2->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/eta/%s2_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gip2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip2->GetX()[i],gip2->GetY()[i],gip2->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/eta/%s2_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gip2->GetN(); i++) {
    if(gip2->GetX()[i]<0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip2->GetX()[i],gip2->GetY()[i],gip2->GetEY()[i]);
    if(gim2->GetX()[i]>0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip2->GetX()[i],gip2->GetY()[i],gip2->GetEY()[i]);
  }
  fclose(out);

  out = fopen(Form("%s/eta/%s3_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gim3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gim3->GetX()[i],gim3->GetY()[i],gim3->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/eta/%s3_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gip3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip3->GetX()[i],gip3->GetY()[i],gip3->GetEY()[i]);
  fclose(out);
  out = fopen(Form("%s/eta/%s2_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  for(int i = 0; i<gip3->GetN(); i++) {
    if(gip3->GetX()[i]<0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip3->GetX()[i],gip3->GetY()[i],gip3->GetEY()[i]);
    if(gim3->GetX()[i]>0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip3->GetX()[i],gip3->GetY()[i],gip3->GetEY()[i]);
  }
  fclose(out);

  return;
}

void showCase(string anal="N2"){
  string danal = "2015_"+anal+"SUB";
  int epord = 2;
  string mid;
  double scale ;
  double lowscale;
  double verylowscale;
  double highscale;
  double vmin = 0;
  if(anal=="N2") {
    mid="trackmid2";
    scale = 0.6;
    lowscale = 0.4;
    verylowscale = 0.2;
    highscale = scale;
  }
  if(anal=="N3") {
    mid="trackmid3";
    scale = 0.2;
    lowscale = 0.2;
    verylowscale = 0.2;
    highscale = scale;
  }
  if(anal=="N4") {
    mid="trackmid4";
    scale = 0.2;
    lowscale = 0.2;
    verylowscale = 0.2;
    highscale = scale;
  }
  if(anal=="N42"||anal=="D24") {
    mid="trackmid2";
    scale = 0.1;
    lowscale = 0.06;
    verylowscale = 0.04;
    highscale = scale;
  }
  if(anal=="N5") {
    mid="trackmid5";
    scale = 0.1;
    lowscale = 0.1;
    verylowscale = 0.1;
    veryhighscale=0.2;
    highscale = scale;
  }
  if(anal=="N6") {
    mid="trackmid6";
    scale = 0.08;
    lowscale = 0.08;
    verylowscale = 0.08;
    veryhighscale=0.2;
    highscale = scale;
  }
  if(anal=="N523") {
    mid="trackmid2";
    scale = 0.1;
    lowscale = 0.06;
    verylowscale = 0.04;
    vmin = -0.02;
    highscale = scale;
  }
  if(anal=="N63") {
    mid="trackmid3";
    scale = 0.1;
    lowscale = 0.06;
    verylowscale = 0.04;
    vmin = -0.02;
    highscale = scale;
  }
  if(anal=="D34") {
    mid="trackmid3";
    scale = 0.1;
    lowscale = 0.06;
    verylowscale = 0.04;
    vmin = -0.02;
    highscale = scale;
  }
  if(anal=="N723") {
    mid="trackmid2";
    scale = 0.04;
    lowscale = 0.02;
    verylowscale = 0.02;
    vmin = -0.02;
    highscale = scale;
  }

  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"0-5", vmin,verylowscale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"5-10", vmin, verylowscale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"10-15", vmin, lowscale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"15-20", vmin, lowscale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"20-25", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"25-30", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"30-35", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"35-40", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"40-45", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"45-50", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"50-60", vmin, scale);
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"60-70", vmin, highscale);
}
