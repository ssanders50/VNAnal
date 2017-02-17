string datadir = "../PbPb_2015/vnanal_2015/results/data/eta/";
string theorydir = "../PbPb_2015/vnanal_2015/results/theory/";
string chidir = "data/chi/";
//string datadir = "SP_results/data/eta/";
//string theorydir = "SP_results/theory/";
TGraph * theory(string chi = "chi4", string calc="0.08 Tfo=150MeV", int color = kGreen);
TGraph * ampt(string chi);
TGraph * amptMin(string chi);
TGraph * amptMax(string chi);
TGraph * amptShade(string chi);

int centbins[]{0,5,10,15,20,25,30,35,40,45,50,60,70};
TGraphErrors * chistat(string chi, int style=20, int color=kBlue){
  string cname = chidir+chi+".dat";
  TGraphErrors * g = new TGraphErrors(cname.data(),"%lf %lf %lf %lf");
  int nbins = g->GetN();
  if(chi=="chi723A") nbins--;
  double x[20];
  double y[20];
  double ey[20];
  for(int i = 0; i<nbins; i++) {
    x[i] = g->GetX()[i];
    y[i] = g->GetY()[i];
    ey[i] = g->GetEX()[i];
  };
  TGraphErrors * gret = new TGraphErrors(nbins,x,y,0,ey);
  gret->SetMarkerStyle(style);
  gret->SetMarkerColor(color);
  gret->SetLineColor(color);
  gret->SetLineWidth(2);
  return gret;
};
TGraphErrors * chisys(string chi){
  string cname = chidir+chi+".dat";
  TGraphErrors * g = new TGraphErrors(cname.data(),"%lf %lf %lf %lf");
  int nbins = g->GetN();
  if(chi=="chi723A") nbins--;
  double x[20];
  double y[20];
  double ex[20];
  double ey[20];
  for(int i = 0; i<nbins; i++) {
    x[i] = g->GetX()[i];
    y[i] = g->GetY()[i];
    ex[i] = 1.5;
    ey[i] = g->GetEY()[i];
  };
  TGraphErrors * gret = new TGraphErrors(nbins,x,y,ex,ey);

  gret->SetFillColorAlpha(12,0.30);
  return gret;
};
void chi(){
  TCanvas * c = new TCanvas("c","c",1250,330);
  c->Divide(5,1,0,0);
  c->cd(1);
  TH1D * h422 = new TH1D("h","h",100,0,79.99);
  h422->SetMinimum(-1);
  h422->SetMaximum(3.999);
  h422->GetXaxis()->SetTitleFont(43);
  h422->GetXaxis()->SetTitleSize(22);
  h422->GetYaxis()->SetTitleFont(43);
  h422->GetYaxis()->SetTitleSize(22);
  h422->GetXaxis()->SetLabelFont(43);
  h422->GetXaxis()->SetLabelSize(18);
  h422->GetYaxis()->SetLabelFont(43);
  h422->GetYaxis()->SetLabelSize(18);
  h422->GetYaxis()->SetTitleOffset(1.2);
  h422->GetYaxis()->CenterTitle();
  h422->GetXaxis()->CenterTitle();
  h422->SetXTitle("Centrality (\%)");
  h422->SetYTitle("#chi");
  h422->Draw();
  chisys("chi42")->Draw("[] 2");
  chistat("chi42")->Draw("p");
  TGraph * ghydro = theory("chi4","0.08 Tfo=150MeV",kGreen);
  TGraph * gAMPT = ampt("chi4");
  ghydro->Draw("l");
  amptMin("chi4")->Draw("l");
  amptMax("chi4")->Draw("l");
  amptShade("chi4")->Draw("f");
  gAMPT->Draw("p");
  TLegend * leg = new TLegend(0.2,0.6,0.85,0.8);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(22);
  leg->AddEntry(ghydro,"hydro #eta/s = 0.08","l");
  leg->AddEntry(gAMPT,"AMPT","lp");
  leg->Draw();
  TLatex * reac = new TLatex(5,3.2,"PbPb #sqrt{s_{NN}}=5.02 TeV");
  reac->SetTextFont(43);
  reac->SetTextSize(24);
  reac->Draw();
  TLatex * l422 = new TLatex(10,-0.5,"#chi_{422}");
  l422->SetTextFont(43);
  l422->SetTextSize(32);
  l422->Draw();
  c->cd(2);
  TH1D * h523 = new TH1D("h523","h",100,0.1,79.99);
  h523->SetMinimum(-1);
  h523->SetMaximum(4);
  h523->GetXaxis()->SetTitleFont(43);
  h523->GetXaxis()->SetTitleSize(22);
  h523->GetYaxis()->SetTitleFont(43);
  h523->GetYaxis()->SetTitleSize(22);
  h523->SetXTitle("Centrality (\%)");
  h523->SetYTitle("#chi_{5}");
  h523->GetXaxis()->SetLabelFont(43);
  h523->GetXaxis()->SetLabelSize(18);
  h523->GetYaxis()->SetLabelFont(43);
  h523->GetYaxis()->SetLabelSize(18);
  h523->GetYaxis()->SetTitleOffset(1.0);
  h523->GetYaxis()->CenterTitle();
  h523->GetXaxis()->CenterTitle();
  h523->Draw();
  chisys("chi523A")->Draw("[] 2");
  chistat("chi523A")->Draw("p");
  theory("chi5","0.08 Tfo=150MeV",kGreen)->Draw("l");
  amptMin("chi5")->Draw("l");
  amptMax("chi5")->Draw("l");
  amptShade("chi5")->Draw("f");
  ampt("chi5")->Draw("p");
  TLatex * l523 = new TLatex(10,-0.5,"#chi_{523}");
  l523->SetTextFont(43);
  l523->SetTextSize(32);
  l523->Draw();
  c->cd(3);
  TH1D * h6222 = new TH1D("h6222","h",100,0.01,79.99);
  h6222->SetMinimum(-1);
  h6222->SetMaximum(4);
  h6222->GetXaxis()->SetTitleFont(43);
  h6222->GetXaxis()->SetTitleSize(22);
  h6222->GetYaxis()->SetTitleFont(43);
  h6222->GetYaxis()->SetTitleSize(22);
  h6222->SetXTitle("Centrality (\%)");
  h6222->SetYTitle("#chi_{62}");
  h6222->GetXaxis()->SetLabelFont(43);
  h6222->GetXaxis()->SetLabelSize(18);
  h6222->GetYaxis()->SetLabelFont(43);
  h6222->GetYaxis()->SetLabelSize(18);
  h6222->GetYaxis()->SetTitleOffset(1.0);
  h6222->GetYaxis()->CenterTitle();
  h6222->GetXaxis()->CenterTitle();
  h6222->Draw();
  chisys("chi62")->Draw("[] 2");
  chistat("chi62")->Draw("p");
  theory("chi62","0.08 Tfo=150MeV",kGreen)->Draw("l");
  amptMin("chi62")->Draw("l");
  amptMax("chi62")->Draw("l");
  amptShade("chi62")->Draw("f");
  ampt("chi62")->Draw("p");
  TLatex * l6222 = new TLatex(10,-0.5,"#chi_{6222}");
  l6222->SetTextFont(43);
  l6222->SetTextSize(32);
  l6222->Draw();
  c->cd(4);
  TH1D * h633 = new TH1D("h633","h",100,0.01,79.99);
  h633->SetMinimum(-1);
  h633->SetMaximum(4);
  h633->GetXaxis()->SetTitleFont(43);
  h633->GetXaxis()->SetTitleSize(22);
  h633->GetYaxis()->SetTitleFont(43);
  h633->GetYaxis()->SetTitleSize(22);
  h633->SetXTitle("Centrality (\%)");
  h633->SetYTitle("#chi_{63}");
  h633->GetXaxis()->SetLabelFont(43);
  h633->GetXaxis()->SetLabelSize(18);
  h633->GetYaxis()->SetLabelFont(43);
  h633->GetYaxis()->SetLabelSize(18);
  h633->GetYaxis()->SetTitleOffset(1.0);
  h633->GetYaxis()->CenterTitle();
  h633->GetXaxis()->CenterTitle();
  h633->Draw();
  chisys("chi63")->Draw("[] 2");
  chistat("chi63")->Draw("p");
  theory("chi63","0.08 Tfo=150MeV",kGreen)->Draw("l");
  amptMin("chi63")->Draw("l");
  amptMax("chi63")->Draw("l");
  amptShade("chi63")->Draw("f");
  ampt("chi63")->Draw("p");
  TLatex * l633 = new TLatex(10,-0.5,"#chi_{633}");
  l633->SetTextFont(43);
  l633->SetTextSize(32);
  l633->Draw();
  c->cd(5);
  TH1D * h723 = new TH1D("h","h",100,0.01,79.99);
  h723->SetMinimum(-1);
  h723->SetMaximum(4);
  h723->GetXaxis()->SetTitleFont(43);
  h723->GetXaxis()->SetTitleSize(22);
  h723->GetYaxis()->SetTitleFont(43);
  h723->GetYaxis()->SetTitleSize(22);
  h723->SetXTitle("Centrality (\%)");
  h723->SetYTitle("#chi_{7}");
  h723->GetXaxis()->SetLabelFont(43);
  h723->GetXaxis()->SetLabelSize(18);
  h723->GetYaxis()->SetLabelFont(43);
  h723->GetYaxis()->SetLabelSize(18);
  h723->GetYaxis()->SetTitleOffset(1.0);
  h723->GetYaxis()->CenterTitle();
  h723->GetXaxis()->CenterTitle();
  h723->Draw();
  chisys("chi723A")->Draw("[] 2");
  chistat("chi723A")->Draw("p");
  theory("chi7","0.08 Tfo=150MeV",kGreen)->Draw("l");
  amptMin("chi7")->Draw("l");
  amptMax("chi7")->Draw("l");
  amptShade("chi7")->Draw("f");
  ampt("chi7")->Draw("p");
  TLatex * l7223 = new TLatex(10,-0.5,"#chi_{7223}");
  l7223->SetTextFont(43);
  l7223->SetTextSize(32);
  l7223->Draw();
  c->Print("resp.pdf","pdf");
  
}


TGraph * theory(string chi, string calc, int color){
  string fname = theorydir+"CMS-"+chi+".dat";
  cout<<fname<<endl;
  FILE * fin=fopen(fname.data(),"r");
  char buf[100];
  double x[20];
  double y[20];
  int npnt = 0;
  while(fgets(buf,100,fin)!=NULL){
    std::string chk = buf;
    if(chk.find(calc)==std::string::npos) continue;
    fgets(buf,100,fin);
    bool getnext = true;
    while(getnext){
      string nchk = fgets(buf,100,fin);
      if(nchk.length()==1) {
	getnext=false;
	continue;
      }
      double fx,fy;
      sscanf(buf,"%lf %lf",&fx,&fy);
      x[npnt] = fx;
      y[npnt] = fy;
      ++npnt;
    }
    break;
  }
  TGraph *g = new TGraph(npnt,x,y);
  g->SetLineColor(kRed);
  g->SetLineStyle(2);
  g->SetLineWidth(2);
  g->SetName(chi.data());
  return g;
}


TGraph * ampt(string chi){
  string calc="AMPT";
  string fname = theorydir+"CMS-"+chi+".dat";
  cout<<fname<<endl;
  FILE * fin=fopen(fname.data(),"r");
  char buf[100];
  double x[20];
  double y[20];
  int npnt = 0;
  while(fgets(buf,100,fin)!=NULL){
    std::string chk = buf;
    if(chk.find(calc)==std::string::npos) continue;
    fgets(buf,100,fin);
    bool getnext = true;
    while(getnext){
      string nchk = fgets(buf,100,fin);
      if(nchk.length()==1) {
	getnext=false;
	continue;
      }
      double fx,fy,fmin,fmax;
      sscanf(buf,"%lf %lf %lf %lf",&fx,&fy,&fmin,&fmax);
      x[npnt] = fx;
      y[npnt] = fy;
      ++npnt;
    }
    break;
  }
  TGraph *g = new TGraph(npnt,x,y);
  g->SetLineColor(kGreen);
  g->SetMarkerColor(kGreen);
  g->SetMarkerSize(0.5);
  g->SetLineWidth(2);
  g->SetName(chi.data());
  fclose(fin);
  return g;
}


TGraph * amptMin(string chi){
  string calc="AMPT";
  string fname = theorydir+"CMS-"+chi+".dat";
  FILE * fin=fopen(fname.data(),"r");
  char buf[100];
  double x[20];
  double y[20];
  int npnt = 0;
  while(fgets(buf,100,fin)!=NULL){
    std::string chk = buf;
    if(chk.find(calc)==std::string::npos) continue;
    fgets(buf,100,fin);
    bool getnext = true;
    while(getnext){
      string nchk = fgets(buf,100,fin);
      if(nchk.length()==1) {
	getnext=false;
	continue;
      }
      double fx,fy,fmin,fmax;
      sscanf(buf,"%lf %lf %lf %lf",&fx,&fy,&fmin,&fmax);
      x[npnt] = fx;
      y[npnt] = fmin;
      ++npnt;
    }
    break;
  }
  TGraph *g = new TGraph(npnt,x,y);
  g->SetLineColor(kGreen);
  g->SetLineWidth(1);
  return g;
}

TGraph * amptMax(string chi){
  string calc="AMPT";
  string fname = theorydir+"CMS-"+chi+".dat";
  FILE * fin=fopen(fname.data(),"r");
  char buf[100];
  double x[20];
  double y[20];
  int npnt = 0;
  while(fgets(buf,100,fin)!=NULL){
    std::string chk = buf;
    if(chk.find(calc)==std::string::npos) continue;
    fgets(buf,100,fin);
    bool getnext = true;
    while(getnext){
      string nchk = fgets(buf,100,fin);
      if(nchk.length()==1) {
	getnext=false;
	continue;
      }
      double fx,fy,fmin,fmax;
      sscanf(buf,"%lf %lf %lf %lf",&fx,&fy,&fmin,&fmax);
      x[npnt] = fx;
      y[npnt] = fmax;
      ++npnt;
    }
    break;
  }
  TGraph *g = new TGraph(npnt,x,y);
  g->SetLineColor(kGreen);
  g->SetLineWidth(1);
  return g;
}

TGraph * amptShade(string chi){
  string calc="AMPT";
  string fname = theorydir+"CMS-"+chi+".dat";
  FILE * fin=fopen(fname.data(),"r");
  char buf[100];
  double x[20];
  double y[20];
  double ymax[20];
  double ymin[20];
  int npnt = 0;
  while(fgets(buf,100,fin)!=NULL){
    std::string chk = buf;
    if(chk.find(calc)==std::string::npos) continue;
    fgets(buf,100,fin);
    bool getnext = true;
    while(getnext){
      string nchk = fgets(buf,100,fin);
      if(nchk.length()==1) {
	getnext=false;
	continue;
      }
      double fx,fy,fmin,fmax;
      sscanf(buf,"%lf %lf %lf %lf",&fx,&fy,&fmin,&fmax);
      x[npnt] = fx;
      y[npnt] = fy;
      ymin[npnt] = fmin;
      ymax[npnt] = fmax;
      ++npnt;
    }
    break;
  }
  TGraph *g = new TGraph(2*npnt);
  for(int i =0; i<npnt; i++) {
    g->SetPoint(i,x[i],ymax[i]);
    g->SetPoint(npnt+i,x[npnt-i-1],ymin[npnt-i-1]);
  }
  g->SetLineColor(kGreen);
  g->SetFillStyle(3013);
  g->SetFillColor(kGreen);
  return g;
}
