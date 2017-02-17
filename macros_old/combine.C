#include "getPtGraph.C"
#include "getIntGraph.C"
FILE * fout;
void DoIt(string reac = "PbPb", string dir="../PbPb_2015/vnanal_2015/results/", string mid = "trackmid2",string scase="2015_N2SUB",string rnge = "25-30"){
  double cent = 0;
  if(rnge=="0-5") cent = 2.5;
  if(rnge=="5-10") cent = 7.5;
  if(rnge=="10-15") cent = 12.5;
  if(rnge=="15-20") cent = 17.5;
  if(rnge=="20-25") cent=22.5;
  if(rnge=="25-30") cent=27.5;
  if(rnge=="30-35") cent=32.5;
  if(rnge=="35-40") cent=37.5;
  if(rnge=="40-45") cent=42.5;
  if(rnge=="45-50") cent=47.5;
  if(rnge=="50-60") cent=55;
  if(rnge=="60-70") cent=65;
  string redstr = scase.substr(scase.find("_")+1,scase.length()-scase.find("_")-1);
  string scalc;
  int vn;
  int epord;
  bool chi = false;
  string scase2;
  int vn2;
  int epord2;
  string mid2;
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
  if(redstr == "chi42SUB") {
    scalc = "#chi_{4} {#Psi_{2}}";
    scase = "2015_N42SUB"; 
    vn = 4;
    epord = 4; 
    scase2 = "2015_D24SUB";
    vn2 = 2;
    epord2 = 2;
    mid2="trackmid2";
    chi = true;

  }
  if(redstr == "chi63SUB") {
    scalc = "#chi_{6} {#Psi_{3}}";
    scase = "2015_N63SUB"; 
    vn = 6;
    epord = 6; 
    scase2 = "2015_D34SUB";
    vn2 = 3;
    epord2 = 3;
    mid2="trackmid3";
    chi = true;

  }
  if(redstr == "chi62SUB") {
    scalc = "#chi_{6} {#Psi_{2}}";
    scase = "2015_N62SUB"; 
    vn = 6;
    epord = 6; 
    scase2 = "2015_D26SUB";
    vn2 = 2;
    epord2 = 2;
    mid2="trackmid2";
    chi = true;

  }
  if(redstr == "chi523SUB" ) {
    scalc = "#chi_{5} {#Psi_{23}}";
    scase = "2015_N523SUB"; 
    vn = 5;
    epord = 2; 
    scase2 = "2015_D2232SUB";
    vn2 = 2;
    epord2 = 2;
    mid2="trackmid2";
    chi = true;

  }
  if(redstr == "chi523ASUB" ) {
    scalc = "#chi_{5} {#Psi_{2_{#alpha}3_{#beta}}}";
    scase = "2015_N523ASUB"; 
    vn = 5;
    epord = 2; 
    scase2 = "2015_D2232ASUB";
    vn2 = 2;
    epord2 = 2;
    mid2="trackmid2";
    chi = true;

  }
  if(redstr == "chi723ASUB" ) {
    scalc = "#chi_{7} {#Psi_{2_{#alpha}3_{#beta}}}";
    scase = "2015_N723ASUB"; 
    vn = 7;
    epord = 2; 
    scase2 = "2015_D2432ASUB";
    vn2 = 2;
    epord2 = 2;
    mid2="trackmid2";
    chi = true;

  }
  if(redstr == "D24SUB") {
    scalc = "v_{2}^{2} {#Psi_{2}}"; 
    vn = 2;
    epord = 2; 
  }
  if(redstr == "D26SUB") {
    scalc = "v_{2}^{3} {#Psi_{2}}"; 
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
  if(redstr == "N62SUB" ) {
    scalc = "v_{6} {#Psi_{2}}";
    vn = 6;
    epord = 6; 
  }
  if(redstr == "N63SUB" ) {
    scalc = "v_{6} {#Psi_{3}}";
    vn = 6;
    epord = 6; 
  }
  if(redstr == "D34SUB" ) {
    scalc = "#sqrt{v_{3}^{2}} {#Psi_{3}}";
    vn = 3;
    epord = 3; 
  }
  if(redstr == "D2232SUB" ) {
    scalc = "#sqrt{v_{2}^{2} v_{3}^{2}} {#Psi_{23}}";
    vn = 2;
    epord = 2; 
  }
  if(redstr == "D2432SUB" ) {
    scalc = "#sqrt{v_{2}^{4} v_{3}^{2}} {#Psi_{23}}";
    vn = 2;
    epord = 2; 
  }
  if(redstr == "D2232ASUB" ) {
    scalc = "#sqrt{v_{2}^{2} v_{3}^{2}} {#Psi_{2_{a}3_{b}}}";
    vn = 2;
    epord = 2; 
  }
  if(redstr == "D2432ASUB" ) {
    scalc = "#sqrt{v_{2}^{4} v_{3}^{2}} {#Psi_{2_{a}3_{b}}}";
    vn = 2;
    epord = 2; 
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
  if(redstr == "N523ASUB" ) {
    scalc = "v_{5} {#Psi_{2_{a}3_{b}}}"; 
    vn = 5;
    epord = 2;
  } 
  if(redstr == "N723SUB" ) {
    scalc = "v_{7} {#Psi_{23}}"; 
    vn = 7;
    epord = 2;
  } 
  if(redstr == "N723ASUB" ) {
    scalc = "v_{7} {#Psi_{2_{a}3_{b}}}"; 
    vn = 7;
    epord = 2;
  } 
  string mname = dir+"results_"+mid+"/"+scase+"2/"+reac+"_"+scase+"2_v"+to_string(vn)+"_HFm"+to_string(vn)+".root";

  std::string urnge = rnge;
  urnge.replace(urnge.find("-"),1, "_");
    
  //TGraphErrors * gpt = getPtGraph(mname,rnge);
  TGraphErrors * gim = getIntGraph(mname,rnge);
  TGraphErrors * chipt=0;
  TGraphErrors * chiim=0;
  string cmname;
  if(chi) {
    cmname = dir+"results_"+mid2+"/"+scase2+"2/"+reac+"_"+scase2+"2_v"+to_string(vn2)+"_HFm"+to_string(vn2)+".root";    
    //chipt = getPtGraph(cmname,rnge);
    chiim = getIntGraph(cmname,rnge);
  }
  cout<<"============ "<<rnge<<" ==============="<<endl;
  cout<<mname<<endl;
  double avn = 0;
  double avn2 = 0;
  double avw = 0;
  double sysav = 0;
  int nsysav = 0;
  for(int i = 0; i<gim->GetN(); i++) {
    cout<<gim->GetX()[i]<<"\t"<<gim->GetY()[i]<<"\t"<<gim->GetEY()[i]<<"\t"<<gim->GetEX()[i]<<endl;
    double val = gim->GetY()[i];
    double err = gim->GetEY()[i];
    double w = 1./pow(err,2.);
    double xval = gim->GetX()[i];
    
    if(fabs(xval)<1.6&& val>0) {
      avn+=val*w;
      avn2+=pow(val*w,2.);
      avw+=w;
      sysav+=gim->GetEX()[i];
      ++nsysav;
    }
  }
  double nmean = avn/avw;
  double nvar = sqrt(1/avw);
  double nsys = sysav/nsysav;
  
  cout<<"mean: "<<avn/avw<<endl;
  cout<<"var:  "<<sqrt(1/avw)<<endl;
  cout<<"sys: "<<sysav/nsysav<<endl;
  cout<<cmname<<endl;
  avn = 0;
  avn2 = 0;
  avw = 0;
  sysav = 0;
  nsysav = 0;
  for(int i = 0; i<chiim->GetN(); i++) {
    cout<<chiim->GetX()[i]<<"\t"<<chiim->GetY()[i]<<"\t"<<chiim->GetEY()[i]<<"\t"<<chiim->GetEX()[i]<<endl;
    double val = chiim->GetY()[i];
    double err = chiim->GetEY()[i];
    double w = 1./pow(err,2.);
    double xval = gim->GetX()[i];
    if(fabs(xval)<1.6 && val>0) {      
      avn+=val*w;
      avn2+=pow(val*w,2.);
      avw+=w;
      sysav+=gim->GetEX()[i];
      ++nsysav;
    }
  }
  double dmean = avn/avw;
  double dvar = sqrt(1/avw);
  double dsys = sysav/nsysav;
  cout<<"mean: "<<avn/avw<<endl;;
  cout<<"var:  "<<sqrt(1/avw)<<endl;
  cout<<"sys: "<<sysav/nsysav<<endl;
  cout<<cmname<<endl;

  double chiv = nmean/dmean;
  double chistat = (nvar/nmean)*chiv;
  double chisys = (nsys/nmean)*chiv;
  cout<<"FINAL: "<<chiv<<"\t"<<chistat<<"\t"<<chisys<<endl;
  fprintf(fout,"%7.1f\t%8.5f\t%8.5f\t%8.5f\n",cent,chiv,chistat,chisys);
  // npts = 0;
  // for(int i = 0; i<im3->GetNbinsX(); i++) {
  //   x[i] = im3->GetBinCenter(i+1);
    
  //   double y1 = im3->GetBinContent(i+1);
  //   double e1 = im3->GetBinError(i+1);
  //   val = y1;
  //   err = e1;
  //   if(chi) {
  //     double y2 = cim3->GetBinContent(i+1);
  //     double e2 = cim3->GetBinError(i+1);
  //     val = y1/y2;
  //     err = val*sqrt(pow(e1/y1,2)+pow(e2/y2,2));
  //   }
  //   y[i] = val;
  //   ey[i] = err; 
  //   ++npts;
  // }
  // TGraphErrors * gim3 = new TGraphErrors(npts,x,y,0,ey);
  // gim3->SetMarkerStyle(20);
  // gim3->SetMarkerColor(kRed);
  // gim3->SetLineColor(kRed);

  // npts = 0;
  // for(int i = 0; i<ip3->GetNbinsX(); i++) {
  //   x[i] = ip3->GetBinCenter(i+1);
    
  //   double y1 = ip3->GetBinContent(i+1);
  //   double e1 = ip3->GetBinError(i+1);
  //   val = y1;
  //   err = e1;
  //   if(chi) {
  //     double y2 = cip3->GetBinContent(i+1);
  //     double e2 = cip3->GetBinError(i+1);
  //     val = y1/y2;
  //     err = val*sqrt(pow(e1/y1,2)+pow(e2/y2,2));
  //   }
  //   y[i] = val;
  //   ey[i] = err; 

  //   ++npts;
  // }
  // TGraphErrors * gip3 = new TGraphErrors(npts,x,y,0,ey);
  // gip3->SetMarkerStyle(24);
  // gip3->SetMarkerColor(kBlue);
  // gip3->SetLineColor(kBlue);



  // string datadir = dir+"data";
  // if((figout = fopen(datadir.data(),"r")) != NULL) {
  //   fclose(figout);
  // } else {
  //   system(Form("mkdir %s",datadir.data()));
  // }

  // string ptdir = dir+"data/pt";
  // if((figout = fopen(ptdir.data(),"r")) != NULL) {
  //   fclose(figout);
  // } else {
  //   system(Form("mkdir %s",datadir.data()));
  // }

  // string etadir = dir+"data/eta";
  // if((figout = fopen(etadir.data(),"r")) != NULL) {
  //   fclose(figout);
  // } else {
  //   system(Form("mkdir %s",etadir.data()));
  // }

  // c->Print(Form("%s/%s_%s.pdf",figdir.data(),redstr.data(),rnge.data()),"pdf");
  // FILE * out ;

  // out = fopen(Form("%s/pt/%s2_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gptm2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptm2->GetX()[i],gptm2->GetY()[i],gptm2->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/pt/%s2_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gptp2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptp2->GetX()[i],gptp2->GetY()[i],gptp2->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/pt/%s2_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gpt2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gpt2->GetX()[i],gpt2->GetY()[i],gpt2->GetEY()[i]);
  // fclose(out);

  // out = fopen(Form("%s/pt/%s3_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gptm3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptm3->GetX()[i],gptm3->GetY()[i],gptm3->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/pt/%s3_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gptp3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gptp3->GetX()[i],gptp3->GetY()[i],gptp3->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/pt/%s3_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gpt3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gpt3->GetX()[i],gpt3->GetY()[i],gpt3->GetEY()[i]);
  // fclose(out);


  // out = fopen(Form("%s/eta/%s2_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gim2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gim2->GetX()[i],gim2->GetY()[i],gim2->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/eta/%s2_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gip2->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip2->GetX()[i],gip2->GetY()[i],gip2->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/eta/%s2_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gip2->GetN(); i++) {
  //   if(gip2->GetX()[i]<0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip2->GetX()[i],gip2->GetY()[i],gip2->GetEY()[i]);
  //   if(gim2->GetX()[i]>0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip2->GetX()[i],gip2->GetY()[i],gip2->GetEY()[i]);
  // }
  // fclose(out);

  // out = fopen(Form("%s/eta/%s3_hfm_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gim3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gim3->GetX()[i],gim3->GetY()[i],gim3->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/eta/%s3_hfp_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gip3->GetN(); i++) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip3->GetX()[i],gip3->GetY()[i],gip3->GetEY()[i]);
  // fclose(out);
  // out = fopen(Form("%s/eta/%s2_%s.dat",datadir.data(),redstr.data(),rnge.data()),"w");
  // for(int i = 0; i<gip3->GetN(); i++) {
  //   if(gip3->GetX()[i]<0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip3->GetX()[i],gip3->GetY()[i],gip3->GetEY()[i]);
  //   if(gim3->GetX()[i]>0) fprintf(out,"%5.2lf\t%9.6lf\t%9.6lf\n",gip3->GetX()[i],gip3->GetY()[i],gip3->GetEY()[i]);
  // }
  // fclose(out);

  return;
}

void combine(string anal="N2"){
  string danal = "2015_"+anal+"SUB";
  string foutn = "data/chi/"+anal+".dat";
  fout = fopen(foutn.data(),"w");
  int epord = 2;
  string mid;
  double scale = 1.; ;
  double lowscale;
  double verylowscale;
  double highscale;
  double ptmax=1.;
  double ptmin = 10.;
  double etamin = 10;
  double etamax=0.4;
  if(anal=="chi42" || anal=="chi523" || anal=="chi523A" || anal=="chi723A" || anal=="chi62"){
    mid="trackmid2";
  }
  if(anal=="chi63" ){
    mid="trackmid3";
  }
  if(anal=="N2") {
    mid="trackmid2";
  }
  if(anal=="N3") {
    mid="trackmid3";
  }
  if(anal=="N4") {
    mid="trackmid4";
  }
  if(anal=="N42"||anal=="D24"||anal=="D26") {
    mid="trackmid2";
  }
  if(anal=="D2232" || anal=="D2232A") {
    mid="trackmid2";
  }
  if(anal=="D2432" || anal=="D2432A") {
    mid="trackmid2";
  }
  if(anal=="N5") {
    mid="trackmid5";
  }
  if(anal=="N6") {
    mid="trackmid6";
  }
  if(anal=="N523" || anal=="N523A") {
    mid="trackmid2";
  }
  if(anal=="N63") {
    mid="trackmid3";
  }
  if(anal=="N62") {
    mid="trackmid2";
  }
  if(anal=="D34") {
    mid="trackmid3";
  }
  if(anal=="N723" || anal=="N723A") {
    mid="trackmid2";
  }

  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"0-5");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"5-10");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"10-15");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"15-20");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"20-25");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"25-30");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"30-35");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"35-40");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"40-45");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"45-50");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"50-60");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/", mid,danal,"60-70");
}
