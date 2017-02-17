void DoIt(string reac = "PbPb", string dir="../PbPb_2015/vnanal_2015/results/data/", string danal = "N24",string dnrm="D24",string rnge = "25-30"){
  string scalc;
  int vn;
  int epord;
  if(danal == "N24") {
    scalc = "v_{4} {#Psi_{2}}"; 
    vn = 4;
    epord = 4; 
  }
  string ptname2 = dir+"pt/"+danal+"SUB2_"+rnge+".dat";
  string nrmname2 =  dir+"pt/"+dnrm+"SUB2_"+rnge+".dat";
  TGraphErrors * gpt = new TGraphErrors(ptname2.data(),"%lg %lg %lg");
  gpt->SetMarkerStyle(20);
  gpt->SetMarkerColor(kBlue);
  TGraphErrors * gnrm = new TGraphErrors(nrmname2.data(),"%lg %lg %lg");
  gnrm->SetMarkerStyle(24);
  gnrm->SetMarkerColor(kRed);
  TH1D * h = new TH1D("h","h",100,0,8);
  h->SetMinimum(0);
  h->SetMaximum(0.1);
  TH1D * h2 = new TH1D("h2","h2",100,0,8);
  h2->SetMinimum(0);
  h2->SetMaximum(2);

  TCanvas * c = new TCanvas("c","c",1000,800);
  c->Divide(2,2);
  c->cd(1);
  h->Draw();
  gpt->Draw("p");
  gnrm->Draw("p");
  c->cd(2);
  TGraphErrors * chipt2 = (TGraphErrors *) gpt->Clone("chipt2");
  for(int i = 0; i<chipt2->GetN(); i++) {
    chipt2->GetY()[i] = chipt2->GetY()[i]/gnrm->GetY()[i];
    chipt2->GetEY()[i] = chipt2->GetEY()[i]/gnrm->GetY()[i];
  }
  h2->Draw();
  chipt2->Draw("p");
  return;
}

void nonLinear(string anal="N42", string norm="D24"){
  string danal = anal;
  string dnrm = norm;
  //DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"0-5");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"5-10");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"10-15");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"15-20");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"20-25");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"25-30");
  DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"30-35");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"35-40");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"40-45");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"45-50");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"50-60");
  // DoIt("PbPb", "../PbPb_2015/vnanal_2015/results/data/", danal,dnrm,"60-70");
}
