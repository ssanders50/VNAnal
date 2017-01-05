TH1D * trkbins;
TH1D * centbins;
GenSupport * genm[30][20];
GenSupport * genp[30][20];
TH2D * qxav[30];
TH2D * qyav[30];


void Setup(string analtype, GetEventInfo * info){
  if(info->isPbPb2015pixel()) trkoff = false;
  cout<<"setup with epord_: "<<epord_<<" ANAL: "<<AnalNames[ANAL].data()<<endl;
  int nbins = 0;
  if(trkoff) {
    nbins = ntrkbins;
    cout<<"Use NtrkOff binning"<<endl;
  } else {
    nbins = ncentbins;
    cout<<"Use Centrality binning"<<endl;
  }
  if(RECENTERTRACKS) {
    cout<<"Track recentering is ON"<<endl;
  } else {
    cout<<"Track recentering is OFF"<<endl;
  }

  trkbins = new TH1D("trkbins","trkbins",ntrkbins,trkBins);
  trkbins->Sumw2();
  centbins = new TH1D("centbins","centbins",ncentbins,centBins);
  centbins->Sumw2();
  trkbins->SetDirectory(0);
  centbins->SetDirectory(0);
  int ep2a = epa;
  int ep2b = epb;
  int ep2c = epc;
  for(int i = 0; i<nbins; i++) {
      qxav[i]=0;
      qyav[i]=0;
  }
  int minbin = 0;
  int maxbin = 0;
  for(int i = 0; i<nbins; i++) {
    if(trkoff){ 
      minbin = trkBins[i];
      maxbin = trkBins[i+1];
    } else {
      minbin = centBins[i];
      maxbin = centBins[i+1];
    }
    genp[i][0] = new GenSupport(epord_, info, ep2a, ep2b, ep2c, ep3a, ep3b, ep3c, minbin, maxbin,0);
    genm[i][0] = new GenSupport(epord_, info, ep2b, ep2a, ep2c, ep3b, ep3a, ep3c, minbin, maxbin,0);
    qxav[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav%d",i));
    qyav[i] = (TH2D *) info->getTemplate()->Clone(Form("qyav%d",i));
    qxav[i]->Reset();
    qyav[i]->Reset();
    for(int j = 1; j<=10; j++) {
      genp[i][j] = new GenSupport(epord_, info, ep2a, ep2b, ep2c, ep3a, ep3b, ep3c, minbin, maxbin,j);
      genm[i][j] = new GenSupport(epord_, info, ep2b, ep2a, ep2c, ep3b, ep3a, ep3c, minbin, maxbin,j);
    }
  }
  return;
}

void OutputResults(TString reac, TString midn, TString trig){
  cout<<"OutputResults"<<endl;
  cout<<"   reac:     "<<reac.Data()<<endl;
  cout<<"   midn:     "<<midn.Data()<<endl;
  cout<<"   trig:     "<<trig.Data()<<endl;
  cout<<"   ANAL:     "<<ANAL<<endl;
  cout<<"   AnalName: "<<AnalNames[ANAL].data()<<endl;
  cout<<"   epord_:   "<<epord_<<endl;
  string strig="dummy";
  int nbins = 0;
  if(trkoff) {
    nbins = ntrkbins;
    cout<<"Use NtrkOff binning"<<endl;
  } else {
    nbins = ncentbins;
    cout<<"Use Centrality binning"<<endl;
  }

  if(trig.Contains("2011ppreco")) strig="2011ppreco";
  if(trig.Contains("2015")) strig="2015";
  if(trig.Contains("MB")) strig="MB";
  if(trig.Contains("HM100")) strig="HM100";
  if(trig.Contains("HM130")) strig="HM130";
  if(trig.Contains("HM160")) strig="HM160";
  if(trig.Contains("HM190")) strig="HM190";
  if(trig.Contains("HM220")) strig="HM220";
  if(trig.Contains("HFp")) strig+="_pSide";
  if(trig.Contains("HFm")) strig+="_PbSide";
  strig+="_";
  strig+=AnalNames[ANAL];
  if(!fopen(Form("results/results_%s",midn.Data()),"r")) system(Form("mkdir results/results_%s",midn.Data()));
  if(!fopen(Form("results/results_%s/%s",midn.Data(),strig.data()),"r")) system(Form("mkdir results/results_%s/%s",midn.Data(),strig.data()));
  TFile * savep[11];
  savep[0] = new TFile  (Form("results/results_%s/%s/%s_%s_v%d_HFp%d.root",midn.Data(),strig.data(),reac.Data(),strig.data(),epord_,epord_),"recreate");
  for(int j = 1; j<=10; j++) { 
    savep[j] = new TFile(Form("results/results_%s/%s/%s_%s_v%d_HFp%d_%d.root", midn.Data(),strig.data(),reac.Data(),strig.data(),epord_,epord_,j),"recreate"); 
  } 
  cout<<"Output results to: "<<savep[0]->GetName()<<endl;
  trkbins->Write();
  centbins->Write();
  for(int i = 0; i<nbins; i++) {
    if(trkoff) {
      cout<<"****** genp track"<<trkBins[i]<<"\t"<<trkBins[i+1]<<endl;
    } else {
      cout<<"****** genp cent"<<centBins[i]<<"\t"<<centBins[i+1]<<endl;
    }
    for(int j = 0; j<=10; j++) {
      if(genp[i][j]->getCount()>10) {
	savep[j]->cd();
	savep[j]->mkdir(Form("%d_%d",genp[i][j]->getMinTrk(),genp[i][j]->getMaxTrk()))->cd();
	if(genp[i][j]->getRescor()) genp[i][j]->getRescor()->Write("rescor");
	if(genp[i][j]->getRescorAv()) genp[i][j]->getRescorAv()->Write("rescorav");
	if(genp[i][j]->getEPVNobs()) genp[i][j]->getEPVNobs()->Write("epVnObs");
	if(genp[i][j]->getEPVN()) genp[i][j]->getEPVN()->Write("epVn");
		if(genp[i][j]->getSPVN()) genp[i][j]->getSPVN()->Write("spVn");
	if(genp[i][j]->getQAB()) genp[i][j]->getQAB()->Write("QAB");
	if(genp[i][j]->getQAC()) genp[i][j]->getQAC()->Write("QAC");
	if(genp[i][j]->getQBC()) genp[i][j]->getQBC()->Write("QBC");
	if(genp[i][j]->getQnA()) genp[i][j]->getQnA()->Write("QnA");
	if(genp[i][j]->getSpecRaw()) genp[i][j]->getSpecRaw()->Write("specRaw");
	if(genp[i][j]->getSpecEffCorrected()) genp[i][j]->getSpecEffCorrected()->Write("specEffCorrected");
	if(genp[i][j]->getEPInt(0.3, 3.0)) genp[i][j]->getEPInt(0.3, 3.0)->Write("epInt");
	if(genp[i][j]->getSPInt(0.3, 3.0)) genp[i][j]->getSPInt(0.3, 3.0)->Write("spInt");
	if(genp[i][j]->getAvPt()) genp[i][j]->getAvPt()->Write("avpt");
	if(genp[i][j]->getAvNoff()) genp[i][j]->getAvNoff()->Write("avnoff");
	if(genp[i][j]->getSPEPRatio()) genp[i][j]->getSPEPRatio()->Write("spepratio");
	if(genp[i][j]->getSPEPIntRatio(0.3, 3.0)) genp[i][j]->getSPEPIntRatio(0.3, 3.0)->Write("spepintratio");
	if(genp[i][j]->getQnFraction()) genp[i][j]->getQnFraction()->Write("QnFraction");
	if(j==0) {
	  if(qxav[i]) qxav[i]->Write("qxav");
	  if(qyav[i]) qyav[i]->Write("qyav");
	}
	TDirectory * vn1d = gDirectory->mkdir("vnObs");
	TDirectory * vn2d = gDirectory->mkdir("vnObs2D");
	vn1d->cd();
	for(int k = 0; k<16; k++) {
	  genp[i][j]->getVnObs(k)->Write(Form("vnObs_%d",k));
	}
	vn2d->cd();
	for(int k = 0; k<16; k++) {
	  genp[i][j]->getVnObs2D(k)->Write(Form("vnObs2D_%d",k));
	}

      }
    }
  }
  TFile * savem[11];
  savem[0] =    new TFile(Form("results/results_%s/%s/%s_%s_v%d_HFm%d.root", midn.Data(),strig.data(),reac.Data(),strig.data(),epord_,epord_),"recreate");
  for(int j = 1; j<=10; j++) { 
    savem[j] =  new TFile(Form("results/results_%s/%s/%s_%s_v%d_HFm%d_%d.root",midn.Data(),strig.data(),reac.Data(),strig.data(),epord_,epord_,j),"recreate"); 
  } 
  
  //cout<<"Output results to: "<<savem[0]->GetName()<<endl;
  trkbins->Write();
  centbins->Write();
  for(int i = 0; i<nbins; i++) {
    if(trkoff) {
      cout<<"****** genm track"<<trkBins[i]<<"\t"<<trkBins[i+1]<<endl;
    } else {
      cout<<"****** genm cent"<<centBins[i]<<"\t"<<centBins[i+1]<<endl;
    }
    for(int j = 0; j<=10; j++) {
      if(genm[i][j]->getCount()>10) {
	savem[j]->cd();
	savem[j]->mkdir(Form("%d_%d",genm[i][j]->getMinTrk(),genm[i][j]->getMaxTrk()))->cd();
	if(genm[i][j]->getRescor()) genm[i][j]->getRescor()->Write("rescor");
	if(genm[i][j]->getRescorAv()) genm[i][j]->getRescorAv()->Write("rescorav");
	if(genm[i][j]->getEPVNobs()) genm[i][j]->getEPVNobs()->Write("epVnObs");
	if(genm[i][j]->getEPVN()) genm[i][j]->getEPVN()->Write("epVn");
	if(genm[i][j]->getSPVN()) genm[i][j]->getSPVN()->Write("spVn");
	if(genm[i][j]->getQAB()) genm[i][j]->getQAB()->Write("QAB");
	if(genm[i][j]->getQAC()) genm[i][j]->getQAC()->Write("QAC");
	if(genm[i][j]->getQBC()) genm[i][j]->getQBC()->Write("QBC");
	if(genm[i][j]->getQnA()) genm[i][j]->getQnA()->Write("QnA");
	if(genm[i][j]->getSpecRaw()) genm[i][j]->getSpecRaw()->Write("specRaw");
	if(genm[i][j]->getSpecEffCorrected()) genm[i][j]->getSpecEffCorrected()->Write("specEffCorrected");
	if(genm[i][j]->getEPInt(0.3, 3.0)) genm[i][j]->getEPInt(0.3, 3.0)->Write("epInt");
	if(genm[i][j]->getSPInt(0.3, 3.0)) genm[i][j]->getSPInt(0.3, 3.0)->Write("spInt");
	if(genm[i][j]->getAvPt()) genm[i][j]->getAvPt()->Write("avpt");
	if(genm[i][j]->getAvNoff()) genm[i][j]->getAvNoff()->Write("avnoff");
	if(genm[i][j]->getSPEPRatio()) genm[i][j]->getSPEPRatio()->Write("spepratio");
	if(genm[i][j]->getSPEPIntRatio(0.3, 3.0)) genm[i][j]->getSPEPIntRatio(0.3, 3.0)->Write("spepintratio");
	if(genm[i][j]->getQnFraction()) genm[i][j]->getQnFraction()->Write("QnFraction");
	if(j==0) {
	  if(qxav[i]) qxav[i]->Write("qxav");
	  if(qyav[i]) qyav[i]->Write("qyav");
	}
      	TDirectory * vn1d = gDirectory->mkdir("vnObs");
      	TDirectory * vn2d = gDirectory->mkdir("vnObs2D");
      	vn1d->cd();
      	for(int k = 0; k<16; k++) {
      	  genm[i][j]->getVnObs(k)->Write(Form("vnObs_%d",k));
      	}
      	vn2d->cd();
      	for(int k = 0; k<16; k++) {
      	  genm[i][j]->getVnObs2D(k)->Write(Form("vnObs2D_%d",k));
      	}
      }
    }
  }
}
