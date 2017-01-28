if(reac.Contains("pPb")){
  prename = "rpflat";
  ep3a = hi::HFp3;
  ep3b = hi::HFm3;
  ep3c = hi::trackmid3;
  mid3n = "trackmid3";
  for(int i = 0; i<nanals; i++) if(AnalNames[i] == anal) {ANAL = i; break;}
  if(anal=="N2SUB2" || anal=="N2SUB3" || anal=="D24SUB2" || anal=="D24SUB3") {
    epord_ = 2;
    trig  = "crab_PbPb_2015_v2_A";
    trig2 = "crab_PbPb_2015_v2_B";
    trig3 = "crab_PbPb_2015_v2_C";
    mid2n = "trackmid2";
    epa = hi::HFp2;
    epb = hi::HFm2;
    epc = hi::trackmid2;
  }
  if(anal=="N3SUB2" || anal=="N3SUB3" || anal=="D34SUB2" || anal == "D34SUB3") {
    epord_ = 3;
    trig  = "crab_PbPb_2015_v3_A";
    trig2 = "crab_PbPb_2015_v3_B";
    trig3 = "crab_PbPb_2015_v3_C";
    mid2n = "trackmid3";
    epa = hi::HFp3;
    epb = hi::HFm3;
    epc = hi::trackmid3;
  }
  if(anal=="N4SUB2" || anal=="N4SUB3" ) {
    epord_ = 4;
    trig  = "crab_PbPb_2015_v4_A";
    trig2 = "crab_PbPb_2015_v4_B";
    trig3 = "crab_PbPb_2015_v4_C";
    mid2n = "trackmid4";
    epa = hi::HFp4;
    epb = hi::HFm4;
    epc = hi::trackmid4;
  }
  if(anal=="N5SUB2" || anal=="N5SUB3") {
    epord_ = 5;
    trig  = "crab_PbPb_2015_v5_A";
    trig2 = "crab_PbPb_2015_v5_B";
    trig3 = "crab_PbPb_2015_v5_C";
    mid2n = "trackmid5";
    epa = hi::HFp5;
    epb = hi::HFm5;
    epc = hi::trackmid5;
  }
  if(anal=="N6SUB2" || anal=="N6SUB3") {
    epord_ = 6;
    trig  = "crab_PbPb_2015_v6_A";
    trig2 = "crab_PbPb_2015_v6_B";
    trig3 = "crab_PbPb_2015_v6_C";
    mid2n = "trackmid6";
    epa = hi::HFp6;
    epb = hi::HFm6;
    epc = hi::trackmid6;
  }
  if(anal=="N7SUB2" || anal=="N7SUB3") {
    epord_ = 7;
    trig  = "crab_PbPb_2015_v7_A";
    trig2 = "crab_PbPb_2015_v7_B";
    trig3 = "crab_PbPb_2015_v7_C";
    mid2n = "trackmid7";
    epa = hi::HFp7;
    epb = hi::HFm7;
    epc = hi::trackmid7;
  }
  if(anal=="N42SUB2" || anal=="N42SUB3" ) {
    epord_ = 4;
    trig  = "crab_PbPb_2015_v4_A";
    trig2 = "crab_PbPb_2015_v4_B";
    trig3 = "crab_PbPb_2015_v4_C";
    mid2n = "trackmid2";
    epa = hi::HFp2;
    epb = hi::HFm2;
    epc = hi::trackmid2;
  }
  if(anal=="N523SUB2" || anal=="N523SUB3" ) {
    epord_ = 5;
    trig  = "crab_PbPb_2015_v5_A";
    trig2 = "crab_PbPb_2015_v5_B";
    trig3 = "crab_PbPb_2015_v5_C";
    mid2n = "trackmid2";
    epa = hi::HFp2;
    epb = hi::HFm2;
    epc = hi::trackmid2;
  }
  if(anal=="N723SUB2" || anal=="N723SUB3" ) {
    epord_ = 7;
    trig  = "crab_PbPb_2015_v7_A";
    trig2 = "crab_PbPb_2015_v7_B";
    trig3 = "crab_PbPb_2015_v7_C";
    mid2n = "trackmid2";
    epa = hi::HFp2;
    epb = hi::HFm2;
    epc = hi::trackmid2;
  }
  if(anal=="N63SUB2" || anal=="N63SUB3" ) {
    epord_ = 6;
    trig  = "crab_PbPb_2015_v6_A";
    trig2 = "crab_PbPb_2015_v6_B";
    trig3 = "crab_PbPb_2015_v6_C";
    mid2n = "trackmid3";
    epa = hi::HFp3;
    epb = hi::HFm3;
    epc = hi::trackmid3;
  }
  if(anal=="D2232SUB2" || anal=="D2232SUB3" ) {
    epord_ = 2;
    trig  = "crab_PbPb_2015_v2_v3_A";
    trig2 = "crab_PbPb_2015_v2_v3_B";
    trig3 = "crab_PbPb_2015_v2_v3_C";
    mid2n = "trackmid2";
    epa = hi::HFp2;
    epb = hi::HFm2;
    epc = hi::trackmid2;
  }
 }
