#include "global.h"

void getval(double npass, double nall, double &p, double &pe){
  p = npass/nall;
  pe = sqrt(p*(1-p)/nall);
}

int main(int argc, char *argv[]){
  
  int mode = atoi(argv[1]);
  int bpid = atoi(argv[2]);
  int nbchild = atoi(argv[3]);
  bool bkgd = mode==30000000;
  int l0 = 0;
  int hlt1 = 0;
  char str[100];

  sprintf(str,"plots.%d.root",mode);
  TFile fout(str,"recreate");
  const int nptbins = 7;
  double ptbins[8] = {2,4,6,8,10,15,25,50};
  TH1F *hpt = new TH1F("hpt","",nptbins,ptbins);
  TH1F *hptbdt1 = new TH1F("hptbdt1","",nptbins,ptbins);
  TH1F *hptbdt2 = new TH1F("hptbdt2","",nptbins,ptbins);
  TH1F *hptmn = new TH1F("hptmn","",nptbins,ptbins);
  TH1F *hptada = new TH1F("hptada","",nptbins,ptbins);
  TH1F *hptfor = new TH1F("hptfor","",nptbins,ptbins);
  TH1F *hptugb = new TH1F("hptugb","",nptbins,ptbins);

  const int ntbins =6;
  double tbins[7] = {0,0.25,0.5,1,2.5,5,10};
  TH1F *ht = new TH1F("ht","",ntbins,tbins);
  TH1F *htbdt1 = new TH1F("htbdt1","",ntbins,tbins);
  TH1F *htbdt2 = new TH1F("htbdt2","",ntbins,tbins);
  TH1F *htmn = new TH1F("htmn","",ntbins,tbins);
  TH1F *htada = new TH1F("htada","",ntbins,tbins);
  TH1F *htfor = new TH1F("htfor","",ntbins,tbins);
  TH1F *htugb = new TH1F("htugb","",ntbins,tbins);

  TH2F *hdp = new TH2F("hdp","",10,0.5,25,10,0.5,25);
  TH2F *hdpbdt1 = new TH2F("hdpbdt1","",10,0.5,25,10,0.5,25);
  TH2F *hdpmn = new TH2F("hdpmn","",10,0.5,25,10,0.5,25);

  int tidx=2;
  if(mode == 11114001 || mode == 11124001 || mode == 11296013 
     || mode == 11874042 || mode == 12103035 || mode == 12165106 ||
     mode == 12265042) tidx = 3;
  sprintf(str,"data;%d",tidx);
  TChain *t = new TChain(str);
  
  //sprintf(str,"Hlt.MC15.MD.%d.0.140818.00.root",mode);
  sprintf(str,"data.new/Hlt.MC15.MD.%d.0.140824.00.root",mode);
  t->Add(str);
  //sprintf(str,"Hlt.MC15.MU.%d.0.140818.00.root",mode);
  sprintf(str,"data.new/Hlt.MC15.MU.%d.0.140824.00.root",mode);
  t->Add(str);
  
  float bdt_mn[10000];
  t->SetBranchAddress("bdt_MN",&bdt_mn);
  float bdt_ada[10000];
  t->SetBranchAddress("bdt_ada",&bdt_ada);
  float bdt_forest[10000];
  t->SetBranchAddress("bdt_forest",&bdt_forest);
  float bdt_ugb[10000];
  t->SetBranchAddress("bdt_ugb_fl",&bdt_ugb);
  int train;
  if(tidx == 3) t->SetBranchAddress("train_test_label",&train);

  // cuts
  double mncut = 0.925795006942;
  double adacut = 0.600167561647;
  double forcut = 0.874676230445;
  double ugbcut = 0.777668282603;

  //t->Add("output.root");
  int nent = t->GetEntries();
  DATA.Init(t);
  SVR svr;
  vector<int> hlt1trks;

  int nl0=0,nhlt1=0,nhlt2_run1=0,nhlt2_pre=0,nhlt2_run1mod=0;
  int nmn=0,nada=0,nfor=0,nugb=0;
  //cout << nent << endl;
  for(int e = 0; e < nent; e++){
    DATA.GetEntry(e);
    if(tidx == 3 && train == 1) continue;
    double bpt,btau,bm12,bm13;
    if(!bkgd && !goodGenB(bpid,nbchild,bpt,btau,bm12,bm13)) continue;
    // now btau is the true lifetime of the b
    //cout << btau << endl;
    if(!l0Pass(l0)) continue;
    nl0++;
    //if(!hlt1Pass(hlt1,hlt1trks)) continue;
    hlt1Pass(hlt1,hlt1trks);
    nhlt1++;
    hpt->Fill(bpt);
    ht->Fill(btau); 
    hdp->Fill(bm12*bm12,bm13*bm13);

    int nsvr = DATA.hlt2_svr_chi2->size();
    int npassrun1=0,npassrun1mod=0,npasspre=0;
    int nmnpass=0,nadapass=0,nforpass=0,nugbpass=0;
    int npasspretrue=0;
    //cout << "---" << endl;
    for(int s = 0; s < nsvr; s++){
      setSVR(bpid,nbchild,svr, s, hlt1trks);
      //cout << s << " " << bdt_mn[s] << endl;
      if(svrPassRun1(svr)) npassrun1++;
      if(svrPassRun1Modified(svr)) npassrun1mod++;
      if(svrPreSel(svr)){
	npasspre++;
	if(svr.truematch > 0) npasspretrue++;
	if(bdt_mn[s] > mncut) nmnpass++;
	if(bdt_ada[s] > adacut) nadapass++;
	if(bdt_forest[s] > forcut) nforpass++;
	if(bdt_ugb[s] > ugbcut) nugbpass++;
      }
    }
    //cout << npasspre << " " << npasspretrue << endl;

    if(npassrun1 > 0) {
      nhlt2_run1++;
      hptbdt1->Fill(bpt);
      htbdt1->Fill(btau);
      hdpbdt1->Fill(bm12*bm12,bm13*bm13);
    }
    if(npassrun1mod > 0) {
      nhlt2_run1mod++;
      hptbdt2->Fill(bpt);
      htbdt2->Fill(btau);
    }
    if(npasspre > 0) nhlt2_pre++;
    //if(nsvr > 0) nhlt2_pre++;
    if(nmnpass > 0) {
      nmn++;
      hptmn->Fill(bpt);
      htmn->Fill(btau);
      hdpmn->Fill(bm12*bm12,bm13*bm13);
    }
    if(nadapass > 0) {
      nada++;
      hptada->Fill(bpt);
      htada->Fill(btau);
    }
    if(nforpass > 0){
      nfor++;
      hptfor->Fill(bpt);
      htfor->Fill(btau);
    }
    if(nugbpass > 0){
      nugb++;
      hptugb->Fill(bpt);
      htugb->Fill(btau);
    }
  }
  double l0_rate = 1e6; 
  double rate_hlt1 =  (nhlt1/(double)nl0)*l0_rate;
  double prun1,prun1e,prun1mod,prun1mode,ppre,ppree;
  getval(nhlt2_run1,nhlt1,prun1,prun1e);
  getval(nhlt2_run1mod,nhlt1,prun1mod,prun1mode);
  getval(nhlt2_pre,nhlt1,ppre,ppree);
  double pmn,pmne,pada,padae,pfor,pfore,pugb,pugbe;
  getval(nmn,nhlt1,pmn,pmne);
  getval(nada,nhlt1,pada,padae);
  getval(nfor,nhlt1,pfor,pfore);
  getval(nugb,nhlt1,pugb,pugbe);

  cout << nl0 << " " << nhlt1 << endl;
  if(bkgd){
    cout << rate_hlt1 << " " << prun1*rate_hlt1 << "+-" << prun1e*rate_hlt1 
	 << " " << prun1mod*rate_hlt1 << "+-" << prun1mode*rate_hlt1 
	 << " " << ppre*rate_hlt1 << "+-" << ppree*rate_hlt1 
	 << endl;
    cout << pmn*rate_hlt1 << "+-" << pmne*rate_hlt1 << endl;
    cout << pada*rate_hlt1 << "+-" << padae*rate_hlt1 << endl;
    cout << pfor*rate_hlt1 << "+-" << pfore*rate_hlt1 << endl;
    cout << pugb*rate_hlt1 << "+-" << pugbe*rate_hlt1 << endl;
  }
  else{
    cout << prun1 << "+-" << prun1e 
	 << " " << prun1mod << "+-" << prun1mode 
	 << " " << ppre << "+-" << ppree 
	 << endl;
    cout << pmn << "+-" << pmne << endl;
    cout << pada << "+-" << padae << endl;
    cout << pfor << "+-" << pfore << endl;
    cout << pugb << "+-" << pugbe << endl;
  }

  fout.cd();
  hpt->Write();
  hptbdt1->Write();
  hptbdt2->Write();
  hptmn->Write();
  hptada->Write();
  hptfor->Write();
  hptugb->Write();

  ht->Write();
  htbdt1->Write();
  htbdt2->Write();
  htmn->Write();
  htada->Write();
  htfor->Write();
  htugb->Write();

  hdp->Write();
  hdpbdt1->Write();
  hdpmn->Write();

  return 0;
}

