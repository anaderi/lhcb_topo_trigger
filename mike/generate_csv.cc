#include "global.h"
#include <fstream>
#include <string>
#include <iostream>

using std::cout;

int main(int argc, char *argv[]){

  int mode = atoi(argv[1]);
  std::string path("data");
  if (argc >= 5) {
    path = std::string(argv[4]);
  }

  int bpid = atoi(argv[2]);
  int nbchild = atoi(argv[3]);

  std::string dest_path("csv_data");
  if (argc == 6) {
    dest_path = std::string(argv[5]);
  }
  bool bkgd = mode==30000000;
  int l0 = 0;
  int hlt1 = 0;
  char str[100];
  TChain *t = new TChain("data");
  printf("%d\n", mode);
  //sprintf(str,"Hlt.MC15.MD.%d.0.140818.00.root",mode);
  sprintf(str,"%s/Hlt.MC15.MD.%d.0.140824.00.root", path.c_str(), mode);
  t->Add(str);
  //sprintf(str,"Hlt.MC15.MU.%d.0.140818.00.root",mode);
  sprintf(str,"%s/Hlt.MC15.MU.%d.0.140824.00.root", path.c_str(), mode);
  t->Add(str);

  //t->Add("output.root");
  int nent = t->GetEntries();
  DATA.Init(t);
  SVR svr;
  vector<int> hlt1trks;

  sprintf(str, "%s/mod_%d.csv", dest_path.c_str(), mode);
  ofstream out(str);
  out << "# event_id\tsvr\tsignal\tpresel\ttruematch\tsumpt\tm\tmcor\tipchi2\tvchi2\tsumipchi2\tfdr\tnlt16\tminpt\tptau\tq\teta\tpt\tnhlt1\tnmu\tn\tbdt\tfdchi2\tNUMBER\tbtau\tis_passed_l0_goodGenB\n";

  int nl0=0,nhlt1=0,nhlt2_run1=0,nhlt2_pre=0,nhlt2_run1mod=0;
  int nmn=0,nada=0,nfor=0,nugb=0;
  //cout << nent << endl;
  int empty_good = 0;
  int empty_bad = 0;
  for(int e = 0; e < nent; e++){
    DATA.GetEntry(e);
    bool is_passed = false;
    double bpt,btau,bm12,bm13;
    btau = -1;
    if (!(!bkgd && !goodGenB(bpid,nbchild,bpt,btau,bm12,bm13))
       && l0Pass(l0)) {
      nl0++;
      nhlt1++;
      is_passed = true;
    }
    // now btau is the true lifetime of the b
    //cout << btau << endl;
    //if(!hlt1Pass(hlt1,hlt1trks)) continue;
    hlt1Pass(hlt1,hlt1trks);
    int nsvr = DATA.hlt2_svr_chi2->size();
    if (nsvr == 0) {
      if (is_passed) {
         empty_good += 1;
      } else {
         empty_bad += 1;
      }
    }
    for(int s = 0; s < nsvr; s++){
      setSVR(bpid,nbchild,svr, s, hlt1trks);

      if (is_passed) {
        out << mode << "_" << nhlt1; // printing event_id: file_id + "_" + event_number;
      } else {
        out << mode << "_" << "not_passed";
      }
      out << "\t" << s;
      out << "\t" << int(bkgd == 0); // printing is_signal
      out << "\t" << int(svrPreSel(svr) == true); // passed signal
      out << "\t" << svr.truematch;
      out << "\t" << svr.sumpt
          << "\t" << svr.m
          << "\t" << svr.mcor
          << "\t" << svr.ipchi2
          << "\t" << svr.vchi2
          << "\t" << svr.sumipchi2
          << "\t" << svr.fdr
          << "\t" << svr.nlt16
          << "\t" << svr.minpt
          << "\t" << svr.ptau
          << "\t" << svr.q
          << "\t" << svr.eta
          << "\t" << svr.pt
          << "\t" << svr.nhlt1
          << "\t" << svr.nmu
          << "\t" << svr.n
          << "\t" << svr.bdt
          << "\t" << svr.fdchi2
          << "\t" << e
          << "\t" << btau
	  << "\t" << is_passed << "\n";
    }
  }
  cout << "Total number: " << nent << std::endl;
  cout << "Passed events: " << nhlt1 << std::endl;
  cout << "Empty passed: " << empty_good << std::endl;
  cout << "Empty not passed: " << empty_bad << std::endl;
  return 0;
}

