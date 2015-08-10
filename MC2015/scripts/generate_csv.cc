#include "types.h"

#include <fstream>
#include <string>

void addTrackBranches(TTree *t, Track &track){
  t->Branch("TrkPT",&track.pt);
  t->Branch("TrkP",&track.p);
  t->Branch("TrkIP",&track.ip);
  t->Branch("TrkIPCHI2",&track.ipchi2);
  t->Branch("TrkNVELO",&track.nvelo);
  t->Branch("TrkNT",&track.nt);
  t->Branch("TrkCHI2",&track.chi2);
  t->Branch("TrkISMU",&track.ismu);
  t->Branch("TrkGood",&track.good);
}

void addSVBranches(TTree *t, SV &sv){
  t->Branch("SVSumPT",&sv.sumpt);
  t->Branch("SVM",&sv.m);
  t->Branch("SVMCor",&sv.mcor);
  t->Branch("SVIPChi2",&sv.ipchi2);
  t->Branch("SVChi2",&sv.chi2);
  t->Branch("SVSumIPChi2",&sv.sumipchi2);
  t->Branch("SVFDR",&sv.fdr);
  t->Branch("SVNLT16",&sv.nlt16);
  t->Branch("SVMinPT",&sv.minpt);
  t->Branch("SVEta",&sv.eta);
  t->Branch("SVPT",&sv.pt);
  t->Branch("SVNMu",&sv.nmu);
  t->Branch("SVN",&sv.n);
  t->Branch("SVFDChi2",&sv.fdchi2);
  t->Branch("SVNGood",&sv.ngood);
}

int main(int argc, char *argv[]){

  int mode = atoi(argv[1]);
  int gpid = atoi(argv[2]);
  int nchild = atoi(argv[3]);
  std::string base_path(argv[4]);

  bool bkgd = mode==30000000;
  char str[1000];
  TChain *t = new TChain("data");
  if(!bkgd){
    sprintf(str,"%s/skims/small.%d.root", base_path.c_str(), mode);
    t->Add(str);
  }
  else{
    sprintf(str, "%s/skims/small.30000000.root", base_path.c_str());
    t->Add(str);
  }

  ofstream out_track, out_body;
  sprintf(str, "%s/prepared_hlt_track/mod_%d.csv", base_path.c_str(), mode);
  out_track.open(str);
  out_track << "unique\tmode\tevent_number\ttrk_number\tpass_trk\tsignal\tpt\tp\tip\tipchi2\tnvelo\tnt\tchi2\tismu\tgood\tsig\n";
  sprintf(str, "%s/prepared_hlt_body/mod_%d.csv", base_path.c_str(), mode);
  out_body.open(str);
  out_body << "unique\tmode\tevent_number\tsv_number\tpass_2body\tpass_nbody\tsignal\tsumpt\tm\tmcor\tipchi2\tchi2\tsumipchi2\tfdr\tnlt16\tminpt\teta\tpt\tnmu\tn\tfdchi2\tmaxtchi2\tngood\tnmu1\tmupt\tn1trk\tsig\tidx\n";

  Track track;
  SV sv;

  int nent = t->GetEntries();
  // if(!bkgd && nent > 1000) nent=1000;
  // else if(bkgd) nent=10000;

  DATA.Init(t);

  int ptrk=0,psv1=0,psv2=0;
  int no_trcks = 0;
  int no_sv = 0;
  for(int e = 0; e < nent; e++){
    DATA.GetEntry(e);
    int gen = -1;
    if(!bkgd) gen = goodGen(gpid,nchild);

    int ntrks1 = DATA.hlt1_trk1_ip->size();
    bool passtrk=false;
    if (ntrks1 == 0) {
      no_trcks++;
    }
    for(int i=0; i<ntrks1; i++){
      setTrack(i,track);
      bool cur_passed = true;
      if(!bkgd && !track.sig) cur_passed=false;
      if(track.pt < 500) cur_passed=false;
      if(track.chi2 > 3) cur_passed=false;
      if(track.ipchi2 < 4) cur_passed=false;
      out_track << mode << "_" << e << "\t"
                << mode << "\t"                  // mode
                << e << "\t"                     // event number
//                << mode << "_" << e << "\t"      // event id
                << i << "\t"                     // trk number
                << int(cur_passed) << "\t"
                << int(bkgd == 0) << "\t"
                << track.pt << "\t"
                << track.p << "\t"
                << track.ip << "\t"
                << track.ipchi2 << "\t"
                << track.nvelo << "\t"
                << track.nt << "\t"
                << track.chi2 << "\t"
                << track.ismu << "\t"
                << track.good << "\t"
                << track.sig << "\n";
      if (cur_passed) {
        passtrk = true;
      }
    }
    if(passtrk) ptrk++;

    int nsv = DATA.hlt2_svr_n_trks->size();
    if (nsv == 0) {
      no_sv++;
    }
    bool passsv_2b=false;
    bool passsv_nb=false;
    for(int i=0; i<nsv; i++){
      setSV(i,sv);
      bool cur_passed = true;
      if(!bkgd && !sv.sig) cur_passed=false;
      // if(sv.n != 2) cur_passed=false;  -- see below
      // if(sv.minpt < 500) cur_passed=false; -- totally removed
      if(sv.eta < 2 || sv.eta > 5) cur_passed=false;
      if(sv.chi2 > 10) cur_passed=false;
      if(sv.maxtchi2 > 3) cur_passed=false;
//      if(sv.mcor < 1000 || sv.mcor > 10e3) cur_passed=false;
      if(sv.mcor < 1000) cur_passed=false;
      //      if(sv.m > 7000) cout << "m: " << sv.m << endl;
      bool cur_passed_2body = cur_passed && (sv.n == 2) && (sv.minpt >= 500);
      bool cur_passed_nbody = cur_passed && (sv.nlt16 < 2);
      out_body << mode << "_" << e << "\t" 
             << mode << "\t"                  // mode
             << e << "\t"                     // event number
  //           << mode << "_" << e << "\t"      // event id
             << i << "\t"                     // trk number
             << int(cur_passed_2body) << "\t"
             << int(cur_passed_nbody) << "\t"
             << int(bkgd == 0) << "\t"
             << sv.sumpt << "\t"
              << sv.m << "\t"
              << sv.mcor << "\t"
              << sv.ipchi2 << "\t"
              << sv.chi2 << "\t"
              << sv.sumipchi2 << "\t"
              << sv.fdr << "\t"
              << sv.nlt16 << "\t"
              << sv.minpt << "\t"
              << sv.eta << "\t"
              << sv.pt << "\t"
              << sv.nmu << "\t"
              << sv.n << "\t"
              << sv.fdchi2 << "\t"
              << sv.maxtchi2 << "\t"
              << sv.ngood << "\t"
              << sv.nmu1 << "\t"
              << sv.mupt << "\t"
              << sv.n1trk << "\t"
              << sv.sig << "\t"
              << sv.idx << "\n";

      if (cur_passed_2body) {
        passsv_2b=true;
      }
      if (cur_passed_nbody) {
        passsv_nb=true;
      }

    }
    if(passsv_2b) psv1++;
    if(passsv_nb) psv2++;
  }
  cout << "Mode: " << mode << endl;
  cout << "Events: " << nent << "; No tracks: " << no_trcks << "; No sv: " << no_sv << endl;
  cout << ptrk/(double)nent << " " << psv1/(double)nent << " " << psv2/(double)nent << endl;
  return 0;
}
