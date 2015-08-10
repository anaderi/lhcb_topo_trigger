#include <iostream>

#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
#define data_cxx 1
#include "data.h"
#include "data.C"

data DATA;

bool l0Cut(double spd, double elec, double phot, double had, double mu, double dimu){
  bool l0dimu = (DATA.l0_spd_mult < 900) && (sqrt((DATA.l0_mu1_pt/1000.) * (DATA.l0_mu2_pt/1000.)) > dimu/1000.);
  if(!l0dimu && DATA.l0_spd_mult > spd) return false;
  bool l0elec = DATA.l0_elc_et > elec;
  bool l0phot = DATA.l0_pho_et > phot;
  bool l0had = DATA.l0_had_et > had;
  bool l0mu = DATA.l0_mu1_pt > mu;
  //cout << DATA.l0_spd_mult << " " << DATA.l0_mu1_pt << " " << DATA.l0_elc_et << " " << DATA.l0_pho_et << " " <<  DATA.l0_had_et << " " << endl;
  return (l0elec || l0phot || l0had || l0mu || l0dimu);
}

// implement various BW division results here (there will only be a few)
bool l0Pass(int which=0){
  return l0Cut(300,153*20,177*20,202*20,50*40,951*40);
  //return l0Cut(600,179*20,246*20,248*20,62*40,951*40);
}

/*
struct Track {
  double pt,p,ip,ipchi2,nvelo,nt,chi2,ismu;
};

void voidTrack(Track &t){
  t.pt=-1; t.p=-1; t.ip=-1; t.ipchi2=-1; t.nvelo=-1; t.nt=-1; t.chi2=-1; t.ismu=-1;
}

void setTrack(int i, Track &t){
  int ntrks1 = DATA.hlt1_trk1_ip->size();
  int ntrks2 = DATA.hlt2_trk_px->size();
  voidTrack(t);
  if(i < 0 || i >= ntrks1) return;
  t.pt =DATA.hlt1_trk1_pt->at(i);
  t.p = DATA.hlt1_trk1_p->at(i);
  t.ip = DATA.hlt1_trk1_ip->at(i);
  t.ipchi2 = DATA.hlt1_trk1_ip_chi2->at(i);
  t.nvelo = DATA.hlt1_trk1_velo_hits->at(i);
  t.nt = DATA.hlt1_trk1_t_hits->at(i);
  t.chi2 = DATA.hlt1_trk1_chi2_ndof->at(i);
  t.ismu=0;
  int g = DATA.hlt1_trk1_idx_gen->at(i);
  int gpid = DATA.gen_pid->at(g);
  if(fabs(gpid) == 13) t.ismu=1;
}
*/

// is truth-level b-hadron in the "good" region?
bool goodGen(int gpid, int nchild){
  int ngood = 0;
  int ngen = DATA.gen_pid->size();
  int npv = DATA.hlt2_pvr_z->size();
  for(int g = 0; g < ngen; g++){
    int pid = fabs(DATA.gen_pid->at(g));
    if(pid != gpid) continue;
    TLorentzVector p4(DATA.gen_px->at(g),DATA.gen_py->at(g),DATA.gen_pz->at(g),
		      DATA.gen_e->at(g));
    double pt = p4.Pt()/1000.;
    if(p4.Pt() < 2000) continue;
    if(p4.PseudoRapidity() < 2) continue;
    if(p4.PseudoRapidity() > 5) continue;
    int gf = DATA.gen_idx_childf->at(g);
    int gl = DATA.gen_idx_childl->at(g);
    if(gf < 0 || gf >= ngen) continue;
    if(gl < 0 || gl >= ngen) continue;
    if(gl-gf != nchild) continue;
    double x = DATA.gen_x->at(g), y = DATA.gen_y->at(g), z = DATA.gen_z->at(g);
    double fdrmax=0,fdmin = 1e6;
    TLorentzVector p1,p2,p3;
    int pidx=0;
    for(int gg = gf; gg <= gl; gg++){
      double xx = DATA.gen_x->at(gg), yy = DATA.gen_y->at(gg),
	zz = DATA.gen_z->at(gg);
      TVector3 fd3d(xx-x,yy-y,zz-z);
      double fdr = fd3d.Perp();
      if(fdr > fdrmax) fdrmax=fdr;
      if(fd3d.Mag() > 0 && fd3d.Mag() < fdmin) fdmin = fd3d.Mag();
      if(pidx==0) p1.SetPxPyPzE(DATA.gen_px->at(gg),DATA.gen_py->at(gg),DATA.gen_pz->at(gg),DATA.gen_e->at(gg));
      if(pidx==1) p2.SetPxPyPzE(DATA.gen_px->at(gg),DATA.gen_py->at(gg),DATA.gen_pz->at(gg),DATA.gen_e->at(gg));
      if(pidx==2) p3.SetPxPyPzE(DATA.gen_px->at(gg),DATA.gen_py->at(gg),DATA.gen_pz->at(gg),DATA.gen_e->at(gg));
      pidx++;
    }
    double tau = (fdmin/(3e11))*(p4.M()/p4.P()); // tau in s
    tau *= 1e12;
    if(tau < 0.2) continue;
    ngood++;
  }
  return ngood > 0;
}



int main(int argc, char *argv[]){

  int mode = atoi(argv[1]);
  int gpid = atoi(argv[2]);
  int nchild = atoi(argv[3]);
  std::string base_path(argv[4]);
  bool bkgd = mode==30000000;

  char str[100];
  TChain *t = new TChain("data");
  if(!bkgd){
    sprintf(str,"%s/Hlt.MC15.MD.%d.1.150204.00.root",base_path.c_str(), mode);
    t->Add(str);
    sprintf(str,"%s/Hlt.MC15.MU.%d.1.150204.00.root",base_path.c_str(), mode);
    t->Add(str);
  }
  else{
    sprintf(str, "%s/Hlt.MC15.MD.30000000.0.150212.00.root", base_path.c_str());
    t->Add(str);
    sprintf(str, "%s/Hlt.MC15.MU.30000000.0.150212.00.root", base_path.c_str());
    t->Add(str);
  }

  sprintf(str,"%s/skims/small.%d.root", base_path.c_str(), mode);
  TFile fout(str,"recreate");
  TTree *tout = t->CloneTree(0);

  int nent = t->GetEntries();
  cout << "n(ent): " << nent << endl;
  DATA.Init(t);

  int count=0;
  for(int e = 0; e < nent; e++){
    DATA.GetEntry(e);
    if(!bkgd && !goodGen(gpid,nchild)) continue;
    if(!l0Pass()) continue;
    //    cout << "---" << endl;
    /*
    int ntrks1 = DATA.hlt1_trk1_ip->size();
    Track track;
    for(int i=0; i<ntrks1; i++){
      setTrack(i,track);
    }
    */
    tout->Fill();
    count++;
    // if(!bkgd && count >= 1e3) break;
    // if(bkgd && count >= 1e4) break;
  }
  cout << count << endl;
  fout.cd();
  tout->Write();

  return 0;
}
