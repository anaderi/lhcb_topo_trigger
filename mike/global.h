#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
#define data_cxx 1
#include "data.h"
#include "data.C"

data DATA;

bool l0Cut(double spd, double elec, double phot, double had, double mu, double dimu){
  bool l0dimu = (DATA.l0_spd_mult < 900) && (sqrt((DATA.l0_mu1_pt/1000.) * (DATA.l0_mu2_pt/1000.)) > dimu/1000.);//(DATA.l0_dmu_pt > dimu);//
  if(!l0dimu && DATA.l0_spd_mult > spd) return false;
  bool l0elec = DATA.l0_elc_et > elec;
  bool l0phot = DATA.l0_pho_et > phot;
  bool l0had = DATA.l0_had_et > had;
  bool l0mu = DATA.l0_mu1_pt > mu;
  //  return l0phot;
  return (l0elec || l0phot || l0had || l0mu || l0dimu);
}

// implement various BW division results here (there will only be a few)
bool l0Pass(int which=0){
  if(which == 0) return l0Cut(600,4000,4000,4200,2400,6160);
  return true;
}

double trackMatch12(int trk1, int trk2){
  TVector3 p(DATA.hlt2_trk_px->at(trk2),DATA.hlt2_trk_py->at(trk2),
	     DATA.hlt2_trk_pz->at(trk2));
  double dpt = fabs(DATA.hlt1_trk1_pt->at(trk1)-p.Pt());
  double dp = fabs(DATA.hlt1_trk1_p->at(trk1)-p.Mag());
  return sqrt(dpt*dpt + dp*dp)/p.Mag();
}

bool hlt1Cut(double ip, double nvelo, double nt, double chi2ndof, 
	     double ipchi2, double p, double pt, 
	     vector<int> &pass_trks){
  int ntrks = DATA.hlt1_trk1_ip->size();
  int ntrks2 = DATA.hlt2_trk_px->size();
  for(int i = 0; i < ntrks; i++){
    // find hlt2 track
    int hlt2idx = -1; double best = 1e10;
    for(int j = 0; j < ntrks2; j++){
      double match = trackMatch12(i,j);
      if(match < best){
	best = match;
	hlt2idx = j;
      }
    }
    bool is_mu = hlt2idx >=0 && DATA.hlt2_trk_is_mu->at(hlt2idx);    
    if(ip > 0 && DATA.hlt1_trk1_ip->at(i) < ip) continue;
    if(nvelo > 0 && DATA.hlt1_trk1_velo_hits->at(i) < nvelo) continue;
    if(nt > 0 && DATA.hlt1_trk1_t_hits->at(i) < nt) continue;
    if(chi2ndof > 0 && DATA.hlt1_trk1_chi2_ndof->at(i) > chi2ndof) continue;
    if(DATA.hlt1_trk1_ip_chi2->at(i) < ipchi2) continue;
    if(DATA.hlt1_trk1_p->at(i) < p) continue;
    if(DATA.hlt1_trk1_pt->at(i) < pt) continue;
    pass_trks.push_back(i);
  }
  return pass_trks.size() > 0;
}

// implement various hlt1 scenarios here
bool hlt1Pass(int which, vector<int> &pass_trks){
  pass_trks.clear();
  if(which == 0) return hlt1Cut(-1,9,16,2,25,3000,2500,pass_trks);
  return true;
}

struct SVR {
  double sumpt, m, mcor, ipchi2, vchi2, sumipchi2, fdr, nlt16, minpt,
    ptau, q, eta, pt, nhlt1, nmu, n, bdt, fdchi2;
  double truematch;
};

void setSVR(int bpid, int nbchild, SVR &svr, int s, vector<int> &hlt1trks){
  TLorentzVector p4svr(0,0,0,0),p4trk;
  int ntrk = DATA.hlt2_svr_n_trks->at(s);
  svr.n = ntrk;
  svr.minpt = 1e12;
  svr.sumpt = 0;
  svr.nlt16 = 0;
  svr.q = 0;
  svr.nmu = 0;
  svr.nhlt1 = 0;
  svr.truematch = 0;
  int nhlt1 = hlt1trks.size();
  int trkmax = DATA.hlt2_trk_px->size();
  double sumip=0;
  int ntrkmatchb=0;
  double genpt[4] = {0,0,0,0};
  for(int t = 0; t < ntrk; t++){
    int trk=-1;
    if(t==0) trk = DATA.hlt2_svr_idx_trk_1->at(s);
    if(t==1) trk = DATA.hlt2_svr_idx_trk_2->at(s);
    if(t==2) trk = DATA.hlt2_svr_idx_trk_3->at(s);
    if(t==3) trk = DATA.hlt2_svr_idx_trk_4->at(s);
    if(trk < 0 || trk >= trkmax) {
      cout << "Error! " << svr.n << "-body vertex but track " << t 
	   << " not found." << endl;
      continue; 
    }
    if(DATA.hlt2_trk_dr_sig->at(trk) < 0.05) ntrkmatchb++;
    int gdx = DATA.hlt2_trk_idx_gen->at(trk);
    /*
    if(gdx >= 0 && gdx < DATA.gen_pid->size()){
      if(DATA.gen_idx_mother->at(gdx) >= 0) ntrkmatchb++;
    }
    */
    p4trk.SetXYZM(DATA.hlt2_trk_px->at(trk),
		  DATA.hlt2_trk_py->at(trk),
		  DATA.hlt2_trk_pz->at(trk),139.57);
    double pt = p4trk.Pt();
    genpt[t]=TVector3(DATA.gen_px->at(gdx),DATA.gen_py->at(gdx),DATA.gen_pz->at(gdx)).Pt();
    p4svr += p4trk;
    svr.sumpt += pt;
    double ip = DATA.hlt2_trk_ip_chi2_min->at(trk); 
    sumip += ip;
    if(pt < svr.minpt) svr.minpt = pt;
    if(DATA.hlt2_trk_ip_chi2_min->at(trk) < 16) svr.nlt16++;
    svr.q += DATA.hlt2_trk_charge->at(trk);
    if(DATA.hlt2_trk_is_mu->at(trk) > 0) svr.nmu++;

    /*
    // did this track pass hlt1?
    for(int t1 = 0; t1 < nhlt1; t1++){
      if(trackMatch12(hlt1trks[t1],trk) < 0.1){
	svr.nhlt1++;
	break;
      }
    }
    */
    if(pt > 1700 && ip > 16) svr.nhlt1++;
    
  }
  svr.m = p4svr.M();
  svr.mcor = DATA.hlt2_svr_m_cor->at(s);
  svr.ipchi2 = DATA.hlt2_svr_ip_chi2_min->at(s);
  svr.vchi2 = DATA.hlt2_svr_chi2->at(s);
  svr.sumipchi2 = DATA.hlt2_svr_ip_chi2_sum->at(s);
  svr.fdr = DATA.hlt2_svr_fdt_min->at(s);
  svr.ptau = DATA.hlt2_svr_fd->at(s)/p4svr.P()*1000;
  int pv = DATA.hlt2_svr_idx_pvr->at(s);
  TVector3 dir(DATA.hlt2_svr_x->at(s) - DATA.hlt2_pvr_x->at(pv),
	       DATA.hlt2_svr_y->at(s) - DATA.hlt2_pvr_y->at(pv),
	       DATA.hlt2_svr_z->at(s) - DATA.hlt2_pvr_z->at(pv));
  svr.eta = dir.PseudoRapidity();
  svr.pt = p4svr.Pt();
  svr.bdt = DATA.hlt2_svr_bdt->at(s);
  svr.fdchi2 = DATA.hlt2_svr_fd_chi2->at(s);

  svr.truematch = 0;
  if(ntrkmatchb == svr.n) svr.truematch = 1;
}

bool svrPassRun1(SVR &svr){
  double bdtcut = 0.4;
  if(svr.n == 4) bdtcut = 0.3;
  if(svr.nmu > 0) bdtcut = 0.1;
  if(svr.bdt < bdtcut) return false;
  if(svr.fdchi2 < 100) return false;
  if(svr.sumpt < 4000) return false;
  if(svr.n > 2 && svr.sumpt < 4500) return false;
  if(svr.minpt < 500) return false;
  if(svr.nhlt1 == 0) return false;
  return true;
}

bool svrPassRun1Modified(SVR &svr){
  double bdtcut = 0.3;
  if(svr.n == 4) bdtcut = 0.3;
  if(svr.nmu > 0) bdtcut = 0.1;
  if(svr.bdt < bdtcut) return false;
  if(svr.fdchi2 < 32) return false;
  return true;
}

// pre-BDT selection
bool svrPreSel(SVR &svr) {
  if(svr.eta < 2 || svr.eta > 5) return false;
  if(svr.fdr < 0.1) return false;
  if(svr.mcor > 10e3) return false;
  if(svr.pt < 2e3) return false;
  if(svr.sumpt < 3e3) return false;
  if(svr.sumipchi2 < svr.n*25) return false;
  if(svr.nlt16 > 1) return false;
  //  if(svr.ptau > 2) return false;
  return true;
}

// is truth-level b-hadron in the "good" region?
bool goodGenB(int bpid, int nbchild, double &pt, double &tau, double &m12, double &m13, double ptcut=2000, double taucut=0.2){
  int ngood = 0;
  int ngen = DATA.gen_pid->size();
  int npv = DATA.hlt2_pvr_z->size();
  for(int g = 0; g < ngen; g++){
    int pid = fabs(DATA.gen_pid->at(g));
    bool isb = (pid==511 || pid==521 || pid==531 || pid==5122);
    if(!isb) continue;
    if(pid != bpid) continue;
    TLorentzVector p4(DATA.gen_px->at(g),DATA.gen_py->at(g),DATA.gen_pz->at(g),
		      DATA.gen_e->at(g));
    pt = p4.Pt()/1000.;
    if(p4.Pt() < ptcut) continue;
    if(p4.PseudoRapidity() < 2) continue;
    if(p4.PseudoRapidity() > 5) continue;    
    int npt500=0;
    int gf = DATA.gen_idx_childf->at(g);
    int gl = DATA.gen_idx_childl->at(g);
    if(gf < 0 || gf >= ngen) continue;
    if(gl < 0 || gl >= ngen) continue;
    if(gl-gf != nbchild) continue;
    double x = DATA.gen_x->at(g), y = DATA.gen_y->at(g), z = DATA.gen_z->at(g);
    double fdrmax=0,fdmin = 1e6;
    TLorentzVector p1,p2,p3;
    int pidx=0;
    for(int gg = gf; gg <= gl; gg++){
      if(TVector3(DATA.gen_px->at(gg),DATA.gen_py->at(gg),DATA.gen_pz->at(gg)).Pt() > 500) npt500++;
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
    //if(npt500 < 2) continue;
    tau = (fdmin/(3e11))*(p4.M()/p4.P()); // tau in s
    tau *= 1e12;
    m12 = (p1+p2).M()/1000.;
    m13 = (p1+p3).M()/1000.;
    if(tau < taucut) continue;
    //if(fdrmax < fdrcut) continue;
    ngood++;
    //    cout << tau << " " << fdmin << endl;
  }
  return ngood > 0;
}

