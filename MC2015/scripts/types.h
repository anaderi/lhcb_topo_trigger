#include <iostream>

#include "TChain.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
#define data_cxx 1
#include "data.h"
#include "data.C"

data DATA;


struct Track {
  bool sig;
  float pt,p,ip,ipchi2,nvelo,nt,chi2,ismu,good;
};

void voidTrack(Track &t){
  t.pt=-1; t.p=-1; t.ip=-1; t.ipchi2=-1; t.nvelo=-1; t.nt=-1; t.chi2=-1; t.ismu=-1; t.sig=false; t.good=-1;
}

struct SV {
  int idx;
  bool sig;
  float sumpt, m, mcor, ipchi2, chi2, sumipchi2, fdr, nlt16, minpt,
    eta, pt, nmu, n, fdchi2, maxtchi2, ngood, nmu1,mupt, n1trk;
};

void voidSV(SV &sv){
  sv.idx=-1;
  sv.sig=false;
  sv.sumpt=-1; sv.m=-1; sv.mcor=-1; sv.ipchi2=-1; sv.chi2=-1; sv.sumipchi2=-1; sv.fdr=-1; sv.nlt16=-1; sv.minpt=-1;sv.eta=-1; sv.pt=-1; sv.nmu=-1; sv.n=-1; sv.chi2=-1; sv.ngood=-1; sv.nmu=-1; sv.nmu1=-1; sv.mupt=-1; sv.n1trk=-1;
}

bool isMu(int i){
  int ntrks = DATA.hlt1_muon_pt->size();
  int g = DATA.hlt1_trk1_idx_gen->at(i);
  for(int j=0; j<ntrks; j++){
    if(DATA.hlt1_muon_idx_gen->at(j) == g){
      return true;
    }
  }
  return false;
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
  t.ismu=isMu(i);
  t.sig = DATA.hlt1_trk1_idx_sig->at(i) >= 0;
  t.good = (t.chi2 < 3 && t.nvelo>=9 && t.nt>=16);
  /*  
  int g = DATA.hlt1_trk1_idx_gen->at(i);
  int gpid = DATA.gen_pid->at(g);
  if(fabs(gpid) == 13) t.ismu=1;
  */
}

void setSV(int i, SV &sv){
  TLorentzVector p4svr(0,0,0,0),p4trk;
  sv.idx=i;
  sv.n = DATA.hlt2_svr_n_trks->at(i);
  sv.sumpt = 0;
  sv.sumipchi2 = 0;
  sv.minpt = 1e9;
  sv.nlt16 = 0;
  sv.nmu = 0;
  sv.sig=true;
  sv.ngood = 0;
  sv.maxtchi2 = 0;
  sv.nmu1=0;
  sv.mupt=-1;
  sv.n1trk=0;
  for(int t = 0; t < sv.n; t++){
    int trk=-1;
    if(t==0) trk = DATA.hlt2_svr_idx_trk_1->at(i);
    if(t==1) trk = DATA.hlt2_svr_idx_trk_2->at(i);
    if(t==2) trk = DATA.hlt2_svr_idx_trk_3->at(i);
    if(t==3) trk = DATA.hlt2_svr_idx_trk_4->at(i);
    if(trk < 0) {
      cout << "Error! " << sv.n << "-body vertex but track " << t 
	   << " not found." << endl;
      continue; 
    }
    p4trk.SetXYZM(DATA.hlt2_trk_px->at(trk),
		  DATA.hlt2_trk_py->at(trk),
		  DATA.hlt2_trk_pz->at(trk),139.57);
    double pt = p4trk.Pt();
    p4svr += p4trk;
    sv.sumpt += pt;
    double ip = DATA.hlt2_trk_ip_chi2_min->at(trk); 
    sv.sumipchi2 += ip;
    if(pt < sv.minpt) sv.minpt = pt;
    if(DATA.hlt2_trk_ip_chi2_min->at(trk) < 16) sv.nlt16++;
    if(DATA.hlt2_trk_is_mu->at(trk) > 0) {
      sv.nmu++;
      if(pt > sv.mupt) sv.mupt = pt;
      if(pt > 600 && ip > 6) sv.nmu1++;
    }
    if(pt > 1000 && ip > 16) sv.n1trk++;
    if(DATA.hlt2_trk_idx_sig->at(trk) < 0) sv.sig=false;
    double tchi2 = DATA.hlt2_trk_chi2->at(trk)/DATA.hlt2_trk_ndof->at(trk);
    if(tchi2 > sv.maxtchi2) sv.maxtchi2 = tchi2;
    bool good = DATA.hlt2_trk_hits_velo->at(trk) >= 9 && 
      (DATA.hlt2_trk_hits_tt->at(trk)+DATA.hlt2_trk_hits_it->at(trk)+DATA.hlt2_trk_hits_ot->at(trk));
    if(good) sv.ngood++;
  }
  sv.m = p4svr.M();
  sv.mcor = DATA.hlt2_svr_m_cor->at(i);
  sv.ipchi2 = DATA.hlt2_svr_ip_chi2_min->at(i);
  sv.chi2 = DATA.hlt2_svr_chi2->at(i)/DATA.hlt2_svr_ndof->at(i); 
  sv.fdr = DATA.hlt2_svr_fdt_min->at(i);
  int pv = DATA.hlt2_svr_idx_pvr->at(i);
  TVector3 dir(DATA.hlt2_svr_x->at(i) - DATA.hlt2_pvr_x->at(pv),
	       DATA.hlt2_svr_y->at(i) - DATA.hlt2_pvr_y->at(pv),
	       DATA.hlt2_svr_z->at(i) - DATA.hlt2_pvr_z->at(pv));
  sv.eta = dir.PseudoRapidity();
  sv.pt = p4svr.Pt();
  sv.fdchi2 = DATA.hlt2_svr_fd_chi2->at(i);
}

bool tosLinkSV(SV &sv2, int i){
  if(i == sv2.idx) return true;
  int tfound=0;
  for(int t = 0; t < sv2.n; t++){
    int trk=-1;
    if(t==0) trk = DATA.hlt2_svr_idx_trk_1->at(sv2.idx);
    if(t==1) trk = DATA.hlt2_svr_idx_trk_2->at(sv2.idx);
    if(t==2) trk = DATA.hlt2_svr_idx_trk_3->at(sv2.idx);
    if(t==3) trk = DATA.hlt2_svr_idx_trk_4->at(sv2.idx);
    for(int tt=0; tt<2; tt++){
      int ttrk=-1;
      if(tt==0) ttrk = DATA.hlt2_svr_idx_trk_1->at(i);
      if(tt==1) ttrk = DATA.hlt2_svr_idx_trk_2->at(i);
      if(ttrk == trk){
	tfound++;
	break;
      }
    }
  }
  return tfound == 2;
}

bool tosLinkTrk(SV &sv2, int i){
  double pt_trk = DATA.hlt1_trk1_pt->at(i);
  for(int t = 0; t < sv2.n; t++){
    int trk=-1;
    if(t==0) trk = DATA.hlt2_svr_idx_trk_1->at(sv2.idx);
    if(t==1) trk = DATA.hlt2_svr_idx_trk_2->at(sv2.idx);
    if(t==2) trk = DATA.hlt2_svr_idx_trk_3->at(sv2.idx);
    if(t==3) trk = DATA.hlt2_svr_idx_trk_4->at(sv2.idx);
    TVector3 p3(DATA.hlt2_trk_px->at(trk),
		DATA.hlt2_trk_py->at(trk),
		DATA.hlt2_trk_pz->at(trk));
    if(fabs(p3.Perp()-pt_trk)/pt_trk < 0.01) return true;
  }
  return false;
}

// is truth-level b-hadron in the "good" region?
int goodGen(int gpid, int nchild){
  int good = -1;
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
    good = g;
  }
  return good;
}


