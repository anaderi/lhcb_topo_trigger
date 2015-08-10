//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 12 18:15:43 2015 by ROOT version 5.34/05
// from TTree data/data
// found on file: Hlt.MC15.MD.11114001.1.150204.00.root
//////////////////////////////////////////////////////////

#ifndef data_h
#define data_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        l0_spd_mult;
   Double_t        l0_elc_et;
   Double_t        l0_pho_et;
   Double_t        l0_had_et;
   Double_t        l0_mu1_pt;
   Double_t        l0_mu2_pt;
   Double_t        l0_dmu_pt;
   vector<double>  *hlt1_trk0_ip;
   vector<double>  *hlt1_trk0_ip_chi2;
   vector<double>  *hlt1_trk0_chi2_ndof;
   vector<double>  *hlt1_trk0_t_hits;
   vector<double>  *hlt1_trk0_velo_hits;
   vector<double>  *hlt1_trk0_velo_q;
   vector<double>  *hlt1_trk1_idx_gen;
   vector<double>  *hlt1_trk1_idx_sig;
   vector<double>  *hlt1_trk1_dr_gen;
   vector<double>  *hlt1_trk1_dr_sig;
   vector<double>  *hlt1_trk1_pt;
   vector<double>  *hlt1_trk1_p;
   vector<double>  *hlt1_trk1_ip;
   vector<double>  *hlt1_trk1_ip_chi2;
   vector<double>  *hlt1_trk1_chi2_ndof;
   vector<double>  *hlt1_trk1_t_hits;
   vector<double>  *hlt1_trk1_velo_hits;
   vector<double>  *hlt1_trk1_velo_q;
   vector<double>  *hlt1_muon_idx_gen;
   vector<double>  *hlt1_muon_idx_sig;
   vector<double>  *hlt1_muon_dr_gen;
   vector<double>  *hlt1_muon_dr_sig;
   vector<double>  *hlt1_muon_pt;
   vector<double>  *hlt1_muon_p;
   vector<double>  *hlt1_muon_ip;
   vector<double>  *hlt1_muon_ip_chi2;
   vector<double>  *hlt1_muon_chi2_ndof;
   vector<double>  *hlt1_muon_t_hits;
   vector<double>  *hlt1_muon_velo_hits;
   vector<double>  *hlt1_muon_velo_q;
   vector<double>  *hlt2_pvr_x;
   vector<double>  *hlt2_pvr_y;
   vector<double>  *hlt2_pvr_z;
   vector<double>  *hlt2_pvr_chi2;
   vector<double>  *hlt2_pvr_ndof;
   vector<double>  *hlt2_trk_idx_gen;
   vector<double>  *hlt2_trk_idx_sig;
   vector<double>  *hlt2_trk_idx_pvr;
   vector<double>  *hlt2_trk_dr_gen;
   vector<double>  *hlt2_trk_dr_sig;
   vector<double>  *hlt2_trk_charge;
   vector<double>  *hlt2_trk_type;
   vector<double>  *hlt2_trk_is_mu;
   vector<double>  *hlt2_trk_px;
   vector<double>  *hlt2_trk_py;
   vector<double>  *hlt2_trk_pz;
   vector<double>  *hlt2_trk_x;
   vector<double>  *hlt2_trk_y;
   vector<double>  *hlt2_trk_z;
   vector<double>  *hlt2_trk_ndof;
   vector<double>  *hlt2_trk_chi2;
   vector<double>  *hlt2_trk_prob_chi2;
   vector<double>  *hlt2_trk_hits_velo;
   vector<double>  *hlt2_trk_hits_tt;
   vector<double>  *hlt2_trk_hits_it;
   vector<double>  *hlt2_trk_hits_ot;
   vector<double>  *hlt2_trk_hits_muon;
   vector<double>  *hlt2_trk_ip_min;
   vector<double>  *hlt2_trk_ip_chi2_min;
   vector<double>  *hlt2_svr_idx_pvr;
   vector<double>  *hlt2_svr_idx_trk_1;
   vector<double>  *hlt2_svr_idx_trk_2;
   vector<double>  *hlt2_svr_idx_trk_3;
   vector<double>  *hlt2_svr_idx_trk_4;
   vector<double>  *hlt2_svr_x;
   vector<double>  *hlt2_svr_y;
   vector<double>  *hlt2_svr_z;
   vector<double>  *hlt2_svr_chi2;
   vector<double>  *hlt2_svr_ndof;
   vector<double>  *hlt2_svr_trk_doca_min;
   vector<double>  *hlt2_svr_ip_min;
   vector<double>  *hlt2_svr_ip_chi2_min;
   vector<double>  *hlt2_svr_fd;
   vector<double>  *hlt2_svr_fd_chi2;
   vector<double>  *hlt2_svr_px;
   vector<double>  *hlt2_svr_py;
   vector<double>  *hlt2_svr_pz;
   vector<double>  *hlt2_svr_m_cor;
   vector<double>  *hlt2_svr_ip_chi2_sum;
   vector<double>  *hlt2_svr_fdt_min;
   vector<double>  *hlt2_svr_n_trks;
   vector<double>  *hlt2_svr_bdt;
   vector<double>  *hlt2_svr_trg;
   vector<double>  *gen_idx_mother;
   vector<double>  *gen_idx_childf;
   vector<double>  *gen_idx_childl;
   vector<double>  *gen_pid;
   vector<double>  *gen_px;
   vector<double>  *gen_py;
   vector<double>  *gen_pz;
   vector<double>  *gen_e;
   vector<double>  *gen_x;
   vector<double>  *gen_y;
   vector<double>  *gen_z;

   // List of branches
   TBranch        *b_l0_spd_mult;   //!
   TBranch        *b_l0_elc_et;   //!
   TBranch        *b_l0_pho_et;   //!
   TBranch        *b_l0_had_et;   //!
   TBranch        *b_l0_mu1_pt;   //!
   TBranch        *b_l0_mu2_pt;   //!
   TBranch        *b_l0_dmu_pt;   //!
   TBranch        *b_hlt1_trk0_ip;   //!
   TBranch        *b_hlt1_trk0_ip_chi2;   //!
   TBranch        *b_hlt1_trk0_chi2_ndof;   //!
   TBranch        *b_hlt1_trk0_t_hits;   //!
   TBranch        *b_hlt1_trk0_velo_hits;   //!
   TBranch        *b_hlt1_trk0_velo_q;   //!
   TBranch        *b_hlt1_trk1_idx_gen;   //!
   TBranch        *b_hlt1_trk1_idx_sig;   //!
   TBranch        *b_hlt1_trk1_dr_gen;   //!
   TBranch        *b_hlt1_trk1_dr_sig;   //!
   TBranch        *b_hlt1_trk1_pt;   //!
   TBranch        *b_hlt1_trk1_p;   //!
   TBranch        *b_hlt1_trk1_ip;   //!
   TBranch        *b_hlt1_trk1_ip_chi2;   //!
   TBranch        *b_hlt1_trk1_chi2_ndof;   //!
   TBranch        *b_hlt1_trk1_t_hits;   //!
   TBranch        *b_hlt1_trk1_velo_hits;   //!
   TBranch        *b_hlt1_trk1_velo_q;   //!
   TBranch        *b_hlt1_muon_idx_gen;   //!
   TBranch        *b_hlt1_muon_idx_sig;   //!
   TBranch        *b_hlt1_muon_dr_gen;   //!
   TBranch        *b_hlt1_muon_dr_sig;   //!
   TBranch        *b_hlt1_muon_pt;   //!
   TBranch        *b_hlt1_muon_p;   //!
   TBranch        *b_hlt1_muon_ip;   //!
   TBranch        *b_hlt1_muon_ip_chi2;   //!
   TBranch        *b_hlt1_muon_chi2_ndof;   //!
   TBranch        *b_hlt1_muon_t_hits;   //!
   TBranch        *b_hlt1_muon_velo_hits;   //!
   TBranch        *b_hlt1_muon_velo_q;   //!
   TBranch        *b_hlt2_pvr_x;   //!
   TBranch        *b_hlt2_pvr_y;   //!
   TBranch        *b_hlt2_pvr_z;   //!
   TBranch        *b_hlt2_pvr_chi2;   //!
   TBranch        *b_hlt2_pvr_ndof;   //!
   TBranch        *b_hlt2_trk_idx_gen;   //!
   TBranch        *b_hlt2_trk_idx_sig;   //!
   TBranch        *b_hlt2_trk_idx_pvr;   //!
   TBranch        *b_hlt2_trk_dr_gen;   //!
   TBranch        *b_hlt2_trk_dr_sig;   //!
   TBranch        *b_hlt2_trk_charge;   //!
   TBranch        *b_hlt2_trk_type;   //!
   TBranch        *b_hlt2_trk_is_mu;   //!
   TBranch        *b_hlt2_trk_px;   //!
   TBranch        *b_hlt2_trk_py;   //!
   TBranch        *b_hlt2_trk_pz;   //!
   TBranch        *b_hlt2_trk_x;   //!
   TBranch        *b_hlt2_trk_y;   //!
   TBranch        *b_hlt2_trk_z;   //!
   TBranch        *b_hlt2_trk_ndof;   //!
   TBranch        *b_hlt2_trk_chi2;   //!
   TBranch        *b_hlt2_trk_prob_chi2;   //!
   TBranch        *b_hlt2_trk_hits_velo;   //!
   TBranch        *b_hlt2_trk_hits_tt;   //!
   TBranch        *b_hlt2_trk_hits_it;   //!
   TBranch        *b_hlt2_trk_hits_ot;   //!
   TBranch        *b_hlt2_trk_hits_muon;   //!
   TBranch        *b_hlt2_trk_ip_min;   //!
   TBranch        *b_hlt2_trk_ip_chi2_min;   //!
   TBranch        *b_hlt2_svr_idx_pvr;   //!
   TBranch        *b_hlt2_svr_idx_trk_1;   //!
   TBranch        *b_hlt2_svr_idx_trk_2;   //!
   TBranch        *b_hlt2_svr_idx_trk_3;   //!
   TBranch        *b_hlt2_svr_idx_trk_4;   //!
   TBranch        *b_hlt2_svr_x;   //!
   TBranch        *b_hlt2_svr_y;   //!
   TBranch        *b_hlt2_svr_z;   //!
   TBranch        *b_hlt2_svr_chi2;   //!
   TBranch        *b_hlt2_svr_ndof;   //!
   TBranch        *b_hlt2_svr_trk_doca_min;   //!
   TBranch        *b_hlt2_svr_ip_min;   //!
   TBranch        *b_hlt2_svr_ip_chi2_min;   //!
   TBranch        *b_hlt2_svr_fd;   //!
   TBranch        *b_hlt2_svr_fd_chi2;   //!
   TBranch        *b_hlt2_svr_px;   //!
   TBranch        *b_hlt2_svr_py;   //!
   TBranch        *b_hlt2_svr_pz;   //!
   TBranch        *b_hlt2_svr_m_cor;   //!
   TBranch        *b_hlt2_svr_ip_chi2_sum;   //!
   TBranch        *b_hlt2_svr_fdt_min;   //!
   TBranch        *b_hlt2_svr_n_trks;   //!
   TBranch        *b_hlt2_svr_bdt;   //!
   TBranch        *b_hlt2_svr_trg;   //!
   TBranch        *b_gen_idx_mother;   //!
   TBranch        *b_gen_idx_childf;   //!
   TBranch        *b_gen_idx_childl;   //!
   TBranch        *b_gen_pid;   //!
   TBranch        *b_gen_px;   //!
   TBranch        *b_gen_py;   //!
   TBranch        *b_gen_pz;   //!
   TBranch        *b_gen_e;   //!
   TBranch        *b_gen_x;   //!
   TBranch        *b_gen_y;   //!
   TBranch        *b_gen_z;   //!

   data(TTree *tree=0);
   virtual ~data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef data_cxx
data::data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Hlt.MC15.MD.11114001.1.150204.00.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Hlt.MC15.MD.11114001.1.150204.00.root");
      }
      f->GetObject("data",tree);

   }
*/
   Init(tree);
}

data::~data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hlt1_trk0_ip = 0;
   hlt1_trk0_ip_chi2 = 0;
   hlt1_trk0_chi2_ndof = 0;
   hlt1_trk0_t_hits = 0;
   hlt1_trk0_velo_hits = 0;
   hlt1_trk0_velo_q = 0;
   hlt1_trk1_idx_gen = 0;
   hlt1_trk1_idx_sig = 0;
   hlt1_trk1_dr_gen = 0;
   hlt1_trk1_dr_sig = 0;
   hlt1_trk1_pt = 0;
   hlt1_trk1_p = 0;
   hlt1_trk1_ip = 0;
   hlt1_trk1_ip_chi2 = 0;
   hlt1_trk1_chi2_ndof = 0;
   hlt1_trk1_t_hits = 0;
   hlt1_trk1_velo_hits = 0;
   hlt1_trk1_velo_q = 0;
   hlt1_muon_idx_gen = 0;
   hlt1_muon_idx_sig = 0;
   hlt1_muon_dr_gen = 0;
   hlt1_muon_dr_sig = 0;
   hlt1_muon_pt = 0;
   hlt1_muon_p = 0;
   hlt1_muon_ip = 0;
   hlt1_muon_ip_chi2 = 0;
   hlt1_muon_chi2_ndof = 0;
   hlt1_muon_t_hits = 0;
   hlt1_muon_velo_hits = 0;
   hlt1_muon_velo_q = 0;
   hlt2_pvr_x = 0;
   hlt2_pvr_y = 0;
   hlt2_pvr_z = 0;
   hlt2_pvr_chi2 = 0;
   hlt2_pvr_ndof = 0;
   hlt2_trk_idx_gen = 0;
   hlt2_trk_idx_sig = 0;
   hlt2_trk_idx_pvr = 0;
   hlt2_trk_dr_gen = 0;
   hlt2_trk_dr_sig = 0;
   hlt2_trk_charge = 0;
   hlt2_trk_type = 0;
   hlt2_trk_is_mu = 0;
   hlt2_trk_px = 0;
   hlt2_trk_py = 0;
   hlt2_trk_pz = 0;
   hlt2_trk_x = 0;
   hlt2_trk_y = 0;
   hlt2_trk_z = 0;
   hlt2_trk_ndof = 0;
   hlt2_trk_chi2 = 0;
   hlt2_trk_prob_chi2 = 0;
   hlt2_trk_hits_velo = 0;
   hlt2_trk_hits_tt = 0;
   hlt2_trk_hits_it = 0;
   hlt2_trk_hits_ot = 0;
   hlt2_trk_hits_muon = 0;
   hlt2_trk_ip_min = 0;
   hlt2_trk_ip_chi2_min = 0;
   hlt2_svr_idx_pvr = 0;
   hlt2_svr_idx_trk_1 = 0;
   hlt2_svr_idx_trk_2 = 0;
   hlt2_svr_idx_trk_3 = 0;
   hlt2_svr_idx_trk_4 = 0;
   hlt2_svr_x = 0;
   hlt2_svr_y = 0;
   hlt2_svr_z = 0;
   hlt2_svr_chi2 = 0;
   hlt2_svr_ndof = 0;
   hlt2_svr_trk_doca_min = 0;
   hlt2_svr_ip_min = 0;
   hlt2_svr_ip_chi2_min = 0;
   hlt2_svr_fd = 0;
   hlt2_svr_fd_chi2 = 0;
   hlt2_svr_px = 0;
   hlt2_svr_py = 0;
   hlt2_svr_pz = 0;
   hlt2_svr_m_cor = 0;
   hlt2_svr_ip_chi2_sum = 0;
   hlt2_svr_fdt_min = 0;
   hlt2_svr_n_trks = 0;
   hlt2_svr_bdt = 0;
   hlt2_svr_trg = 0;
   gen_idx_mother = 0;
   gen_idx_childf = 0;
   gen_idx_childl = 0;
   gen_pid = 0;
   gen_px = 0;
   gen_py = 0;
   gen_pz = 0;
   gen_e = 0;
   gen_x = 0;
   gen_y = 0;
   gen_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("l0_spd_mult", &l0_spd_mult, &b_l0_spd_mult);
   fChain->SetBranchAddress("l0_elc_et", &l0_elc_et, &b_l0_elc_et);
   fChain->SetBranchAddress("l0_pho_et", &l0_pho_et, &b_l0_pho_et);
   fChain->SetBranchAddress("l0_had_et", &l0_had_et, &b_l0_had_et);
   fChain->SetBranchAddress("l0_mu1_pt", &l0_mu1_pt, &b_l0_mu1_pt);
   fChain->SetBranchAddress("l0_mu2_pt", &l0_mu2_pt, &b_l0_mu2_pt);
   fChain->SetBranchAddress("l0_dmu_pt", &l0_dmu_pt, &b_l0_dmu_pt);
   fChain->SetBranchAddress("hlt1_trk0_ip", &hlt1_trk0_ip, &b_hlt1_trk0_ip);
   fChain->SetBranchAddress("hlt1_trk0_ip_chi2", &hlt1_trk0_ip_chi2, &b_hlt1_trk0_ip_chi2);
   fChain->SetBranchAddress("hlt1_trk0_chi2_ndof", &hlt1_trk0_chi2_ndof, &b_hlt1_trk0_chi2_ndof);
   fChain->SetBranchAddress("hlt1_trk0_t_hits", &hlt1_trk0_t_hits, &b_hlt1_trk0_t_hits);
   fChain->SetBranchAddress("hlt1_trk0_velo_hits", &hlt1_trk0_velo_hits, &b_hlt1_trk0_velo_hits);
   fChain->SetBranchAddress("hlt1_trk0_velo_q", &hlt1_trk0_velo_q, &b_hlt1_trk0_velo_q);
   fChain->SetBranchAddress("hlt1_trk1_idx_gen", &hlt1_trk1_idx_gen, &b_hlt1_trk1_idx_gen);
   fChain->SetBranchAddress("hlt1_trk1_idx_sig", &hlt1_trk1_idx_sig, &b_hlt1_trk1_idx_sig);
   fChain->SetBranchAddress("hlt1_trk1_dr_gen", &hlt1_trk1_dr_gen, &b_hlt1_trk1_dr_gen);
   fChain->SetBranchAddress("hlt1_trk1_dr_sig", &hlt1_trk1_dr_sig, &b_hlt1_trk1_dr_sig);
   fChain->SetBranchAddress("hlt1_trk1_pt", &hlt1_trk1_pt, &b_hlt1_trk1_pt);
   fChain->SetBranchAddress("hlt1_trk1_p", &hlt1_trk1_p, &b_hlt1_trk1_p);
   fChain->SetBranchAddress("hlt1_trk1_ip", &hlt1_trk1_ip, &b_hlt1_trk1_ip);
   fChain->SetBranchAddress("hlt1_trk1_ip_chi2", &hlt1_trk1_ip_chi2, &b_hlt1_trk1_ip_chi2);
   fChain->SetBranchAddress("hlt1_trk1_chi2_ndof", &hlt1_trk1_chi2_ndof, &b_hlt1_trk1_chi2_ndof);
   fChain->SetBranchAddress("hlt1_trk1_t_hits", &hlt1_trk1_t_hits, &b_hlt1_trk1_t_hits);
   fChain->SetBranchAddress("hlt1_trk1_velo_hits", &hlt1_trk1_velo_hits, &b_hlt1_trk1_velo_hits);
   fChain->SetBranchAddress("hlt1_trk1_velo_q", &hlt1_trk1_velo_q, &b_hlt1_trk1_velo_q);
   fChain->SetBranchAddress("hlt1_muon_idx_gen", &hlt1_muon_idx_gen, &b_hlt1_muon_idx_gen);
   fChain->SetBranchAddress("hlt1_muon_idx_sig", &hlt1_muon_idx_sig, &b_hlt1_muon_idx_sig);
   fChain->SetBranchAddress("hlt1_muon_dr_gen", &hlt1_muon_dr_gen, &b_hlt1_muon_dr_gen);
   fChain->SetBranchAddress("hlt1_muon_dr_sig", &hlt1_muon_dr_sig, &b_hlt1_muon_dr_sig);
   fChain->SetBranchAddress("hlt1_muon_pt", &hlt1_muon_pt, &b_hlt1_muon_pt);
   fChain->SetBranchAddress("hlt1_muon_p", &hlt1_muon_p, &b_hlt1_muon_p);
   fChain->SetBranchAddress("hlt1_muon_ip", &hlt1_muon_ip, &b_hlt1_muon_ip);
   fChain->SetBranchAddress("hlt1_muon_ip_chi2", &hlt1_muon_ip_chi2, &b_hlt1_muon_ip_chi2);
   fChain->SetBranchAddress("hlt1_muon_chi2_ndof", &hlt1_muon_chi2_ndof, &b_hlt1_muon_chi2_ndof);
   fChain->SetBranchAddress("hlt1_muon_t_hits", &hlt1_muon_t_hits, &b_hlt1_muon_t_hits);
   fChain->SetBranchAddress("hlt1_muon_velo_hits", &hlt1_muon_velo_hits, &b_hlt1_muon_velo_hits);
   fChain->SetBranchAddress("hlt1_muon_velo_q", &hlt1_muon_velo_q, &b_hlt1_muon_velo_q);
   fChain->SetBranchAddress("hlt2_pvr_x", &hlt2_pvr_x, &b_hlt2_pvr_x);
   fChain->SetBranchAddress("hlt2_pvr_y", &hlt2_pvr_y, &b_hlt2_pvr_y);
   fChain->SetBranchAddress("hlt2_pvr_z", &hlt2_pvr_z, &b_hlt2_pvr_z);
   fChain->SetBranchAddress("hlt2_pvr_chi2", &hlt2_pvr_chi2, &b_hlt2_pvr_chi2);
   fChain->SetBranchAddress("hlt2_pvr_ndof", &hlt2_pvr_ndof, &b_hlt2_pvr_ndof);
   fChain->SetBranchAddress("hlt2_trk_idx_gen", &hlt2_trk_idx_gen, &b_hlt2_trk_idx_gen);
   fChain->SetBranchAddress("hlt2_trk_idx_sig", &hlt2_trk_idx_sig, &b_hlt2_trk_idx_sig);
   fChain->SetBranchAddress("hlt2_trk_idx_pvr", &hlt2_trk_idx_pvr, &b_hlt2_trk_idx_pvr);
   fChain->SetBranchAddress("hlt2_trk_dr_gen", &hlt2_trk_dr_gen, &b_hlt2_trk_dr_gen);
   fChain->SetBranchAddress("hlt2_trk_dr_sig", &hlt2_trk_dr_sig, &b_hlt2_trk_dr_sig);
   fChain->SetBranchAddress("hlt2_trk_charge", &hlt2_trk_charge, &b_hlt2_trk_charge);
   fChain->SetBranchAddress("hlt2_trk_type", &hlt2_trk_type, &b_hlt2_trk_type);
   fChain->SetBranchAddress("hlt2_trk_is_mu", &hlt2_trk_is_mu, &b_hlt2_trk_is_mu);
   fChain->SetBranchAddress("hlt2_trk_px", &hlt2_trk_px, &b_hlt2_trk_px);
   fChain->SetBranchAddress("hlt2_trk_py", &hlt2_trk_py, &b_hlt2_trk_py);
   fChain->SetBranchAddress("hlt2_trk_pz", &hlt2_trk_pz, &b_hlt2_trk_pz);
   fChain->SetBranchAddress("hlt2_trk_x", &hlt2_trk_x, &b_hlt2_trk_x);
   fChain->SetBranchAddress("hlt2_trk_y", &hlt2_trk_y, &b_hlt2_trk_y);
   fChain->SetBranchAddress("hlt2_trk_z", &hlt2_trk_z, &b_hlt2_trk_z);
   fChain->SetBranchAddress("hlt2_trk_ndof", &hlt2_trk_ndof, &b_hlt2_trk_ndof);
   fChain->SetBranchAddress("hlt2_trk_chi2", &hlt2_trk_chi2, &b_hlt2_trk_chi2);
   fChain->SetBranchAddress("hlt2_trk_prob_chi2", &hlt2_trk_prob_chi2, &b_hlt2_trk_prob_chi2);
   fChain->SetBranchAddress("hlt2_trk_hits_velo", &hlt2_trk_hits_velo, &b_hlt2_trk_hits_velo);
   fChain->SetBranchAddress("hlt2_trk_hits_tt", &hlt2_trk_hits_tt, &b_hlt2_trk_hits_tt);
   fChain->SetBranchAddress("hlt2_trk_hits_it", &hlt2_trk_hits_it, &b_hlt2_trk_hits_it);
   fChain->SetBranchAddress("hlt2_trk_hits_ot", &hlt2_trk_hits_ot, &b_hlt2_trk_hits_ot);
   fChain->SetBranchAddress("hlt2_trk_hits_muon", &hlt2_trk_hits_muon, &b_hlt2_trk_hits_muon);
   fChain->SetBranchAddress("hlt2_trk_ip_min", &hlt2_trk_ip_min, &b_hlt2_trk_ip_min);
   fChain->SetBranchAddress("hlt2_trk_ip_chi2_min", &hlt2_trk_ip_chi2_min, &b_hlt2_trk_ip_chi2_min);
   fChain->SetBranchAddress("hlt2_svr_idx_pvr", &hlt2_svr_idx_pvr, &b_hlt2_svr_idx_pvr);
   fChain->SetBranchAddress("hlt2_svr_idx_trk_1", &hlt2_svr_idx_trk_1, &b_hlt2_svr_idx_trk_1);
   fChain->SetBranchAddress("hlt2_svr_idx_trk_2", &hlt2_svr_idx_trk_2, &b_hlt2_svr_idx_trk_2);
   fChain->SetBranchAddress("hlt2_svr_idx_trk_3", &hlt2_svr_idx_trk_3, &b_hlt2_svr_idx_trk_3);
   fChain->SetBranchAddress("hlt2_svr_idx_trk_4", &hlt2_svr_idx_trk_4, &b_hlt2_svr_idx_trk_4);
   fChain->SetBranchAddress("hlt2_svr_x", &hlt2_svr_x, &b_hlt2_svr_x);
   fChain->SetBranchAddress("hlt2_svr_y", &hlt2_svr_y, &b_hlt2_svr_y);
   fChain->SetBranchAddress("hlt2_svr_z", &hlt2_svr_z, &b_hlt2_svr_z);
   fChain->SetBranchAddress("hlt2_svr_chi2", &hlt2_svr_chi2, &b_hlt2_svr_chi2);
   fChain->SetBranchAddress("hlt2_svr_ndof", &hlt2_svr_ndof, &b_hlt2_svr_ndof);
   fChain->SetBranchAddress("hlt2_svr_trk_doca_min", &hlt2_svr_trk_doca_min, &b_hlt2_svr_trk_doca_min);
   fChain->SetBranchAddress("hlt2_svr_ip_min", &hlt2_svr_ip_min, &b_hlt2_svr_ip_min);
   fChain->SetBranchAddress("hlt2_svr_ip_chi2_min", &hlt2_svr_ip_chi2_min, &b_hlt2_svr_ip_chi2_min);
   fChain->SetBranchAddress("hlt2_svr_fd", &hlt2_svr_fd, &b_hlt2_svr_fd);
   fChain->SetBranchAddress("hlt2_svr_fd_chi2", &hlt2_svr_fd_chi2, &b_hlt2_svr_fd_chi2);
   fChain->SetBranchAddress("hlt2_svr_px", &hlt2_svr_px, &b_hlt2_svr_px);
   fChain->SetBranchAddress("hlt2_svr_py", &hlt2_svr_py, &b_hlt2_svr_py);
   fChain->SetBranchAddress("hlt2_svr_pz", &hlt2_svr_pz, &b_hlt2_svr_pz);
   fChain->SetBranchAddress("hlt2_svr_m_cor", &hlt2_svr_m_cor, &b_hlt2_svr_m_cor);
   fChain->SetBranchAddress("hlt2_svr_ip_chi2_sum", &hlt2_svr_ip_chi2_sum, &b_hlt2_svr_ip_chi2_sum);
   fChain->SetBranchAddress("hlt2_svr_fdt_min", &hlt2_svr_fdt_min, &b_hlt2_svr_fdt_min);
   fChain->SetBranchAddress("hlt2_svr_n_trks", &hlt2_svr_n_trks, &b_hlt2_svr_n_trks);
   fChain->SetBranchAddress("hlt2_svr_bdt", &hlt2_svr_bdt, &b_hlt2_svr_bdt);
   fChain->SetBranchAddress("hlt2_svr_trg", &hlt2_svr_trg, &b_hlt2_svr_trg);
   fChain->SetBranchAddress("gen_idx_mother", &gen_idx_mother, &b_gen_idx_mother);
   fChain->SetBranchAddress("gen_idx_childf", &gen_idx_childf, &b_gen_idx_childf);
   fChain->SetBranchAddress("gen_idx_childl", &gen_idx_childl, &b_gen_idx_childl);
   fChain->SetBranchAddress("gen_pid", &gen_pid, &b_gen_pid);
   fChain->SetBranchAddress("gen_px", &gen_px, &b_gen_px);
   fChain->SetBranchAddress("gen_py", &gen_py, &b_gen_py);
   fChain->SetBranchAddress("gen_pz", &gen_pz, &b_gen_pz);
   fChain->SetBranchAddress("gen_e", &gen_e, &b_gen_e);
   fChain->SetBranchAddress("gen_x", &gen_x, &b_gen_x);
   fChain->SetBranchAddress("gen_y", &gen_y, &b_gen_y);
   fChain->SetBranchAddress("gen_z", &gen_z, &b_gen_z);
   Notify();
}

Bool_t data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef data_cxx
