//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  9 11:46:32 2016 by ROOT version 5.34/34
// from TChain dune_dst/
//////////////////////////////////////////////////////////

#ifndef make_effs_h
#define make_effs_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <vector>
#include <TString.h>
#include <TEfficiency.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxtrueFSParticles = 100;
const Int_t kMaxrecoFSParticles = 100;

class dune_dst {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           SubRun;
   Int_t           Event;
   Int_t           Detector;

   // Truth variables
   Float_t         Ev;
   Double_t        POTWeight;
   Int_t           sample;
   Int_t           reactionmode;
   Float_t         trueVtxX;
   Float_t         trueVtxY;
   Float_t         trueVtxZ;
   Float_t         trueDistToEdge;
   Float_t         y_true;
   Int_t           target;
   Int_t           NuPDG;
   Bool_t          NuInput;
   Int_t           trueFSParticles;
   Int_t           trueFSParticles_pdg[kMaxtrueFSParticles];   //[trueFSParticles_]
   Int_t           trueFSParticles_id[kMaxtrueFSParticles];   //[trueFSParticles_]
   Float_t         trueFSParticles_energy[kMaxtrueFSParticles];   //[trueFSParticles_]
   Int_t           trueFSParticles_numHits[kMaxtrueFSParticles];   //[trueFSParticles_]
   float   trueFSParticles_dEdx[kMaxtrueFSParticles]; // not used
   Float_t         trueFSParticles_momentum[kMaxtrueFSParticles][3];   //[trueFSParticles_]

   // Reco variables
   Int_t           NuPDG_reco;
   Int_t           sample_reco;
   Int_t           reaction_reco;
   Float_t         Ev_reco;
   Float_t         y_reco;
   Float_t         recoVtxX;
   Float_t         recoVtxY;
   Float_t         recoVtxZ;
   Float_t         sharedE_reco;
   Int_t           recoFSParticles;
   Int_t           recoFSParticles_pdg[kMaxrecoFSParticles];   //[recoFSParticles_]
   Int_t           recoFSParticles_id[kMaxrecoFSParticles];   //[recoFSParticles_]
   Float_t         recoFSParticles_energy[kMaxrecoFSParticles];   //[recoFSParticles_]
   Int_t           recoFSParticles_numHits[kMaxrecoFSParticles];   //[recoFSParticles_]
   float   recoFSParticles_dEdx[kMaxrecoFSParticles]; // not used
   Float_t         recoFSParticles_momentum[kMaxrecoFSParticles][3];   //[recoFSParticles_]

   // List of branches
   //TBranch        *b_genie_mc_truth;   //!
   TBranch        *b_ana_Run;   //!
   TBranch        *b_ana_SubRun;   //!
   TBranch        *b_ana_Event;   //!
   TBranch        *b_ana_Detector;   //!
   TBranch        *b_ana_Ev;   //!
   TBranch        *b_ana_POTWeight;   //!
   TBranch        *b_ana_sample;   //!
   TBranch        *b_ana_reactionmode;   //!
   TBranch        *b_ana_trueVtxX;   //!
   TBranch        *b_ana_trueVtxY;   //!
   TBranch        *b_ana_trueVtxZ;   //!
   TBranch        *b_ana_trueDistToEdge;   //!
   TBranch        *b_ana_y_true;   //!
   TBranch        *b_ana_target;   //!
   TBranch        *b_ana_NuPDG;   //!
   TBranch        *b_ana_NuInput;   //!
   TBranch        *b_ana_trueFSParticles_;   //!
   TBranch        *b_trueFSParticles_pdg;   //!
   TBranch        *b_trueFSParticles_id;   //!
   TBranch        *b_trueFSParticles_energy;   //!
   TBranch        *b_trueFSParticles_numHits;   //!
   TBranch        *b_trueFSParticles_dEdx;   //!
   TBranch        *b_trueFSParticles_momentum;   //!
   TBranch        *b_ana_NuPDG_reco;   //!
   TBranch        *b_ana_sample_reco;   //!
   TBranch        *b_ana_reaction_reco;   //!
   TBranch        *b_ana_Ev_reco;   //!
   TBranch        *b_ana_y_reco;   //!
   TBranch        *b_ana_recoVtxX;   //!
   TBranch        *b_ana_recoVtxY;   //!
   TBranch        *b_ana_recoVtxZ;   //!
   TBranch        *b_ana_sharedE_reco;   //!
   TBranch        *b_ana_recoFSParticles_;   //!
   TBranch        *b_recoFSParticles_pdg;   //!
   TBranch        *b_recoFSParticles_id;   //!
   TBranch        *b_recoFSParticles_energy;   //!
   TBranch        *b_recoFSParticles_numHits;   //!
   TBranch        *b_recoFSParticles_dEdx;   //!
   TBranch        *b_recoFSParticles_momentum;   //!

   dune_dst(TTree *tree=0);
   virtual ~dune_dst();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int n_evt, char* tag,char * fOutFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef make_effs_cxx
dune_dst::dune_dst(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("dune_dst","");

      //FGT files  
       chain->Add("../DUNE_ND_ntps/full/fgt/ndtf_output_nu_fgt.dst.root");
/*      for ( int i = 1; i < 1778; i++) {
	chain->Add(Form("DUNE_ND_ntps/fgt/nu/ndtf_output_nu_fgt_%d.dst.root", i));      
      }
*/      
      //GAr files
      /* chain->Add("../DUNE_ND_ntps/full/gar/ndtf_output_nu_gar.dst.root"); */

      //LAr files
     // chain->Add("../DUNE_ND_ntps/full/lar/ndtf_output_nu_lar.dst.root");
      
      tree = chain;

   }
   Init(tree);
}

dune_dst::~dune_dst()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dune_dst::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dune_dst::LoadTree(Long64_t entry)
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

void dune_dst::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   //genie_mc_truth = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   //fChain->SetBranchAddress("genie_mc_truth", &genie_mc_truth, &b_genie_mc_truth);
   fChain->SetBranchAddress("Run", &Run, &b_ana_Run);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_ana_SubRun);
   fChain->SetBranchAddress("Event", &Event, &b_ana_Event);
   fChain->SetBranchAddress("Detector", &Detector, &b_ana_Detector);
   fChain->SetBranchAddress("Ev", &Ev, &b_ana_Ev);
   fChain->SetBranchAddress("POTWeight", &POTWeight, &b_ana_POTWeight);
   fChain->SetBranchAddress("sample", &sample, &b_ana_sample);
   fChain->SetBranchAddress("reactionmode", &reactionmode, &b_ana_reactionmode);
   fChain->SetBranchAddress("trueVtxX", &trueVtxX, &b_ana_trueVtxX);
   fChain->SetBranchAddress("trueVtxY", &trueVtxY, &b_ana_trueVtxY);
   fChain->SetBranchAddress("trueVtxZ", &trueVtxZ, &b_ana_trueVtxZ);
   fChain->SetBranchAddress("trueDistToEdge", &trueDistToEdge, &b_ana_trueDistToEdge);
   fChain->SetBranchAddress("y_true", &y_true, &b_ana_y_true);
   fChain->SetBranchAddress("target", &target, &b_ana_target);
   fChain->SetBranchAddress("NuPDG", &NuPDG, &b_ana_NuPDG);
   fChain->SetBranchAddress("NuInput", &NuInput, &b_ana_NuInput);
   fChain->SetBranchAddress("trueFSParticles", &trueFSParticles, &b_ana_trueFSParticles_);
   fChain->SetBranchAddress("trueFSParticles.pdg", trueFSParticles_pdg, &b_trueFSParticles_pdg);
   fChain->SetBranchAddress("trueFSParticles.id", trueFSParticles_id, &b_trueFSParticles_id);
   fChain->SetBranchAddress("trueFSParticles.energy", trueFSParticles_energy, &b_trueFSParticles_energy);
   fChain->SetBranchAddress("trueFSParticles.numHits", trueFSParticles_numHits, &b_trueFSParticles_numHits);
   fChain->SetBranchAddress("trueFSParticles.dEdx", trueFSParticles_dEdx, &b_trueFSParticles_dEdx);
   fChain->SetBranchAddress("trueFSParticles.momentum[3]", trueFSParticles_momentum, &b_trueFSParticles_momentum);
   fChain->SetBranchAddress("NuPDG_reco", &NuPDG_reco, &b_ana_NuPDG_reco);
   fChain->SetBranchAddress("sample_reco", &sample_reco, &b_ana_sample_reco);
   fChain->SetBranchAddress("reaction_reco", &reaction_reco, &b_ana_reaction_reco);
   fChain->SetBranchAddress("Ev_reco", &Ev_reco, &b_ana_Ev_reco);
   fChain->SetBranchAddress("y_reco", &y_reco, &b_ana_y_reco);
   fChain->SetBranchAddress("recoVtxX", &recoVtxX, &b_ana_recoVtxX);
   fChain->SetBranchAddress("recoVtxY", &recoVtxY, &b_ana_recoVtxY);
   fChain->SetBranchAddress("recoVtxZ", &recoVtxZ, &b_ana_recoVtxZ);
   fChain->SetBranchAddress("sharedE_reco", &sharedE_reco, &b_ana_sharedE_reco);
   fChain->SetBranchAddress("recoFSParticles", &recoFSParticles, &b_ana_recoFSParticles_);
   fChain->SetBranchAddress("recoFSParticles.pdg", recoFSParticles_pdg, &b_recoFSParticles_pdg);
   fChain->SetBranchAddress("recoFSParticles.id", recoFSParticles_id, &b_recoFSParticles_id);
   fChain->SetBranchAddress("recoFSParticles.energy", recoFSParticles_energy, &b_recoFSParticles_energy);
   fChain->SetBranchAddress("recoFSParticles.numHits", recoFSParticles_numHits, &b_recoFSParticles_numHits);
   fChain->SetBranchAddress("recoFSParticles.dEdx", recoFSParticles_dEdx, &b_recoFSParticles_dEdx);
   fChain->SetBranchAddress("recoFSParticles.momentum[3]", recoFSParticles_momentum, &b_recoFSParticles_momentum);
   Notify();
}

Bool_t dune_dst::Notify()
{
   // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dune_dst::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dune_dst::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef dune_dst_cxx



void doEff( TH2D *h[3] )
{

  double total = h[1]->Integral();
  double lowEff[11] = {0.};
  for( int x = 1; x <= h[0]->GetNbinsX(); ++x ) {
    for( int y = 1; y <= h[0]->GetNbinsY(); ++y ) {
      double num = h[0]->GetBinContent( x, y );
      double denom = h[1]->GetBinContent( x, y );

      if( denom > 0. ){
        
        if ( num == 0. )h[2]->SetBinContent( x, y, -0.0099  );
        else h[2]->SetBinContent( x, y, num/denom );

      }
      else h[2]->SetBinContent( x, y, -1. );

    }
  }

  h[2]->SetMinimum(-0.01); // don't draw bins with no true events, but do draw bins with no acceptance
  h[2]->SetMaximum(1.);
}
const float pi = 4.*atan(1.);


///////p costheta
TH2D * muon_p_costheta_trk[3]; 
TH2D * piplus_p_costheta_trk[3];
TH2D * piminus_p_costheta_trk[3];
TH2D * pi0_p_costheta_trk[3];
TH2D * proton_p_costheta_trk[3];

TH2D * muon_p_costheta_pid[3]; 
TH2D * piplus_p_costheta_pid[3];
TH2D * piminus_p_costheta_pid[3];
TH2D * pi0_p_costheta_pid[3];
TH2D * proton_p_costheta_pid[3];

///////costheta
TH1D * muon_costheta_pid[3]; 
TH1D * piplus_costheta_pid[3];
TH1D * piminus_costheta_pid[3];
TH1D * pi0_costheta_pid[3];
TH1D * proton_costheta_pid[3];

TH1D * muon_costheta_trk[3]; 
TH1D * piplus_costheta_trk[3];
TH1D * piminus_costheta_trk[3];
TH1D * pi0_costheta_trk[3];
TH1D * proton_costheta_trk[3];

///////p
TH1D * muon_p_trk[3]; 
TH1D * piplus_p_trk[3];
TH1D * piminus_p_trk[3];
TH1D * pi0_p_trk[3];
TH1D * proton_p_trk[3];

TH1D * muon_p_pid[3]; 
TH1D * piplus_p_pid[3];
TH1D * piminus_p_pid[3];
TH1D * pi0_p_pid[3];
TH1D * proton_p_pid[3];

//////phi
TH1D * muon_phi_trk[3]; 
TH1D * piplus_phi_trk[3];
TH1D * piminus_phi_trk[3];
TH1D * pi0_phi_trk[3];
TH1D * proton_phi_trk[3];

TH1D * muon_phi_pid[3]; 
TH1D * piplus_phi_pid[3];
TH1D * piminus_phi_pid[3];
TH1D * pi0_phi_pid[3];
TH1D * proton_phi_pid[3];
 
//////phi costheta
TH2D * muon_phi_costheta_trk[3]; 
TH2D * piplus_phi_costheta_trk[3];
TH2D * piminus_phi_costheta_trk[3];
TH2D * pi0_phi_costheta_trk[3];
TH2D * proton_phi_costheta_trk[3];

TH2D * muon_phi_costheta_pid[3]; 
TH2D * piplus_phi_costheta_pid[3];
TH2D * piminus_phi_costheta_pid[3];
TH2D * pi0_phi_costheta_pid[3];
TH2D * proton_phi_costheta_pid[3];

//////phi p
TH2D * muon_phi_p_trk[3]; 
TH2D * piplus_phi_p_trk[3];
TH2D * piminus_phi_p_trk[3];
TH2D * pi0_phi_p_trk[3];
TH2D * proton_phi_p_trk[3];

TH2D * muon_phi_p_pid[3]; 
TH2D * piplus_phi_p_pid[3];
TH2D * piminus_phi_p_pid[3];
TH2D * pi0_phi_p_pid[3];
TH2D * proton_phi_p_pid[3];

//////primary
TH1D * total_int = new TH1D("total", "", 100, 0, 10);
TH1D * argon_int = new TH1D("argon", "", 100, 0, 10);
TH1D * reco_int = new TH1D("reco", "", 100, 0, 10);
TH1D * pid_int = new TH1D("pid", "", 100, 0, 10);

/////Smear
TH2D * had_smear = new TH2D("had_smear",";True E_{had} (GeV);Reconstructed E_{had} (GeV)",100,0,10,100,0,10);
TH2D * muon_smear = new TH2D("muon_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);
TH2D * proton_HM_smear = new TH2D("HM_proton_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);
TH2D * piplus_HM_smear = new TH2D("HM_piplus_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);
TH2D * piminus_HM_smear = new TH2D("HM_piminus_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);
TH2D * proton_SHM_smear = new TH2D("SHM_proton_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);
TH2D * piplus_SHM_smear = new TH2D("SHM_piplus_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);
TH2D * piminus_SHM_smear = new TH2D("SHM_piminus_smear",";True E_{#mu} (GeV);Reconstructed E_{#mu} (GeV)",100,0,10,100,0,10);




void fillHists(char* type, int pdg, int ind, float p, float costheta, float phi) {
    if(pdg == 2212){
      if (type == "trk"){
	proton_p_costheta_trk[ind]->Fill(p, costheta);
	proton_phi_costheta_trk[ind]->Fill(phi, costheta);
	proton_p_trk[ind]->Fill(p);
	proton_costheta_trk[ind]->Fill(costheta);
	proton_phi_p_trk[ind]->Fill(phi, p);
	proton_phi_trk[ind]->Fill(phi);
      }
      else if (type == "PID") {
	proton_p_costheta_pid[ind]->Fill(p, costheta);
	proton_phi_costheta_pid[ind]->Fill(phi, costheta);
	proton_p_pid[ind]->Fill(p);
	proton_costheta_pid[ind]->Fill(costheta);
	proton_phi_p_pid[ind]->Fill(phi, p);
	proton_phi_pid[ind]->Fill(phi);
      }
    }    
    else if(pdg == 211){
      if (type == "trk"){
	piplus_p_costheta_trk[ind]->Fill(p, costheta);
	piplus_phi_costheta_trk[ind]->Fill(phi, costheta);
	piplus_p_trk[ind]->Fill(p);
	piplus_costheta_trk[ind]->Fill(costheta);
	piplus_phi_p_trk[ind]->Fill(phi, p);
	piplus_phi_trk[ind]->Fill(phi);
      }
      else if (type == "PID") {
	piplus_p_costheta_pid[ind]->Fill(p, costheta);
	piplus_phi_costheta_pid[ind]->Fill(phi, costheta);
	piplus_p_pid[ind]->Fill(p);
	piplus_costheta_pid[ind]->Fill(costheta);
	piplus_phi_p_pid[ind]->Fill(phi, p);
	piplus_phi_pid[ind]->Fill(phi);
      }
    }
    else if(pdg == -211){
      if (type == "trk"){
	piminus_p_costheta_trk[ind]->Fill(p, costheta);
	piminus_phi_costheta_trk[ind]->Fill(phi, costheta);
	piminus_p_trk[ind]->Fill(p);
	piminus_costheta_trk[ind]->Fill(costheta);
	piminus_phi_p_trk[ind]->Fill(phi, p);
	piminus_phi_trk[ind]->Fill(phi);
      }
      else if (type == "PID") {
	piminus_p_costheta_pid[ind]->Fill(p, costheta);
	piminus_phi_costheta_pid[ind]->Fill(phi, costheta);
	piminus_p_pid[ind]->Fill(p);
	piminus_costheta_pid[ind]->Fill(costheta);
	piminus_phi_p_pid[ind]->Fill(phi, p);
	piminus_phi_pid[ind]->Fill(phi);
     }
   }
   else if(pdg == 111){
     if (type == "trk"){
       pi0_p_costheta_trk[ind]->Fill(p, costheta);
       pi0_phi_costheta_trk[ind]->Fill(phi, costheta);
       pi0_p_trk[ind]->Fill(p);
       pi0_costheta_trk[ind]->Fill(costheta);
       pi0_phi_p_trk[ind]->Fill(phi, p);
       pi0_phi_trk[ind]->Fill(phi);
     }
     else if (type == "PID") {
       pi0_p_costheta_pid[ind]->Fill(p, costheta);
       pi0_phi_costheta_pid[ind]->Fill(phi, costheta);
       pi0_p_pid[ind]->Fill(p);
       pi0_costheta_pid[ind]->Fill(costheta);
       pi0_phi_p_pid[ind]->Fill(phi, p);
       pi0_phi_pid[ind]->Fill(phi);
    }
  }

}
