#define make_effs_cxx
#include "make_effs.h"

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

void dune_dst::Loop(int n_evt, char* tag,char * fOutFileName)
{

  if (fChain == 0) return;

  Long64_t nentries = n_evt;
  if(nentries<=0) nentries = fChain->GetEntries();

  Long64_t nbytes = 0, nb = 0;

  //Need array[3] for acc, all, eff
  //Need to decide if vtx info is 
  // 
  // 
  TFile *fout = new TFile(fOutFileName,"RECREATE"); 

  const int nPBins = 25;
  const double PLow = 0;
  const double PUp = 0.1;

  const int nCosthetaBins = 100;
  const double CosthetaLow = -1;
  const double CosthetaUp = 1;

  const int nPhiBins = 100;
  const double PhiLow = 0;
  const double PhiUp = 360;

  for ( int i = 0; i < 3; i++ ) {
    ////p costheta
    muon_p_costheta_trk[i] = new TH2D(Form("muon_p_costheta_trk_%d", i), 
				      "#mu Track Efficiency;p_{#mu} (GeV/c);cos #theta_{#mu}", 
				      nPBins, PLow, PUp,
				      nCosthetaBins, CosthetaLow, CosthetaUp); 
    piplus_p_costheta_trk[i] = new TH2D(Form("piplus_p_costheta_trk_%d", i),
					"#pi^{+} Track Efficiency;p_{#pi} (GeV/c);cos #theta_{#pi}", 
					nPBins, PLow, PUp,
					nCosthetaBins, CosthetaLow, CosthetaUp);
    piminus_p_costheta_trk[i] = new TH2D(Form("piminus_p_costheta_trk_%d", i),
					 "#pi^{-} Track Efficiency;p_{#pi} (GeV/c);cos #theta_{#pi}",
					 nPBins, PLow, PUp,
					 nCosthetaBins, CosthetaLow, CosthetaUp);
    pi0_p_costheta_trk[i] = new TH2D(Form("pi0_p_costheta_trk_%d", i), 
				     "#pi^{0} Track Efficiency;p_{#pi} (GeV/c);cos #theta_{#pi}", 
				     nPBins, PLow, PUp,
				     nCosthetaBins, CosthetaLow, CosthetaUp);
    proton_p_costheta_trk[i] = new TH2D(Form("proton_p_costheta_trk_%d", i), 
					"p Track Efficiency;p_{p} (GeV/c);cos #theta_{p}", 
					nPBins, PLow, PUp,
					nCosthetaBins, CosthetaLow, CosthetaUp);

    muon_p_costheta_pid[i] = new TH2D(Form("muon_p_costheta_pid_%d", i), 
				      "#mu Track Efficiency;p_{#mu} (GeV/c);cos #theta_{#mu}", 
				      nPBins, PLow, PUp,
				      nCosthetaBins, CosthetaLow, CosthetaUp); 

    piplus_p_costheta_pid[i] = new TH2D(Form("piplus_p_costheta_pid_%d", i),
					"#pi^{+} PID Efficiency;p_{#pi} (GeV/c);cos #theta_{#pi}", 
					nPBins, PLow, PUp,
					nCosthetaBins, CosthetaLow, CosthetaUp);
    piminus_p_costheta_pid[i] = new TH2D(Form("piminus_p_costheta_pid_%d", i),
					 "#pi^{-} PID Efficiency;p_{#pi} (GeV/c);cos #theta_{#pi}",
					 nPBins, PLow, PUp,
					 nCosthetaBins, CosthetaLow, CosthetaUp);
    pi0_p_costheta_pid[i] = new TH2D(Form("pi0_p_costheta_pid_%d", i), 
				     "#pi^{0} PID Efficiency;p_{#pi} (GeV/c);cos #theta_{#pi}", 
				     nPBins, PLow, PUp,
				     nCosthetaBins, CosthetaLow, CosthetaUp);
    proton_p_costheta_pid[i] = new TH2D(Form("proton_p_costheta_pid_%d", i), 
					"p PID Efficiency;p_{p} (GeV/c);cos #theta_{p}", 
					nPBins, PLow, PUp,
					nCosthetaBins, CosthetaLow, CosthetaUp);


    //phi costheta
    muon_phi_costheta_trk[i] = new TH2D(Form("muon_phi_costheta_trk_%d", i), 
					"#mu Track Efficiency;#phi_{#mu};cos #theta_{#mu}", 
					nPhiBins, PhiLow, PhiUp,
					nCosthetaBins, CosthetaLow, CosthetaUp); 
    piplus_phi_costheta_trk[i] = new TH2D(Form("piplus_phi_costheta_trk_%d", i),
					  "#pi^{+} Track Efficiency;#phi_{#pi};cos #theta_{#pi}", 
					  nPhiBins, PhiLow, PhiUp,
					  nCosthetaBins, CosthetaLow, CosthetaUp);
    piminus_phi_costheta_trk[i] = new TH2D(Form("piminus_phi_costheta_trk_%d", i),
					   "#pi^{-} Track Efficiency;#phi_{#pi};cos #theta_{#pi}",
					   nPhiBins, PhiLow, PhiUp,
					   nCosthetaBins, CosthetaLow, CosthetaUp);
    pi0_phi_costheta_trk[i] = new TH2D(Form("pi0_phi_costheta_trk_%d", i), 
				       "#pi^{0} Track Efficiency;#phi_{#pi};cos #theta_{#pi}", 
				       nPhiBins, PhiLow, PhiUp,
				       nCosthetaBins, CosthetaLow, CosthetaUp);
    proton_phi_costheta_trk[i] = new TH2D(Form("proton_phi_costheta_trk_%d", i), 
					  "p Track Efficiency;#phi_{p};cos #theta_{p}", 
					  nPhiBins, PhiLow, PhiUp,
					  nCosthetaBins, CosthetaLow, CosthetaUp);
    muon_phi_costheta_pid[i] = new TH2D(Form("muon_phi_costheta_pid_%d", i), 
					"#mu Track Efficiency;#phi_{#mu};cos #theta_{#mu}", 
					nPhiBins, PhiLow, PhiUp,
					nCosthetaBins, CosthetaLow, CosthetaUp); 

    piplus_phi_costheta_pid[i] = new TH2D(Form("piplus_phi_costheta_pid_%d", i),
					  "#pi^{+} PID Efficiency;#phi_{#pi};cos #theta_{#pi}", 
					  nPhiBins, PhiLow, PhiUp,
					  nCosthetaBins, CosthetaLow, CosthetaUp);
    piminus_phi_costheta_pid[i] = new TH2D(Form("piminus_phi_costheta_pid_%d", i),
					   "#pi^{-} PID Efficiency;#phi_{#pi};cos #theta_{#pi}",
					   nPhiBins, PhiLow, PhiUp,
					   nCosthetaBins, CosthetaLow, CosthetaUp);
    pi0_phi_costheta_pid[i] = new TH2D(Form("pi0_phi_costheta_pid_%d", i), 
				       "#pi^{0} PID Efficiency;#phi_{#pi};cos #theta_{#pi}", 
				       nPhiBins, PhiLow, PhiUp,
				       nCosthetaBins, CosthetaLow, CosthetaUp);
    proton_phi_costheta_pid[i] = new TH2D(Form("proton_phi_costheta_pid_%d", i), 
					  "p PID Efficiency;#phi_{p};cos #theta_{p}", 
					  nPhiBins, PhiLow, PhiUp,
					  nCosthetaBins, CosthetaLow, CosthetaUp);

    //phi p
    muon_phi_p_trk[i] = new TH2D(Form("muon_phi_p_trk_%d", i), 
				 "#mu Track Efficiency;#phi_{#mu};p_{#mu} (GeV/c)", 
				 nPhiBins, PhiLow, PhiUp,
				 nPBins, PLow, PUp); 
    piplus_phi_p_trk[i] = new TH2D(Form("piplus_phi_p_trk_%d", i),
					  "#pi^{+} Track Efficiency;#phi_{#pi};p_{#pi} (GeV/c)", 
					  nPhiBins, PhiLow, PhiUp,
					  nPBins, PLow, PUp);
    piminus_phi_p_trk[i] = new TH2D(Form("piminus_phi_p_trk_%d", i),
					   "#pi^{-} Track Efficiency;#phi_{#pi};p_{#pi} (GeV/c)",
					   nPhiBins, PhiLow, PhiUp,
					   nPBins, PLow, PUp);
    pi0_phi_p_trk[i] = new TH2D(Form("pi0_phi_p_trk_%d", i), 
				       "#pi^{0} Track Efficiency;#phi_{#pi};p_{#pi} (GeV/c)", 
				       nPhiBins, PhiLow, PhiUp,
				       nPBins, PLow, PUp);
    proton_phi_p_trk[i] = new TH2D(Form("proton_phi_p_trk_%d", i), 
					  "p Track Efficiency;#phi_{p};p_{p} (GeV/c)", 
					  nPhiBins, PhiLow, PhiUp,
					  nPBins, PLow, PUp);
    muon_phi_p_pid[i] = new TH2D(Form("muon_phi_p_pid_%d", i), 
					"#mu Track Efficiency;#phi_{#mu};p_{#mu} (GeV/c)", 
					nPhiBins, PhiLow, PhiUp,
					nPBins, PLow, PUp); 

    piplus_phi_p_pid[i] = new TH2D(Form("piplus_phi_p_pid_%d", i),
					  "#pi^{+} PID Efficiency;#phi_{#pi};p_{#pi} (GeV/c)", 
					  nPhiBins, PhiLow, PhiUp,
					  nPBins, PLow, PUp);
    piminus_phi_p_pid[i] = new TH2D(Form("piminus_phi_p_pid_%d", i),
					   "#pi^{-} PID Efficiency;#phi_{#pi};p_{#pi} (GeV/c)",
					   nPhiBins, PhiLow, PhiUp,
					   nPBins, PLow, PUp);
    pi0_phi_p_pid[i] = new TH2D(Form("pi0_phi_p_pid_%d", i), 
				       "#pi^{0} PID Efficiency;#phi_{#pi};p_{#pi} (GeV/c)", 
				       nPhiBins, PhiLow, PhiUp,
				       nPBins, PLow, PUp);
    proton_phi_p_pid[i] = new TH2D(Form("proton_phi_p_pid_%d", i), 
					  "p PID Efficiency;#phi_{p};p_{p} (GeV/c)", 
					  nPhiBins, PhiLow, PhiUp,
					  nPBins, PLow, PUp);

    //costheta
    muon_costheta_trk[i] = new TH1D(Form("muon_costheta_trk_%d", i), 
				    "#mu Track Efficiency;cos #theta_{#mu}",
				    nCosthetaBins, CosthetaLow, CosthetaUp); 
    piplus_costheta_trk[i] = new TH1D(Form("piplus_costheta_trk_%d", i),
				      "#pi^{+} Track Efficiency;cos #theta_{#pi}",
				      nCosthetaBins, CosthetaLow, CosthetaUp);
    piminus_costheta_trk[i] = new TH1D(Form("piminus_costheta_trk_%d", i),
				       "#pi^{-} Track Efficiency;cos #theta_{#pi}",
				       nCosthetaBins, CosthetaLow, CosthetaUp);
    pi0_costheta_trk[i] = new TH1D(Form("pi0_costheta_trk_%d", i), 
				   "#pi^{0} Track Efficiency;cos #theta_{#pi}",
				   nCosthetaBins, CosthetaLow, CosthetaUp);
    proton_costheta_trk[i] = new TH1D(Form("proton_costheta_trk_%d", i), 
				      "p Track Efficiency;cos #theta_{p}",
				      nCosthetaBins, CosthetaLow, CosthetaUp);

    muon_costheta_pid[i] = new TH1D(Form("muon_costheta_pid_%d", i), 
				    "#mu PID Efficiency;cos #theta_{#mu}",
				    nCosthetaBins, CosthetaLow, CosthetaUp); 
    piplus_costheta_pid[i] = new TH1D(Form("piplus_costheta_pid_%d", i),
				      "#pi^{+} PID Efficiency;cos #theta_{#pi}",
				      nCosthetaBins, CosthetaLow, CosthetaUp);
    piminus_costheta_pid[i] = new TH1D(Form("piminus_costheta_pid_%d", i),
				       "#pi^{-} PID Efficiency;cos #theta_{#pi}",
				       nCosthetaBins, CosthetaLow, CosthetaUp);
    pi0_costheta_pid[i] = new TH1D(Form("pi0_costheta_pid_%d", i), 
				   "#pi^{0} PID Efficiency;cos #theta_{#pi}",
				   nCosthetaBins, CosthetaLow, CosthetaUp);
    proton_costheta_pid[i] = new TH1D(Form("proton_costheta_pid_%d", i), 
				      "p PID Efficiency;cos #theta_{p}",
				      nCosthetaBins, CosthetaLow, CosthetaUp);

    //p
    muon_p_trk[i] = new TH1D(Form("muon_p_trk_%d", i), 
			     "#mu Track Efficiency;p_{#mu} (GeV/c)",
			     nPBins, PLow, PUp); 
    piplus_p_trk[i] = new TH1D(Form("piplus_p_trk_%d", i),
			       "#pi^{+} Track Efficiency;p_{#pi} (GeV/c)",
			       nPBins, PLow, PUp);
    piminus_p_trk[i] = new TH1D(Form("piminus_p_trk_%d", i),
				"#pi^{-} Track Efficiency;p_{#pi} (GeV/c)",
				nPBins, PLow, PUp);
    pi0_p_trk[i] = new TH1D(Form("pi0_p_trk_%d", i), 
			    "#pi^{0} Track Efficiency;p_{#pi} (GeV/c)",
			    nPBins, PLow, PUp);
    proton_p_trk[i] = new TH1D(Form("proton_p_trk_%d", i), 
			       "p Track Efficiency;p_{p} (GeV/c)",
			       nPBins, PLow, PUp);
    muon_p_pid[i] = new TH1D(Form("muon_p_pid_%d", i), 
			     "#mu PID Efficiency;p_{#mu} (GeV/c)",
			     nPBins, PLow, PUp); 
    piplus_p_pid[i] = new TH1D(Form("piplus_p_pid_%d", i),
			       "#pi^{+} PID Efficiency;p_{#pi} (GeV/c)",
			       nPBins, PLow, PUp);
    piminus_p_pid[i] = new TH1D(Form("piminus_p_pid_%d", i),
				"#pi^{-} PID Efficiency;p_{#pi} (GeV/c)",
				nPBins, PLow, PUp);
    pi0_p_pid[i] = new TH1D(Form("pi0_p_pid_%d", i), 
			    "#pi^{0} PID Efficiency;p_{#pi} (GeV/c)",
			    nPBins, PLow, PUp);
    proton_p_pid[i] = new TH1D(Form("proton_p_pid_%d", i), 
			       "p PID Efficiency;p_{p} (GeV/c)",
			       nPBins, PLow, PUp);

    //phi
    muon_phi_trk[i] = new TH1D(Form("muon_phi_trk_%d", i), 
			       "#mu Track Efficiency;#phi_{#mu} (degrees)",
			       nPhiBins, PhiLow, PhiUp); 
    piplus_phi_trk[i] = new TH1D(Form("piplus_phi_trk_%d", i),
				 "#pi^{+} Track Efficiency;#phi_{#pi} (degrees)",
				 nPhiBins, PhiLow, PhiUp);
    piminus_phi_trk[i] = new TH1D(Form("piminus_phi_trk_%d", i),
				  "#pi^{-} Track Efficiency;#phi_{#pi} (degrees)",
				  nPhiBins, PhiLow, PhiUp);
    pi0_phi_trk[i] = new TH1D(Form("pi0_phi_trk_%d", i), 
			      "#pi^{0} Track Efficiency;#phi_{#pi} (degrees)",
			      nPhiBins, PhiLow, PhiUp);
    proton_phi_trk[i] = new TH1D(Form("proton_phi_trk_%d", i), 
				 "p Track Efficiency;#phi_{p} (degrees)",
				 nPhiBins, PhiLow, PhiUp);
    muon_phi_pid[i] = new TH1D(Form("muon_phi_pid_%d", i), 
			       "#mu PID Efficiency;#phi_{#mu} (degrees)",
			       nPhiBins, PhiLow, PhiUp); 
    piplus_phi_pid[i] = new TH1D(Form("piplus_phi_pid_%d", i),
				 "#pi^{+} PID Efficiency;#phi_{#pi} (degrees)",
				 nPhiBins, PhiLow, PhiUp);
    piminus_phi_pid[i] = new TH1D(Form("piminus_phi_pid_%d", i),
				  "#pi^{-} PID Efficiency;#phi_{#pi} (degrees)",
				  nPhiBins, PhiLow, PhiUp);
    pi0_phi_pid[i] = new TH1D(Form("pi0_phi_pid_%d", i), 
			      "#pi^{0} PID Efficiency;#phi_{#pi} (degrees)",
			      nPhiBins, PhiLow, PhiUp);
    proton_phi_pid[i] = new TH1D(Form("proton_phi_pid_%d", i), 
				 "p PID Efficiency;#phi_{p} (degrees)",
				 nPhiBins, PhiLow, PhiUp);

  }    


  nentries = 1000; //Testing
  for( Long64_t jentry=0; jentry < nentries; jentry++ ) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
       
    if( jentry % 100000 == 0 ) printf( "Entry %lld of %lld...\n", jentry, nentries );

    Int_t targetZ = target/10000%1000;
    // Only use events on Ar (Z = 18)
    /* if ( targetZ != 18 ) { */
    /*   continue; */
    /* } */

    //Idea: there is info about correctly IDing the particles
    //  so should we include wrong Ereco (due to using the
    //  wrong mass) in this?
    //
    //  
    //See page 16 of the NDTF report!!
    
    total_int->Fill(Ev);

    // Loop over FS particles
    bool is_reconstructed = false;
    bool is_reconstructed_correctly = false;


    float proton_HM = 0.;
    float proton_SHM = 0.;
    float piplus_HM = 0.;
    float piplus_SHM = 0.;
    float piminus_HM = 0.;
    float piminus_SHM = 0.;

    int proton_HM_index = -999;
    int proton_SHM_index = -999;
    int piplus_HM_index = -999;
    int piplus_SHM_index = -999;
    int piminus_HM_index = -999;
    int piminus_SHM_index = -999;


    for( int i = 0; i < trueFSParticles; ++i ) {      
 
      int trackID = trueFSParticles_id[i];
      int pdg = trueFSParticles_pdg[i];
      float px_true = trueFSParticles_momentum[i][0];
      float py_true = trueFSParticles_momentum[i][1];
      float pz_true = trueFSParticles_momentum[i][2];
      float p_true  = sqrt(px_true*px_true + py_true*py_true + pz_true*pz_true);

      float costheta_true = pz_true/p_true; 

      float phi_true = 180.*atan2(py_true, px_true)/pi + 180.;

      fillHists("trk", pdg, 1, p_true, costheta_true, phi_true);



      //Smearing stuff
      float Etrue;
      if(pdg == 2112)Etrue = trueFSParticles_energy[i]-.938;
      else if(pdg == 2212){
        Etrue = trueFSParticles_energy[i] - .938;
        if(p_true > proton_HM){
          proton_HM = p_true;proton_HM_index = i;}
        else if(p_true > proton_SHM){
          proton_SHM = p_true;proton_SHM_index = i;}
      }
      else if(pdg == 211){
        Etrue = trueFSParticles_energy[i];
        if(p_true > piplus_HM){
          piplus_HM = p_true;piplus_HM_index = i;}
        else if(p_true > piplus_SHM){
          piplus_SHM = p_true;piplus_SHM_index = i;}
      }
      else if(pdg == -211){
        Etrue = trueFSParticles_energy[i];
        if(p_true > piminus_HM){
         piminus_HM = p_true;piminus_HM_index = i;}
        else if(p_true > piminus_SHM){
         piminus_SHM = p_true;piminus_SHM_index = i;}      
      }
      else{
        Etrue = trueFSParticles_energy[i];
      }
      ////


      // check for matching reco particle
      for( int j = 0; j < recoFSParticles; ++j ) {
        if( recoFSParticles_id[j] == trackID ) {
          is_reconstructed = true; // something was reconstructed

	  fillHists("trk", pdg, 0, p_true, costheta_true, phi_true);
	  fillHists("PID", pdg, 1, p_true, costheta_true, phi_true);

          int reco_pdg = recoFSParticles_pdg[j];

          if(reco_pdg == 2212){
            had_smear->Fill(Etrue,recoFSParticles_energy[j]-.938);
          }
          else if(reco_pdg == 211 || reco_pdg == -211 || reco_pdg == 111){
            had_smear->Fill(Etrue,recoFSParticles_energy[j]);
          }

          if( reco_pdg == pdg ) {
            is_reconstructed_correctly = true; // reco with correct particle type
	    fillHists("PID", pdg, 0, p_true, costheta_true, phi_true);
            
            if(abs(pdg) == 13) muon_smear->Fill(Etrue,recoFSParticles_energy[j]);
            
          }

          break;

        }
      }

      if(!is_reconstructed){
        if(abs(pdg) != 13) had_smear->Fill(Etrue,0.); 
      }
    }
    
    float proton_HM_etrue = trueFSParticles_energy[proton_HM_index] -.938;

    bool proton_HM_is_reconstructed = false;

/*    for( int j = 0; j < recoFSParticles; ++j ) {

      if(proton_HM_index != -999){
        if( recoFSParticles_id[j] == proton_HM_index ) {
          proton_HM_is_reconstructed == true; 
          int reco_pdg = recoFSParticles_pdg[j];
          if(reco_pdg == 2212){
            proton_HM_smear->Fill(proton_HM_etrue,recoFSParticles_energy[j] - .938);        
          }
          else
          {
            proton_HM_smear->Fill(proton_HM_etrue,recoFSParticles_energy[j]);
          }
        }
      }

    }

    if(!proton_HM_is_reconstructed){
      proton_HM_smear->Fill(proton_HM_etrue,0.);
    }*/


    if ( is_reconstructed ) {
      reco_int->Fill(Ev);
    }
    if ( is_reconstructed_correctly ) {
      pid_int->Fill(Ev);
    }
    if ( targetZ == 18 ) {
      argon_int->Fill(Ev);
    }

    
  }
  fout->cd();

  doEff(muon_p_costheta_trk);
  doEff(piplus_p_costheta_trk);
  doEff(piminus_p_costheta_trk);
  doEff(pi0_p_costheta_trk);
  doEff(proton_p_costheta_trk);

  doEff(muon_phi_costheta_trk);
  doEff(piplus_phi_costheta_trk);
  doEff(piminus_phi_costheta_trk);
  doEff(pi0_phi_costheta_trk);
  doEff(proton_phi_costheta_trk);

  doEff(muon_phi_p_trk);
  doEff(piplus_phi_p_trk);
  doEff(piminus_phi_p_trk);
  doEff(pi0_phi_p_trk);
  doEff(proton_phi_p_trk);

  doEff(muon_p_costheta_pid);
  doEff(piplus_p_costheta_pid);
  doEff(piminus_p_costheta_pid);
  doEff(pi0_p_costheta_pid);
  doEff(proton_p_costheta_pid);

  doEff(muon_phi_costheta_pid);
  doEff(piplus_phi_costheta_pid);
  doEff(piminus_phi_costheta_pid);
  doEff(pi0_phi_costheta_pid);
  doEff(proton_phi_costheta_pid);

  doEff(muon_phi_p_pid);
  doEff(piplus_phi_p_pid);
  doEff(piminus_phi_p_pid);
  doEff(pi0_phi_p_pid);
  doEff(proton_phi_p_pid);  

  TEfficiency* muon_p_trk_eff = new TEfficiency(*muon_p_trk[0], 
					    *muon_p_trk[1]);
  TEfficiency* piplus_p_trk_eff = new TEfficiency(*piplus_p_trk[0], 
					      *piplus_p_trk[1]);
  TEfficiency* piminus_p_trk_eff = new TEfficiency(*piminus_p_trk[0], 
					       *piminus_p_trk[1]);
  TEfficiency* pi0_p_trk_eff = new TEfficiency(*pi0_p_trk[0], 
					   *pi0_p_trk[1]);
  TEfficiency* proton_p_trk_eff = new TEfficiency(*proton_p_trk[0], 
					      *proton_p_trk[1]);

  TEfficiency* muon_costheta_trk_eff = new TEfficiency(*muon_costheta_trk[0], 
						   *muon_costheta_trk[1]);
  TEfficiency* piplus_costheta_trk_eff = new TEfficiency(*piplus_costheta_trk[0], 
						     *piplus_costheta_trk[1]);
  TEfficiency* piminus_costheta_trk_eff = new TEfficiency(*piminus_costheta_trk[0], 
						      *piminus_costheta_trk[1]);
  TEfficiency* pi0_costheta_trk_eff = new TEfficiency(*pi0_costheta_trk[0], 
						  *pi0_costheta_trk[1]);
  TEfficiency* proton_costheta_trk_eff = new TEfficiency(*proton_costheta_trk[0], 
						     *proton_costheta_trk[1]);

  TEfficiency* muon_phi_trk_eff = new TEfficiency(*muon_phi_trk[0], 
						   *muon_phi_trk[1]);
  TEfficiency* piplus_phi_trk_eff = new TEfficiency(*piplus_phi_trk[0], 
						     *piplus_phi_trk[1]);
  TEfficiency* piminus_phi_trk_eff = new TEfficiency(*piminus_phi_trk[0], 
						      *piminus_phi_trk[1]);
  TEfficiency* pi0_phi_trk_eff = new TEfficiency(*pi0_phi_trk[0], 
						  *pi0_phi_trk[1]);
  TEfficiency* proton_phi_trk_eff = new TEfficiency(*proton_phi_trk[0], 
						     *proton_phi_trk[1]);

  TEfficiency* muon_p_pid_eff = new TEfficiency(*muon_p_pid[0], 
					    *muon_p_pid[1]);
  TEfficiency* piplus_p_pid_eff = new TEfficiency(*piplus_p_pid[0], 
					      *piplus_p_pid[1]);
  TEfficiency* piminus_p_pid_eff = new TEfficiency(*piminus_p_pid[0], 
					       *piminus_p_pid[1]);
  TEfficiency* pi0_p_pid_eff = new TEfficiency(*pi0_p_pid[0], 
					   *pi0_p_pid[1]);
  TEfficiency* proton_p_pid_eff = new TEfficiency(*proton_p_pid[0], 
					      *proton_p_pid[1]);
  
  TEfficiency* muon_costheta_pid_eff = new TEfficiency(*muon_costheta_pid[0], 
						   *muon_costheta_pid[1]);
  TEfficiency* piplus_costheta_pid_eff = new TEfficiency(*piplus_costheta_pid[0], 
						     *piplus_costheta_pid[1]);
  TEfficiency* piminus_costheta_pid_eff = new TEfficiency(*piminus_costheta_pid[0], 
						      *piminus_costheta_pid[1]);
  TEfficiency* pi0_costheta_pid_eff = new TEfficiency(*pi0_costheta_pid[0], 
						  *pi0_costheta_pid[1]);
  TEfficiency* proton_costheta_pid_eff = new TEfficiency(*proton_costheta_pid[0], 
						     *proton_costheta_pid[1]);

  TEfficiency* muon_phi_pid_eff = new TEfficiency(*muon_phi_pid[0], 
						   *muon_phi_pid[1]);
  TEfficiency* piplus_phi_pid_eff = new TEfficiency(*piplus_phi_pid[0], 
						     *piplus_phi_pid[1]);
  TEfficiency* piminus_phi_pid_eff = new TEfficiency(*piminus_phi_pid[0], 
						      *piminus_phi_pid[1]);
  TEfficiency* pi0_phi_pid_eff = new TEfficiency(*pi0_phi_pid[0], 
						  *pi0_phi_pid[1]);
  TEfficiency* proton_phi_pid_eff = new TEfficiency(*proton_phi_pid[0], 
						     *proton_phi_pid[1]);
  
  TEfficiency* argon_purity = new TEfficiency(*argon_int, *total_int);
  argon_purity->SetName("argon_purity");
  TEfficiency* reco_purity = new TEfficiency(*reco_int, *total_int);
  reco_purity->SetName("reco_purity");
  TEfficiency* pid_purity = new TEfficiency(*pid_int, *total_int);
  pid_purity->SetName("pid_purity");

  fout->mkdir("1D_Efficiencies");
  fout->cd("1D_Efficiencies");

  fout->mkdir("1D_Efficiencies/p_trk");
  fout->cd("1D_Efficiencies/p_trk");
  muon_p_trk_eff->Write();
  piplus_p_trk_eff->Write();
  piminus_p_trk_eff->Write();
  pi0_p_trk_eff->Write();
  proton_p_trk_eff->Write();

  fout->mkdir("1D_Efficiencies/costheta_trk");
  fout->cd("1D_Efficiencies/costheta_trk");
  muon_costheta_trk_eff->Write();
  piplus_costheta_trk_eff->Write();
  piminus_costheta_trk_eff->Write();
  pi0_costheta_trk_eff->Write();
  proton_costheta_trk_eff->Write();

  fout->mkdir("1D_Efficiencies/phi_trk");
  fout->cd("1D_Efficiencies/phi_trk");
  muon_phi_trk_eff->Write();
  piplus_phi_trk_eff->Write();
  piminus_phi_trk_eff->Write();
  pi0_phi_trk_eff->Write();
  proton_phi_trk_eff->Write();

  fout->mkdir("1D_Efficiencies/p_pid");
  fout->cd("1D_Efficiencies/p_pid");
  muon_p_pid_eff->Write();
  piplus_p_pid_eff->Write();
  piminus_p_pid_eff->Write();
  pi0_p_pid_eff->Write();
  proton_p_pid_eff->Write();

  fout->mkdir("1D_Efficiencies/costheta_pid");
  fout->cd("1D_Efficiencies/costheta_pid");
  muon_costheta_pid_eff->Write();
  piplus_costheta_pid_eff->Write();
  piminus_costheta_pid_eff->Write();
  pi0_costheta_pid_eff->Write();
  proton_costheta_pid_eff->Write();

  fout->mkdir("1D_Efficiencies/phi_pid");
  fout->cd("1D_Efficiencies/phi_pid");
  muon_phi_pid_eff->Write();
  piplus_phi_pid_eff->Write();
  piminus_phi_pid_eff->Write();
  pi0_phi_pid_eff->Write();
  proton_phi_pid_eff->Write();

  fout->cd();
  fout->mkdir("2D_Efficiencies");
  fout->cd("2D_Efficiencies");

  // 2D track efficiencies
  fout->mkdir("2D_Efficiencies/p_costheta_trk");
  fout->cd("2D_Efficiencies/p_costheta_trk");
  muon_p_costheta_trk[2]->Write();
  piplus_p_costheta_trk[2]->Write();
  piminus_p_costheta_trk[2]->Write();
  pi0_p_costheta_trk[2]->Write();
  proton_p_costheta_trk[2]->Write();
    
  fout->mkdir("2D_Efficiencies/phi_costheta_trk");
  fout->cd("2D_Efficiencies/phi_costheta_trk");
  muon_phi_costheta_trk[2]->Write();
  piplus_phi_costheta_trk[2]->Write();
  piminus_phi_costheta_trk[2]->Write();
  pi0_phi_costheta_trk[2]->Write();
  proton_phi_costheta_trk[2]->Write();

  fout->mkdir("2D_Efficiencies/phi_p_trk");
  fout->cd("2D_Efficiencies/phi_p_trk");
  muon_phi_p_trk[2]->Write();
  piplus_phi_p_trk[2]->Write();
  piminus_phi_p_trk[2]->Write();
  pi0_phi_p_trk[2]->Write();
  proton_phi_p_trk[2]->Write();

  // 2D PID efficiecies
  fout->mkdir("2D_Efficiencies/p_costheta_pid");
  fout->cd("2D_Efficiencies/p_costheta_pid");
  muon_p_costheta_pid[2]->Write(); 
  piplus_p_costheta_pid[2]->Write();
  piminus_p_costheta_pid[2]->Write();
  pi0_p_costheta_pid[2]->Write();
  proton_p_costheta_pid[2]->Write();

  fout->mkdir("2D_Efficiencies/phi_costheta_pid");
  fout->cd("2D_Efficiencies/phi_costheta_pid");
  muon_phi_costheta_pid[2]->Write(); 
  piplus_phi_costheta_pid[2]->Write();
  piminus_phi_costheta_pid[2]->Write();
  pi0_phi_costheta_pid[2]->Write();
  proton_phi_costheta_pid[2]->Write();

  fout->mkdir("2D_Efficiencies/phi_p_pid");
  fout->cd("2D_Efficiencies/phi_p_pid");
  muon_phi_p_pid[2]->Write(); 
  piplus_phi_p_pid[2]->Write();
  piminus_phi_p_pid[2]->Write();
  pi0_phi_p_pid[2]->Write();
  proton_phi_p_pid[2]->Write();

  fout->cd();
  fout->mkdir("Event_Rates");
  fout->cd("Event_Rates");
  total_int->Write();
  argon_int->Write();
  reco_int->Write();
  pid_int->Write();
  
  fout->cd();
  fout->mkdir("Smears");
  fout->cd("Smears");
  had_smear->Write();
  muon_smear->Write();

  argon_purity->Write();
  reco_purity->Write();
  pid_purity->Write();

  fout->Close();

}
