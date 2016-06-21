////////////////////////////////////////////////////////////////////////
// $Id: NuWroGen_module.cc,v 1.4 2010/04/27 19:48:46 echurch Exp $
//
//
// NuWro neutrino event generator (parser)
//
// echurch@fnal.gov
//
////////////////////////////////////////////////////////////////////////

#ifndef EVGEN_NUWROGEN_H
#define EVGEN_NUWROGEN_H

#include <sys/stat.h> 
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TBranch.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TStopwatch.h"


// LArSoft includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/SummaryData/RunData.h"

// #include "NWtree.h"

class TH1F;
class TH2F;

// NuWro include - event1dict.h includes event1.h, 
// the definition for the NuWro event class
//#include "event1dict.h"

namespace simb { class MCTruth; }

namespace evgen {


  /// A module to check the results from the Monte Carlo generator
  class NuWroGen : public art::EDProducer {
  public:
    explicit NuWroGen(fhicl::ParameterSet const &pset);
    virtual ~NuWroGen();                        

    void produce(art::Event& evt);  
    void beginJob();
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& p);
    void endJob();

  private:

        std::string ParticleStatus(int StatusCode);
        std::string ReactionChannel(int ccnc,int mode);
    
        void FillHistograms(const simb::MCTruth &mc);

	std::string         fFileName;
	std::string         fTreeName;
	std::ifstream      *fEventFile;

	TStopwatch          fStopwatch;      ///keep track of how long it takes to run the job
	
	std::string fNuWroModuleLabel;
	int fEventNumberOffset;

	std::vector<double> fxsecFluxWtd;
	double fxsecTotal;

	TH1F* fGenerated[6];  ///< Spectra as generated
	
	TH1F* fVertexX;    ///< vertex location of generated events in x
	TH1F* fVertexY;    ///< vertex location of generated events in y
	TH1F* fVertexZ;    ///< vertex location of generated events in z
	
	TH2F* fVertexXY;   ///< vertex location in xy
	TH2F* fVertexXZ;   ///< vertex location in xz
	TH2F* fVertexYZ;   ///< vertex location in yz
	
	TH1F* fDCosX;      ///< direction cosine in x
	TH1F* fDCosY;      ///< direction cosine in y
	TH1F* fDCosZ;      ///< direction cosine in z
	
	TH1F* fMuMomentum; ///< momentum of outgoing muons
	TH1F* fMuDCosX;    ///< direction cosine of outgoing mu in x
	TH1F* fMuDCosY;    ///< direction cosine of outgoing mu in y
	TH1F* fMuDCosZ;    ///< direction cosine of outgoing mu in z
	
	TH1F* fEMomentum;  ///< momentum of outgoing electrons
	TH1F* fEDCosX;     ///< direction cosine of outgoing e in x
	TH1F* fEDCosY;     ///< direction cosine of outgoing e in y
	TH1F* fEDCosZ;     ///< direction cosine of outgoing e in z
	
	TH1F* fCCMode;      ///< CC interaction mode
	TH1F* fNCMode;      ///< CC interaction mode
	TH1F* fDyn;         ///< mode in detail
	TH1F* fWeight;         ///< NuWro Wt
	TH1F* fWeightNW;         ///< NuWro Wt
	TH1F* fDynNew;         
	TH1F* fDynNewThresh;   
	TH2F* f2DynNew;        
	TH2F* f2DynNewThresh;  
	
	TH1F* fDeltaE;     ///< difference in neutrino energy from MCTruth::Enu() vs TParticle
	TH1F* fECons;      ///< histogram to determine if energy is conserved in the event



	  unsigned int countFile;
	  event     *NuWroTTree;
	  TChain          *ch;

	  // here come all the NuWro variables...
	  /*
	  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	  Int_t           fCurrent; //!current Tree number in a TChain
	
	  UInt_t          fUniqueID;
	
	  UInt_t          fBits;
	
	  Bool_t          flag_coh;
     
	  Bool_t          flag_qel;
     
	  Bool_t          flag_dis;
     
	  Bool_t          flag_res;
     
	  Bool_t          flag_nc;
     
	  Bool_t          flag_cc;
     
	  Bool_t          flag_anty;
     
	  Int_t           par_random_seed;
     
	  Int_t           par_number_of_events;
     
	  Int_t           par_number_of_test_events;
     
	  Int_t           par_user_events;
     
	  Int_t           par_beam_type;
     
	  std::string          par_beam_energy_string;
     
	  Int_t           par_beam_particle;
     
	  Double_t        par_beam_direction_x;
     
	  Double_t        par_beam_direction_y;
     
	  Double_t        par_beam_direction_z;
     
	  std::string          par_beam_content_string;
     
	  std::string          par_beam_folder;
     
	  Int_t           par_beam_file_first;
     
	  Int_t           par_beam_file_limit;
     
	  Bool_t          par_beam_weighted;
     
	  std::string          par_beam_file;
     
	  Double_t        par_beam_offset_x;
     
	  Double_t        par_beam_offset_y;
     
	  Double_t        par_beam_offset_z;
     
	  Int_t           par_beam_placement;
     
	  Int_t           par_beam_test_only;
     
	  Int_t           par_nucleus_p;
     
	  Int_t           par_nucleus_n;
     
	  std::string          par_nucleus_density;
     
	  Double_t        par_nucleus_E_b;
     
	  Double_t        par_nucleus_kf;
     
	  Int_t           par_nucleus_target;
     
	  Int_t           par_nucleus_model;
     
	  Int_t           par_target_type;
     
	  std::string          par_target_content_string;
     
	  std::string          par_geo_file;
     
	  std::string          par_geo_name;
     
	  std::string          par_geo_volume;
     
	  Double_t        par_geo_o_x;
     
	  Double_t        par_geo_o_y;
     
	  Double_t        par_geo_o_z;
     
	  Double_t        par_geo_d_x;
     
	  Double_t        par_geo_d_y;
     
	  Double_t        par_geo_d_z;
     
	  Bool_t          par_dyn_qel_cc;
     
	  Bool_t          par_dyn_qel_nc;
     
	  Bool_t          par_dyn_res_cc;
     
	  Bool_t          par_dyn_res_nc;
     
	  Bool_t          par_dyn_dis_cc;
     
	  Bool_t          par_dyn_dis_nc;
     
	  Bool_t          par_dyn_coh_cc;
     
	  Bool_t          par_dyn_coh_nc;
     
	  Int_t           par_qel_vector_ff_set;
     
	  Int_t           par_qel_axial_ff_set;
     
	  Double_t        par_qel_cc_vector_mass;
     
	  Double_t        par_qel_cc_axial_mass;
     
	  Double_t        par_qel_nc_axial_mass;
     
	  Double_t        par_qel_s_axial_mass;
     
	  Int_t           par_qel_strange;
     
	  Int_t           par_qel_strangeEM;
     
	  Double_t        par_delta_s;
     
	  Bool_t          par_flux_correction;
     
	  Int_t           par_delta_FF_set;
     
	  Double_t        par_pion_axial_mass;
     
	  Double_t        par_pion_C5A;
     
	  Double_t        par_res_dis_cut;
     
	  Int_t           par_spp_precision;
     
	  Bool_t          par_coh_mass_correction;
     
	  Bool_t          par_coh_new;
	
	  Bool_t          par_kaskada_on;
	
	  Bool_t          par_kaskada_newangle;
	
	  Bool_t          par_kaskada_redo;
	
	  Bool_t          par_pauli_blocking;
	
	  std::string          par_formation_zone;
	
	  Bool_t          par_first_step;
	
	  Double_t        par_step;
	
	  Int_t           par_xsec;
	
	  Int_t           par_sf_method;
	
	  Bool_t          par_cc_smoothing;
	
	  Bool_t          par_mixed_order;
	
	  Double_t        par_rmin;
	
	  Double_t        par_rmax;
	
	  std::string          par_path_to_data;
	
	  Int_t           in_;
	
	  Double_t        in_t[kMaxin];   //[in_]
	
	  Double_t        in_x[kMaxin];   //[in_]
	
	  Double_t        in_y[kMaxin];   //[in_]
	
	  Double_t        in_z[kMaxin];   //[in_]
	
	  Double_t        in__mass[kMaxin];   //[in_]
	
	  Double_t        in_r_t[kMaxin];   //[in_]
	
	  Double_t        in_r_x[kMaxin];   //[in_]
	
	  Double_t        in_r_y[kMaxin];   //[in_]
	
	  Double_t        in_r_z[kMaxin];   //[in_]
	
	  Int_t           in_pdg[kMaxin];   //[in_]
	
	  Char_t          in_ks[kMaxin];   //[in_]
	
	  Char_t          in_orgin[kMaxin];   //[in_]
	
	  Double_t        in_travelled[kMaxin];   //[in_]
	
	  Int_t           in_id[kMaxin];   //[in_]
	
	  Int_t           in_mother[kMaxin];   //[in_]
	
	  Int_t           in_endproc[kMaxin];   //[in_]
	
	  Double_t        in_fz[kMaxin];   //[in_]
	
	  Double_t        in_pt[kMaxin];   //[in_]
	
	  Int_t           in_mother_pdg[kMaxin];   //[in_]
	
	  Int_t           in_mother_proc[kMaxin];   //[in_]
	
	  Double_t        in_mother_momentum_x[kMaxin];   //[in_]
	
	  Double_t        in_mother_momentum_y[kMaxin];   //[in_]
	
	  Double_t        in_mother_momentum_z[kMaxin];   //[in_]
	
	  Double_t        in_mother_ek[kMaxin];   //[in_]
	
	  Int_t           in_his_nqel[kMaxin];   //[in_]
	
	  Int_t           in_his_nspp[kMaxin];   //[in_]
	
	  Int_t           in_his_ndpp[kMaxin];   //[in_]
	
	  Int_t           in_his_pqel[kMaxin];   //[in_]
	
	  Int_t           in_his_pcex[kMaxin];   //[in_]
	
	  Int_t           in_his_pspp[kMaxin];   //[in_]
	
	  Int_t           in_his_pdpp[kMaxin];   //[in_]
	
	  Int_t           in_his_ptpp[kMaxin];   //[in_]
	
	  Int_t           temp_;
	
	  Double_t        temp_t[kMaxtemp];   //[temp_]
	
	  Double_t        temp_x[kMaxtemp];   //[temp_]
	
	  Double_t        temp_y[kMaxtemp];   //[temp_]
	
	  Double_t        temp_z[kMaxtemp];   //[temp_]
	
	  Double_t        temp__mass[kMaxtemp];   //[temp_]
	
	  Double_t        temp_r_t[kMaxtemp];   //[temp_]
	
	  Double_t        temp_r_x[kMaxtemp];   //[temp_]
	
	  Double_t        temp_r_y[kMaxtemp];   //[temp_]
	
	  Double_t        temp_r_z[kMaxtemp];   //[temp_]
	
	  Int_t           temp_pdg[kMaxtemp];   //[temp_]
	
	  Char_t          temp_ks[kMaxtemp];   //[temp_]
	
	  Char_t          temp_orgin[kMaxtemp];   //[temp_]
	
	  Double_t        temp_travelled[kMaxtemp];   //[temp_]
	
	  Int_t           temp_id[kMaxtemp];   //[temp_]
	
	  Int_t           temp_mother[kMaxtemp];   //[temp_]
	
	  Int_t           temp_endproc[kMaxtemp];   //[temp_]
	
	  Double_t        temp_fz[kMaxtemp];   //[temp_]
	
	  Double_t        temp_pt[kMaxtemp];   //[temp_]
	
	  Int_t           temp_mother_pdg[kMaxtemp];   //[temp_]
	
	  Int_t           temp_mother_proc[kMaxtemp];   //[temp_]
	
	  Double_t        temp_mother_momentum_x[kMaxtemp];   //[temp_]
	
	  Double_t        temp_mother_momentum_y[kMaxtemp];   //[temp_]
	
	  Double_t        temp_mother_momentum_z[kMaxtemp];   //[temp_]
	
	  Double_t        temp_mother_ek[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_nqel[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_nspp[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_ndpp[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_pqel[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_pcex[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_pspp[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_pdpp[kMaxtemp];   //[temp_]
	
	  Int_t           temp_his_ptpp[kMaxtemp];   //[temp_]
	
	  Int_t           out_;
	
	  Double_t        out_t[kMaxout];   //[out_]
	
	  Double_t        out_x[kMaxout];   //[out_]
	
	  Double_t        out_y[kMaxout];   //[out_]
	
	  Double_t        out_z[kMaxout];   //[out_]
	
	  Double_t        out__mass[kMaxout];   //[out_]
	
	  Double_t        out_r_t[kMaxout];   //[out_]
	
	  Double_t        out_r_x[kMaxout];   //[out_]
	
	  Double_t        out_r_y[kMaxout];   //[out_]
	
	  Double_t        out_r_z[kMaxout];   //[out_]
	
	  Int_t           out_pdg[kMaxout];   //[out_]
	
	  Char_t          out_ks[kMaxout];   //[out_]
	
	  Char_t          out_orgin[kMaxout];   //[out_]
	
	  Double_t        out_travelled[kMaxout];   //[out_]
	
	  Int_t           out_id[kMaxout];   //[out_]
	
	  Int_t           out_mother[kMaxout];   //[out_]
	
	  Int_t           out_endproc[kMaxout];   //[out_]
	
	  Double_t        out_fz[kMaxout];   //[out_]
	
	  Double_t        out_pt[kMaxout];   //[out_]
	
	  Int_t           out_mother_pdg[kMaxout];   //[out_]
	
	  Int_t           out_mother_proc[kMaxout];   //[out_]
	
	  Double_t        out_mother_momentum_x[kMaxout];   //[out_]
	
	  Double_t        out_mother_momentum_y[kMaxout];   //[out_]
	
	  Double_t        out_mother_momentum_z[kMaxout];   //[out_]
	
	  Double_t        out_mother_ek[kMaxout];   //[out_]
	
	  Int_t           out_his_nqel[kMaxout];   //[out_]
	
	  Int_t           out_his_nspp[kMaxout];   //[out_]
	
	  Int_t           out_his_ndpp[kMaxout];   //[out_]
	
	  Int_t           out_his_pqel[kMaxout];   //[out_]
	
	  Int_t           out_his_pcex[kMaxout];   //[out_]
	
	  Int_t           out_his_pspp[kMaxout];   //[out_]
	
	  Int_t           out_his_pdpp[kMaxout];   //[out_]
	
	  Int_t           out_his_ptpp[kMaxout];   //[out_]
	
	  Int_t           post_;
	
	  Double_t        post_t[kMaxpost];   //[post_]
	
	  Double_t        post_x[kMaxpost];   //[post_]
	
	  Double_t        post_y[kMaxpost];   //[post_]
	
	  Double_t        post_z[kMaxpost];   //[post_]
	
	  Double_t        post__mass[kMaxpost];   //[post_]
	
	  Double_t        post_r_t[kMaxpost];   //[post_]
	
	  Double_t        post_r_x[kMaxpost];   //[post_]
	
	  Double_t        post_r_y[kMaxpost];   //[post_]
	
	  Double_t        post_r_z[kMaxpost];   //[post_]
	
	  Int_t           post_pdg[kMaxpost];   //[post_]
	
	  Char_t          post_ks[kMaxpost];   //[post_]
	
	  Char_t          post_orgin[kMaxpost];   //[post_]
	
	  Double_t        post_travelled[kMaxpost];   //[post_]
	
	  Int_t           post_id[kMaxpost];   //[post_]
	
	  Int_t           post_mother[kMaxpost];   //[post_]
	
	  Int_t           post_endproc[kMaxpost];   //[post_]
	
	  Double_t        post_fz[kMaxpost];   //[post_]
	
	  Double_t        post_pt[kMaxpost];   //[post_]
	
	  Int_t           post_mother_pdg[kMaxpost];   //[post_]
	
	  Int_t           post_mother_proc[kMaxpost];   //[post_]
	
	  Double_t        post_mother_momentum_x[kMaxpost];   //[post_]
	
	  Double_t        post_mother_momentum_y[kMaxpost];   //[post_]
	
	  Double_t        post_mother_momentum_z[kMaxpost];   //[post_]
	
	  Double_t        post_mother_ek[kMaxpost];   //[post_]
	
	  Int_t           post_his_nqel[kMaxpost];   //[post_]
	
	  Int_t           post_his_nspp[kMaxpost];   //[post_]
	
	  Int_t           post_his_ndpp[kMaxpost];   //[post_]
	
	  Int_t           post_his_pqel[kMaxpost];   //[post_]
	
	  Int_t           post_his_pcex[kMaxpost];   //[post_]
	
	  Int_t           post_his_pspp[kMaxpost];   //[post_]
	
	  Int_t           post_his_pdpp[kMaxpost];   //[post_]
	
	  Int_t           post_his_ptpp[kMaxpost];   //[post_]
	
	  Int_t           all_;
	
	  Double_t        all_t[kMaxall];   //[all_]
	
	  Double_t        all_x[kMaxall];   //[all_]
	
	  Double_t        all_y[kMaxall];   //[all_]
	
	  Double_t        all_z[kMaxall];   //[all_]
	
	  Double_t        all__mass[kMaxall];   //[all_]
	
	  Double_t        all_r_t[kMaxall];   //[all_]
	
	  Double_t        all_r_x[kMaxall];   //[all_]
	
	  Double_t        all_r_y[kMaxall];   //[all_]
	
	  Double_t        all_r_z[kMaxall];   //[all_]
	
	  Int_t           all_pdg[kMaxall];   //[all_]
	
	  Char_t          all_ks[kMaxall];   //[all_]
	
	  Char_t          all_orgin[kMaxall];   //[all_]
	
	  Double_t        all_travelled[kMaxall];   //[all_]
	
	  Int_t           all_id[kMaxall];   //[all_]
	
	  Int_t           all_mother[kMaxall];   //[all_]
	
	  Int_t           all_endproc[kMaxall];   //[all_]
	
	  Double_t        all_fz[kMaxall];   //[all_]
	
	  Double_t        all_pt[kMaxall];   //[all_]
	
	  Int_t           all_mother_pdg[kMaxall];   //[all_]
	
	  Int_t           all_mother_proc[kMaxall];   //[all_]
	
	  Double_t        all_mother_momentum_x[kMaxall];   //[all_]
	
	  Double_t        all_mother_momentum_y[kMaxall];   //[all_]
	
	  Double_t        all_mother_momentum_z[kMaxall];   //[all_]
	
	  Double_t        all_mother_ek[kMaxall];   //[all_]
	
	  Int_t           all_his_nqel[kMaxall];   //[all_]
	
	  Int_t           all_his_nspp[kMaxall];   //[all_]
	
	  Int_t           all_his_ndpp[kMaxall];   //[all_]
	
	  Int_t           all_his_pqel[kMaxall];   //[all_]
	
	  Int_t           all_his_pcex[kMaxall];   //[all_]
	
	  Int_t           all_his_pspp[kMaxall];   //[all_]
	
	  Int_t           all_his_pdpp[kMaxall];   //[all_]
	
	  Int_t           all_his_ptpp[kMaxall];   //[all_]
	
	  Double_t        weight;
	
	  Double_t        norm;
	
	  Double_t        r_x;
	
	  Double_t        r_y;
	
	  Double_t        r_z;
	
	  Double_t        density;
	
	  Int_t           dyn;
	
	  Int_t           nod[12];
	
	  Int_t           place[11][20];
	
	  Double_t        pabsen;
	
	  Bool_t          nopp;
	
	  Bool_t          abs;
	
	  Int_t           nofi;
	
	  Double_t        pen[10];
	
	  Double_t        absr;
	
	  Double_t        odl;
	
	  Double_t        density_hist[50];
	
	  Double_t        radius_hist;
	
	  Int_t           protons_hist;
	
	  Int_t           neutrons_hist;
	
	  Int_t           pr;
	
	  Int_t           nr;
	
	
	  // List of branches
	
	  TBranch        *b_e_fUniqueID;   //!
	
	  TBranch        *b_e_fBits;   //!
	
	  TBranch        *b_e_flag_coh;   //!
	
	  TBranch        *b_e_flag_qel;   //!
	
	  TBranch        *b_e_flag_dis;   //!
	
	  TBranch        *b_e_flag_res;   //!
	
	  TBranch        *b_e_flag_nc;   //!
	
	  TBranch        *b_e_flag_cc;   //!
	
	  TBranch        *b_e_flag_anty;   //!
	
	  TBranch        *b_e_par_random_seed;   //!
	
	  TBranch        *b_e_par_number_of_events;   //!
	
	  TBranch        *b_e_par_number_of_test_events;   //!
	
	  TBranch        *b_e_par_user_events;   //!
	
	  TBranch        *b_e_par_beam_type;   //!
	
	  TBranch        *b_e_par_beam_energy_string;   //!
	
	  TBranch        *b_e_par_beam_particle;   //!
	
	  TBranch        *b_e_par_beam_direction_x;   //!
	
	  TBranch        *b_e_par_beam_direction_y;   //!
	
	  TBranch        *b_e_par_beam_direction_z;   //!
	
	  TBranch        *b_e_par_beam_content_string;   //!
	
	  TBranch        *b_e_par_beam_folder;   //!
	
	  TBranch        *b_e_par_beam_file_first;   //!
	
	  TBranch        *b_e_par_beam_file_limit;   //!
	
	  TBranch        *b_e_par_beam_weighted;   //!
	
	  TBranch        *b_e_par_beam_file;   //!
	
	  TBranch        *b_e_par_beam_offset_x;   //!
	
	  TBranch        *b_e_par_beam_offset_y;   //!
	
	  TBranch        *b_e_par_beam_offset_z;   //!
	
	  TBranch        *b_e_par_beam_placement;   //!
	
	  TBranch        *b_e_par_beam_test_only;   //!
	
	  TBranch        *b_e_par_nucleus_p;   //!
	
	  TBranch        *b_e_par_nucleus_n;   //!
	
	  TBranch        *b_e_par_nucleus_density;   //!
	
	  TBranch        *b_e_par_nucleus_E_b;   //!
	
	  TBranch        *b_e_par_nucleus_kf;   //!
	
	  TBranch        *b_e_par_nucleus_target;   //!
	
	  TBranch        *b_e_par_nucleus_model;   //!
	
	  TBranch        *b_e_par_target_type;   //!
	
	  TBranch        *b_e_par_target_content_string;   //!
	
	  TBranch        *b_e_par_geo_file;   //!
	
	  TBranch        *b_e_par_geo_name;   //!
	
	  TBranch        *b_e_par_geo_volume;   //!
	
	  TBranch        *b_e_par_geo_o_x;   //!
	
	  TBranch        *b_e_par_geo_o_y;   //!
	
	  TBranch        *b_e_par_geo_o_z;   //!
	
	  TBranch        *b_e_par_geo_d_x;   //!
	
	  TBranch        *b_e_par_geo_d_y;   //!
	
	  TBranch        *b_e_par_geo_d_z;   //!
	
	  TBranch        *b_e_par_dyn_qel_cc;   //!
	
	  TBranch        *b_e_par_dyn_qel_nc;   //!
	
	  TBranch        *b_e_par_dyn_res_cc;   //!
	
	  TBranch        *b_e_par_dyn_res_nc;   //!
	
	  TBranch        *b_e_par_dyn_dis_cc;   //!
	
	  TBranch        *b_e_par_dyn_dis_nc;   //!
	
	  TBranch        *b_e_par_dyn_coh_cc;   //!
	
	  TBranch        *b_e_par_dyn_coh_nc;   //!
	
	  TBranch        *b_e_par_qel_vector_ff_set;   //!
	
	  TBranch        *b_e_par_qel_axial_ff_set;   //!
	
	  TBranch        *b_e_par_qel_cc_vector_mass;   //!
	
	  TBranch        *b_e_par_qel_cc_axial_mass;   //!
	
	  TBranch        *b_e_par_qel_nc_axial_mass;   //!
	
	  TBranch        *b_e_par_qel_s_axial_mass;   //!
	
	  TBranch        *b_e_par_qel_strange;   //!
	
	  TBranch        *b_e_par_qel_strangeEM;   //!
	
	  TBranch        *b_e_par_delta_s;   //!
	
	  TBranch        *b_e_par_flux_correction;   //!
	
	  TBranch        *b_e_par_delta_FF_set;   //!
	
	  TBranch        *b_e_par_pion_axial_mass;   //!
	
	  TBranch        *b_e_par_pion_C5A;   //!
	
	  TBranch        *b_e_par_res_dis_cut;   //!
	
	  TBranch        *b_e_par_spp_precision;   //!
	
	  TBranch        *b_e_par_coh_mass_correction;   //!
	
	  TBranch        *b_e_par_coh_new;   //!
	
	  TBranch        *b_e_par_kaskada_on;   //!
	
	  TBranch        *b_e_par_kaskada_newangle;   //!
	
	  TBranch        *b_e_par_kaskada_redo;   //!
	
	  TBranch        *b_e_par_pauli_blocking;   //!
	
	  TBranch        *b_e_par_formation_zone;   //!
	
	  TBranch        *b_e_par_first_step;   //!
	
	  TBranch        *b_e_par_step;   //!
	
	  TBranch        *b_e_par_xsec;   //!
	
	  TBranch        *b_e_par_sf_method;   //!
	
	  TBranch        *b_e_par_cc_smoothing;   //!
	
	  TBranch        *b_e_par_mixed_order;   //!
	
	  TBranch        *b_e_par_rmin;   //!
	
	  TBranch        *b_e_par_rmax;   //!
	
	  TBranch        *b_e_par_path_to_data;   //!
	
	  TBranch        *b_e_in_;   //!
	
	  TBranch        *b_in_t;   //!
	
	  TBranch        *b_in_x;   //!
	
	  TBranch        *b_in_y;   //!
	
	  TBranch        *b_in_z;   //!
	
	  TBranch        *b_in__mass;   //!
	
	  TBranch        *b_in_r_t;   //!
	
	  TBranch        *b_in_r_x;   //!
	
	  TBranch        *b_in_r_y;   //!
	
	  TBranch        *b_in_r_z;   //!
	
	  TBranch        *b_in_pdg;   //!
	
	  TBranch        *b_in_ks;   //!
	
	  TBranch        *b_in_orgin;   //!
	
	  TBranch        *b_in_travelled;   //!
	
	  TBranch        *b_in_id;   //!
	
	  TBranch        *b_in_mother;   //!
	
	  TBranch        *b_in_endproc;   //!
	
	  TBranch        *b_in_fz;   //!
	
	  TBranch        *b_in_pt;   //!
	
	  TBranch        *b_in_mother_pdg;   //!
	
	  TBranch        *b_in_mother_proc;   //!
	
	  TBranch        *b_in_mother_momentum_x;   //!
	
	  TBranch        *b_in_mother_momentum_y;   //!
	
	  TBranch        *b_in_mother_momentum_z;   //!
	
	  TBranch        *b_in_mother_ek;   //!
	
	  TBranch        *b_in_his_nqel;   //!
	
	  TBranch        *b_in_his_nspp;   //!
	
	  TBranch        *b_in_his_ndpp;   //!
	
	  TBranch        *b_in_his_pqel;   //!
	
	  TBranch        *b_in_his_pcex;   //!
	
	  TBranch        *b_in_his_pspp;   //!
	
	  TBranch        *b_in_his_pdpp;   //!
	
	  TBranch        *b_in_his_ptpp;   //!
	
	  TBranch        *b_e_temp_;   //!
	
	  TBranch        *b_temp_t;   //!
	
	  TBranch        *b_temp_x;   //!
	
	  TBranch        *b_temp_y;   //!
	
	  TBranch        *b_temp_z;   //!
	
	  TBranch        *b_temp__mass;   //!
	
	  TBranch        *b_temp_r_t;   //!
	
	  TBranch        *b_temp_r_x;   //!
	
	  TBranch        *b_temp_r_y;   //!
	
	  TBranch        *b_temp_r_z;   //!
	
	  TBranch        *b_temp_pdg;   //!
	
	  TBranch        *b_temp_ks;   //!
	
	  TBranch        *b_temp_orgin;   //!
	
	  TBranch        *b_temp_travelled;   //!
	
	  TBranch        *b_temp_id;   //!
	
	  TBranch        *b_temp_mother;   //!
	
	  TBranch        *b_temp_endproc;   //!
	
	  TBranch        *b_temp_fz;   //!
	
	  TBranch        *b_temp_pt;   //!
	
	  TBranch        *b_temp_mother_pdg;   //!
	
	  TBranch        *b_temp_mother_proc;   //!
	
	  TBranch        *b_temp_mother_momentum_x;   //!
	
	  TBranch        *b_temp_mother_momentum_y;   //!
	
	  TBranch        *b_temp_mother_momentum_z;   //!
	
	  TBranch        *b_temp_mother_ek;   //!
	
	  TBranch        *b_temp_his_nqel;   //!
	
	  TBranch        *b_temp_his_nspp;   //!
	
	  TBranch        *b_temp_his_ndpp;   //!
	
	  TBranch        *b_temp_his_pqel;   //!
	
	  TBranch        *b_temp_his_pcex;   //!
	
	  TBranch        *b_temp_his_pspp;   //!
	
	  TBranch        *b_temp_his_pdpp;   //!
	
	  TBranch        *b_temp_his_ptpp;   //!
	
	  TBranch        *b_e_out_;   //!
	
	  TBranch        *b_out_t;   //!
	
	  TBranch        *b_out_x;   //!
	
	  TBranch        *b_out_y;   //!
	
	  TBranch        *b_out_z;   //!
	
	  TBranch        *b_out__mass;   //!
	
	  TBranch        *b_out_r_t;   //!
	
	  TBranch        *b_out_r_x;   //!
	
	  TBranch        *b_out_r_y;   //!
	
	  TBranch        *b_out_r_z;   //!
	
	  TBranch        *b_out_pdg;   //!
	
	  TBranch        *b_out_ks;   //!
	
	  TBranch        *b_out_orgin;   //!
	
	  TBranch        *b_out_travelled;   //!
	
	  TBranch        *b_out_id;   //!
	
	  TBranch        *b_out_mother;   //!
	
	  TBranch        *b_out_endproc;   //!
	
	  TBranch        *b_out_fz;   //!
	
	  TBranch        *b_out_pt;   //!
	
	  TBranch        *b_out_mother_pdg;   //!
	
	  TBranch        *b_out_mother_proc;   //!
	
	  TBranch        *b_out_mother_momentum_x;   //!
	
	  TBranch        *b_out_mother_momentum_y;   //!
	
	  TBranch        *b_out_mother_momentum_z;   //!
	
	  TBranch        *b_out_mother_ek;   //!
	
	  TBranch        *b_out_his_nqel;   //!
	
	  TBranch        *b_out_his_nspp;   //!
	
	  TBranch        *b_out_his_ndpp;   //!
	
	  TBranch        *b_out_his_pqel;   //!
	
	  TBranch        *b_out_his_pcex;   //!
	
	  TBranch        *b_out_his_pspp;   //!
	
	  TBranch        *b_out_his_pdpp;   //!
	
	  TBranch        *b_out_his_ptpp;   //!
	
	  TBranch        *b_e_post_;   //!
	
	  TBranch        *b_post_t;   //!
	
	  TBranch        *b_post_x;   //!
	
	  TBranch        *b_post_y;   //!
	
	  TBranch        *b_post_z;   //!
	
	  TBranch        *b_post__mass;   //!
	
	  TBranch        *b_post_r_t;   //!
	
	  TBranch        *b_post_r_x;   //!
	
	  TBranch        *b_post_r_y;   //!
	
	  TBranch        *b_post_r_z;   //!
	
	  TBranch        *b_post_pdg;   //!
	
	  TBranch        *b_post_ks;   //!
	
	  TBranch        *b_post_orgin;   //!
	
	  TBranch        *b_post_travelled;   //!
	
	  TBranch        *b_post_id;   //!
	
	  TBranch        *b_post_mother;   //!
	
	  TBranch        *b_post_endproc;   //!
	
	  TBranch        *b_post_fz;   //!
	
	  TBranch        *b_post_pt;   //!
	
	  TBranch        *b_post_mother_pdg;   //!
	
	  TBranch        *b_post_mother_proc;   //!
	
	  TBranch        *b_post_mother_momentum_x;   //!
	
	  TBranch        *b_post_mother_momentum_y;   //!
	
	  TBranch        *b_post_mother_momentum_z;   //!
	
	  TBranch        *b_post_mother_ek;   //!
	
	  TBranch        *b_post_his_nqel;   //!
	
	  TBranch        *b_post_his_nspp;   //!
	
	  TBranch        *b_post_his_ndpp;   //!
	
	  TBranch        *b_post_his_pqel;   //!
	
	  TBranch        *b_post_his_pcex;   //!
	
	  TBranch        *b_post_his_pspp;   //!
	
	  TBranch        *b_post_his_pdpp;   //!
	
	  TBranch        *b_post_his_ptpp;   //!
	
	  TBranch        *b_e_all_;   //!
	
	  TBranch        *b_all_t;   //!
	
	  TBranch        *b_all_x;   //!
	
	  TBranch        *b_all_y;   //!
	
	  TBranch        *b_all_z;   //!
	
	  TBranch        *b_all__mass;   //!
	
	  TBranch        *b_all_r_t;   //!
	
	  TBranch        *b_all_r_x;   //!
	
	  TBranch        *b_all_r_y;   //!
	
	  TBranch        *b_all_r_z;   //!
	
	  TBranch        *b_all_pdg;   //!
	
	  TBranch        *b_all_ks;   //!
	
	  TBranch        *b_all_orgin;   //!
	
	  TBranch        *b_all_travelled;   //!
	
	  TBranch        *b_all_id;   //!
	
	  TBranch        *b_all_mother;   //!
	
	  TBranch        *b_all_endproc;   //!
	
	  TBranch        *b_all_fz;   //!
	
	  TBranch        *b_all_pt;   //!
	
	  TBranch        *b_all_mother_pdg;   //!
	
	  TBranch        *b_all_mother_proc;   //!
	
	  TBranch        *b_all_mother_momentum_x;   //!
	
	  TBranch        *b_all_mother_momentum_y;   //!
	
	  TBranch        *b_all_mother_momentum_z;   //!
	
	  TBranch        *b_all_mother_ek;   //!
	
	  TBranch        *b_all_his_nqel;   //!
	
	  TBranch        *b_all_his_nspp;   //!
	
	  TBranch        *b_all_his_ndpp;   //!
	
	  TBranch        *b_all_his_pqel;   //!
	
	  TBranch        *b_all_his_pcex;   //!
	
	  TBranch        *b_all_his_pspp;   //!
	
	  TBranch        *b_all_his_pdpp;   //!
	
	  TBranch        *b_all_his_ptpp;   //!
	
	  TBranch        *b_e_weight;   //!
	
	  TBranch        *b_e_norm;   //!
	
	  TBranch        *b_e_r_x;   //!
	
	  TBranch        *b_e_r_y;   //!
	
	  TBranch        *b_e_r_z;   //!
	
	  TBranch        *b_e_density;   //!
	
	  TBranch        *b_e_dyn;   //!
	
	  TBranch        *b_e_nod;   //!
	
	  TBranch        *b_e_place;   //!
	
	  TBranch        *b_e_pabsen;   //!
	
	  TBranch        *b_e_nopp;   //!
	
	  TBranch        *b_e_abs;   //!
	
	  TBranch        *b_e_nofi;   //!
	
	  TBranch        *b_e_pen;   //!
	
	  TBranch        *b_e_absr;   //!
	
	  TBranch        *b_e_odl;   //!
	
	  TBranch        *b_e_density_hist;   //!
	
	  TBranch        *b_e_radius_hist;   //!
	
	  TBranch        *b_e_protons_hist;   //!
	
	  TBranch        *b_e_neutrons_hist;   //!
	
	  TBranch        *b_e_pr;   //!
	
	  TBranch        *b_e_nr;   //! 
	  */
	
  };
}


namespace evgen{

  //____________________________________________________________________________
  NuWroGen::NuWroGen(fhicl::ParameterSet const& pset)
  {

    fEventNumberOffset = pset.get< int >("EventNumberOffset",0);

    this->reconfigure(pset); 
    fStopwatch.Start();

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

   }

  //____________________________________________________________________________
  
  NuWroGen::~NuWroGen()
  {  
    fStopwatch.Stop();
  }

  //____________________________________________________________________________
  void NuWroGen::reconfigure(fhicl::ParameterSet const& p)
  {
    fFileName       = p.get< std::string         >("NuWroFile","output.root");
    fTreeName       = p.get< std::string   	  >("TreeName","treeout"); 

    return;
  }
//___________________________________________________________________________
  
  void NuWroGen::beginJob(){
   
    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;

    fGenerated[0] = tfs->make<TH1F>("fGenerated_necc","",  100, 0.0, 20.0);
    fGenerated[1] = tfs->make<TH1F>("fGenerated_nebcc","", 100, 0.0, 20.0);
    fGenerated[2] = tfs->make<TH1F>("fGenerated_nmcc","",  100, 0.0, 20.0);
    fGenerated[3] = tfs->make<TH1F>("fGenerated_nmbcc","", 100, 0.0, 20.0);
    fGenerated[4] = tfs->make<TH1F>("fGenerated_nnc","",   100, 0.0, 20.0);
    fGenerated[5] = tfs->make<TH1F>("fGenerated_nbnc","",  100, 0.0, 20.0);
    
    fDCosX = tfs->make<TH1F>("fDCosX", ";dx/ds", 200, -1., 1.);
    fDCosY = tfs->make<TH1F>("fDCosY", ";dy/ds", 200, -1., 1.);
    fDCosZ = tfs->make<TH1F>("fDCosZ", ";dz/ds", 200, -1., 1.);

    fMuMomentum = tfs->make<TH1F>("fMuMomentum", ";p_{#mu} (GeV/c)", 500, 0., 50.);
    fMuDCosX    = tfs->make<TH1F>("fMuDCosX", ";dx/ds;", 200, -1., 1.);
    fMuDCosY    = tfs->make<TH1F>("fMuDCosY", ";dy/ds;", 200, -1., 1.);
    fMuDCosZ    = tfs->make<TH1F>("fMuDCosZ", ";dz/ds;", 200, -1., 1.);

    fEMomentum  = tfs->make<TH1F>("fEMomentum", ";p_{e} (GeV/c)", 500, 0., 50.);
    fEDCosX     = tfs->make<TH1F>("fEDCosX", ";dx/ds;", 200, -1., 1.);
    fEDCosY     = tfs->make<TH1F>("fEDCosY", ";dy/ds;", 200, -1., 1.);
    fEDCosZ     = tfs->make<TH1F>("fEDCosZ", ";dz/ds;", 200, -1., 1.);

    fCCMode = tfs->make<TH1F>("fCCMode", ";CC Interaction Mode;", 5, 0., 5.);
    fCCMode->GetXaxis()->SetBinLabel(1, "QE");
    fCCMode->GetXaxis()->SetBinLabel(2, "Res");
    fCCMode->GetXaxis()->SetBinLabel(3, "DIS");
    fCCMode->GetXaxis()->SetBinLabel(4, "Coh");
    fCCMode->GetXaxis()->SetBinLabel(5, "kInverseMuDecay");
    fCCMode->GetXaxis()->CenterLabels();

    fNCMode = tfs->make<TH1F>("fNCMode", ";NC Interaction Mode;", 5, 0., 5.);
    fNCMode->GetXaxis()->SetBinLabel(1, "QE");
    fNCMode->GetXaxis()->SetBinLabel(2, "Res");
    fNCMode->GetXaxis()->SetBinLabel(3, "DIS");
    fNCMode->GetXaxis()->SetBinLabel(4, "Coh");
    fNCMode->GetXaxis()->SetBinLabel(5, "kNuElectronElastic");
    fNCMode->GetXaxis()->CenterLabels();

    fWeight = tfs->make<TH1F>("fWeight", ";Weight1;", 20, 1.E-42, 1.E-38);
    fWeightNW = tfs->make<TH1F>("fWeightNW", ";Weight2;", 20, 1.E-42, 1.E-38);
    
    fDyn = tfs->make<TH1F>("fDyn", ";Canonical Interaction Mode;", 13, 0., 13.);
    fDyn->GetXaxis()->SetBinLabel(1, "CCQE");
    fDyn->GetXaxis()->SetBinLabel(2, "NCelastic");
    fDyn->GetXaxis()->SetBinLabel(3, "CCResp2ppi+");
    fDyn->GetXaxis()->SetBinLabel(4, "CCResn2ppi0");
    fDyn->GetXaxis()->SetBinLabel(5, "CCResn2npi+");
    fDyn->GetXaxis()->SetBinLabel(6, "NCResp2ppi0");
    fDyn->GetXaxis()->SetBinLabel(7, "NCResp2npi+");
    fDyn->GetXaxis()->SetBinLabel(8, "NCResn2npi0");
    fDyn->GetXaxis()->SetBinLabel(9, "NCResn2ppi-");
    fDyn->GetXaxis()->SetBinLabel(10, "CC-DIS");
    fDyn->GetXaxis()->SetBinLabel(11, "NC-DIS");
    fDyn->GetXaxis()->SetBinLabel(12, "NC-COH");
    fDyn->GetXaxis()->SetBinLabel(13, "CC-COH");
    fDyn->GetXaxis()->CenterLabels();

    fDynNew = tfs->make<TH1F>("fDynNew", ";New Style Accounting Mode;", 10, 0., 10.);
    fDynNew->GetXaxis()->SetBinLabel(1, "1mu0p0pi");
    fDynNew->GetXaxis()->SetBinLabel(2, "1mu1p0pi");
    fDynNew->GetXaxis()->SetBinLabel(3, "1muge2p0pi");
    fDynNew->GetXaxis()->SetBinLabel(4, "1mu0p1pi");
    fDynNew->GetXaxis()->SetBinLabel(5, "1mu1p1pi");
    fDynNew->GetXaxis()->SetBinLabel(6, "1muge2p1pi");
    fDynNew->GetXaxis()->SetBinLabel(7, "1mu0p2pi");
    fDynNew->GetXaxis()->SetBinLabel(8, "1muge1p2pi");
    fDynNew->GetXaxis()->SetBinLabel(9, "NC");
    fDynNew->GetXaxis()->SetBinLabel(10, "Other");
    fDynNew->GetXaxis()->CenterLabels();

    fDynNewThresh = tfs->make<TH1F>("fDynNewThresh", ";New Style Accounting Mode (Tp>50MeV);", 10, 0., 10.);
    fDynNewThresh->GetXaxis()->SetBinLabel(1, "1mu0p0pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(2, "1mu1p0pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(3, "1muge2p0pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(4, "1mu0p1pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(5, "1mu1p1pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(6, "1muge2p1pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(7, "1mu0p2pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(8, "1muge1p2pi");
    fDynNewThresh->GetXaxis()->SetBinLabel(9, "NC");
    fDynNewThresh->GetXaxis()->SetBinLabel(10, "Other");
    fDynNewThresh->GetXaxis()->CenterLabels();

    f2DynNew = tfs->make<TH2F>("f2DynNew", ";Old vs New Style Accounting Mode;", 13, 0.,13., 10,0.,10.);
 
    f2DynNew->GetXaxis()->SetBinLabel(1, "CCQE");
    f2DynNew->GetXaxis()->SetBinLabel(2, "NCelastic");
    f2DynNew->GetXaxis()->SetBinLabel(3, "CCResp2ppi+");
    f2DynNew->GetXaxis()->SetBinLabel(4, "CCResn2ppi0");
    f2DynNew->GetXaxis()->SetBinLabel(5, "CCResn2npi+");
    f2DynNew->GetXaxis()->SetBinLabel(6, "NCResp2ppi0");
    f2DynNew->GetXaxis()->SetBinLabel(7, "NCResp2npi+");
    f2DynNew->GetXaxis()->SetBinLabel(8, "NCResn2npi0");
    f2DynNew->GetXaxis()->SetBinLabel(9, "NCResn2ppi-");
    f2DynNew->GetXaxis()->SetBinLabel(10, "CC-DIS");
    f2DynNew->GetXaxis()->SetBinLabel(11, "NC-DIS");
    f2DynNew->GetXaxis()->SetBinLabel(12, "NC-COH");
    f2DynNew->GetXaxis()->SetBinLabel(13, "CC-COH");
    f2DynNew->GetXaxis()->CenterLabels();
    f2DynNew->GetYaxis()->SetBinLabel(1, "1mu0p0pi");
    f2DynNew->GetYaxis()->SetBinLabel(2, "1mu1p0pi");
    f2DynNew->GetYaxis()->SetBinLabel(3, "1muge2p0pi");
    f2DynNew->GetYaxis()->SetBinLabel(4, "1mu0p1pi");
    f2DynNew->GetYaxis()->SetBinLabel(5, "1mu1p1pi");
    f2DynNew->GetYaxis()->SetBinLabel(6, "1muge2p1pi");
    f2DynNew->GetYaxis()->SetBinLabel(7, "1mu0p2pi");
    f2DynNew->GetYaxis()->SetBinLabel(8, "1muge1p2pi");
    f2DynNew->GetYaxis()->SetBinLabel(9, "NC");
    f2DynNew->GetYaxis()->SetBinLabel(10, "Other");
    f2DynNew->GetYaxis()->CenterLabels();

    f2DynNewThresh = tfs->make<TH2F>("f2DynNewThresh", ";Old vs New Style Accounting Mode;", 13, 0.,13., 10, 0.,10.);
 
    f2DynNewThresh->GetXaxis()->SetBinLabel(1, "CCQE");
    f2DynNewThresh->GetXaxis()->SetBinLabel(2, "NCelastic");
    f2DynNewThresh->GetXaxis()->SetBinLabel(3, "CCResp2ppi+");
    f2DynNewThresh->GetXaxis()->SetBinLabel(4, "CCResn2ppi0");
    f2DynNewThresh->GetXaxis()->SetBinLabel(5, "CCResn2npi+");
    f2DynNewThresh->GetXaxis()->SetBinLabel(6, "NCResp2ppi0");
    f2DynNewThresh->GetXaxis()->SetBinLabel(7, "NCResp2npi+");
    f2DynNewThresh->GetXaxis()->SetBinLabel(8, "NCResn2npi0");
    f2DynNewThresh->GetXaxis()->SetBinLabel(9, "NCResn2ppi-");
    f2DynNewThresh->GetXaxis()->SetBinLabel(10, "CC-DIS");
    f2DynNewThresh->GetXaxis()->SetBinLabel(11, "NC-DIS");
    f2DynNewThresh->GetXaxis()->SetBinLabel(12, "NC-COH");
    f2DynNewThresh->GetXaxis()->SetBinLabel(13, "CC-COH");
    f2DynNewThresh->GetXaxis()->CenterLabels();
    f2DynNewThresh->GetYaxis()->SetBinLabel(1, "1mu0p0pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(2, "1mu1p0pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(3, "1muge2p0pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(4, "1mu0p1pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(5, "1mu1p1pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(6, "1muge2p1pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(7, "1mu0p2pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(8, "1muge1p2pi");
    f2DynNewThresh->GetYaxis()->SetBinLabel(9, "NC");
    f2DynNewThresh->GetYaxis()->SetBinLabel(10, "Other");
    f2DynNewThresh->GetYaxis()->CenterLabels();
 
    //fDeltaE = tfs->make<TH1F>("fDeltaE", ";#Delta E_{#nu} (GeV);", 200, -1., 1.); 
    fECons  = tfs->make<TH1F>("fECons", ";#Delta E(#nu,lepton);", 500, -5., 5.);

    art::ServiceHandle<geo::Geometry> geo;
    double x = 2.1*geo->DetHalfWidth();
    double y = 2.1*geo->DetHalfHeight();
    double z = 2.*geo->DetLength();
    int xdiv = TMath::Nint(2*x/5.);
    int ydiv = TMath::Nint(2*y/5.);
    int zdiv = TMath::Nint(2*z/5.);
 
    fVertexX = tfs->make<TH1F>("fVertexX", ";x (cm)", xdiv,  -x, x);
    fVertexY = tfs->make<TH1F>("fVertexY", ";y (cm)", ydiv,  -y, y);
    fVertexZ = tfs->make<TH1F>("fVertexZ", ";z (cm)", zdiv, -0.2*z, z);
 
    fVertexXY = tfs->make<TH2F>("fVertexXY", ";x (cm);y (cm)", xdiv,     -x, x, ydiv, -y, y);
    fVertexXZ = tfs->make<TH2F>("fVertexXZ", ";z (cm);x (cm)", zdiv, -0.2*z, z, xdiv, -x, x);
    fVertexYZ = tfs->make<TH2F>("fVertexYZ", ";z (cm);y (cm)", zdiv, -0.2*z, z, ydiv, -y, y);


    TClass::GetClass("line")->GetStreamerInfo(1);
    ch = new TChain(fTreeName.c_str());
    ch->Add(fFileName.c_str());


    std::cout << " Num entries in TTree is " << ch->GetEntries() << std::endl;

    NuWroTTree = new event();
    ch->SetBranchAddress ("e", &NuWroTTree);

    // initiate flux-wtd XSections.
    std::cout << "NuWroGen: Here's the output of the .txt file" << std::endl;
    ifstream xsecTxtFile((fFileName+".txt").c_str());
    unsigned int cntline(0);
    std::string line;
    if (xsecTxtFile.is_open())
      {
	while ( xsecTxtFile.good() )
	  {
	    getline (xsecTxtFile,line);
	    cout << line << endl;
	    if (cntline==0) { cntline++; continue;} // first line is header.
	    stringstream ss(line); // Insert the string into a stream
	    vector<std::string> tokens; // Create vector to hold our words
	    string buf;
	    while (ss >> buf)
	      {
		tokens.push_back(buf);
	      }
	    // want last element in line. That's the xsec.
	    if (tokens.size() && line.length())
	      fxsecFluxWtd.push_back(atof(tokens.back().c_str())); 
	    else  xsecTxtFile.close();
	    tokens.clear();
	  }
      }
    else std::cout << "Unable to open file"; 

    fxsecTotal=0.0;
    for (unsigned int ii=0;ii<fxsecFluxWtd.size();ii++)
      {
	fxsecTotal+=fxsecFluxWtd.at(ii);
      }

    countFile = fEventNumberOffset;
    std::cout << " Start this job on event " << countFile << std::endl;
  }

  //____________________________________________________________________________
  void NuWroGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  void NuWroGen::endJob()
  {
    delete NuWroTTree;
    delete ch;
    fxsecFluxWtd.clear();
  }
  //____________________________________________________________________________
  void NuWroGen::produce(art::Event& evt)
  {

    std::cout << std::endl;
    std::cout<<"------------------------------------------------------------------------------"<<std::endl;
//  std::cout << "run    : " << evt.Header().Run() << std::endl;
//  std::cout << "subrun : " << evt.Header().Subrun() << std::endl;
//  std::cout << "event  : " << evt.Header().Event() << std::endl;
//  std::cout << "event  : " << evt.id().event() << std::endl;  
    
    std::string name, k, dollar;
    int partnumber = 0;
    
    int trackid = -1; // set track id to -i as these are all primary particles and have id <= 0
    std::string primary("primary");
    int FirstMother = -1;

    int Status = -9999;
    

    int ccnc = -9999;
    int mode = -9999;
    int targetnucleusPdg = -9999;
    int hitquarkPdg = -9999;

    TLorentzVector Neutrino;
    TLorentzVector Lepton;
    TLorentzVector Target;
    TLorentzVector q;
    TLorentzVector Hadron4mom;
    double Q2 = -9999;

    int Tpdg = 0;  // for target 
    double Tmass = 0;
    int Tstatus = 11;
    double Tcosx, Tcosy, Tcosz, Tenergy;
    TLorentzVector Tpos;
    
        
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
    simb::MCTruth truth;
    
    /*    
    if (countFile>=ch->GetEntries()) 
      {
	mf::LogInfo("NuWro (parser): You're on row ")  << countFile << " of " << ch->GetEntries()<< ". Moving on ..." <<std::endl;
	return;
      }
    */
    ch->GetEntry(countFile++);
      
    //get the nuwro channel number
	
    //set the interaction type; CC or NC       
    if (NuWroTTree->flag.cc) ccnc = simb::kCC;
    else if (NuWroTTree->flag.nc) ccnc = simb::kNC;
	  
    //set the interaction mode; QE, Res, DIS, Coh, kNuElectronElastic, kInverseMuDecay 
    if ( NuWroTTree->flag.qel )
      mode = simb::kQE; 
    else if ( NuWroTTree->flag.res )
      mode = simb::kRes;
    else if ( NuWroTTree->flag.dis )
      mode = simb::kDIS;
    else if ( NuWroTTree->flag.coh )
      mode = simb::kCoh;
    if(partnumber == -1)
      Status = 0;
    else 
      Status = 1;


    art::ServiceHandle<geo::Geometry> geo;
    double X0 = NuWroTTree->par.geo_o[0] + geo->DetHalfWidth();
    double Y0 = NuWroTTree->par.geo_o[1];
    double Z0 = NuWroTTree->par.geo_o[2] + 0.25*geo->DetLength();
    TLorentzVector pos(X0, Y0, Z0, 0);
    Tpos = pos; // for target


    simb::MCParticle mcpartNu(trackid,
			    NuWroTTree->in[0].pdg,
			    primary,		    
			    FirstMother,
			    NuWroTTree->in[0].m()/1000.,
			    Status
			    );
    TLorentzVector mom(NuWroTTree->in[0].x/1000., 
		       NuWroTTree->in[0].y/1000.,
		       NuWroTTree->in[0].z/1000., 
		       NuWroTTree->in[0].t/1000.);
    mcpartNu.AddTrajectoryPoint(pos,mom);
    truth.Add(mcpartNu);

           
    unsigned int ii(0);
    while(ii<NuWroTTree->post.size())
      {
	// loop over particles in an event	
	    
	simb::MCParticle mcpart(trackid,
				NuWroTTree->post[ii].pdg,
				primary,		    
				FirstMother,
				NuWroTTree->post[ii].m()/1000.,
				Status
				);
	TLorentzVector mom(NuWroTTree->post[ii].x/1000., 
			   NuWroTTree->post[ii].y/1000.,
			   NuWroTTree->post[ii].z/1000., 
			   NuWroTTree->post[ii].t/1000.);
	mcpart.AddTrajectoryPoint(pos,mom);
	truth.Add(mcpart);
	
		
	ii++;
      }// loop over particles in an event
    // Incoming Neutrino is 0th element of in. Outgoing lepton is 0th of out.
    Neutrino.SetPxPyPzE(NuWroTTree->in[0].x/1000., NuWroTTree->in[0].y/1000., NuWroTTree->in[0].z/1000., NuWroTTree->in[0].t/1000. );
    Lepton.SetPxPyPzE(NuWroTTree->out[0].x/1000., NuWroTTree->out[0].y/1000., NuWroTTree->out[0].z/1000., NuWroTTree->out[0].t/1000. );
    
    /////////////////////////////////
    
      Tmass = NuWroTTree->par.nucleus_p + NuWroTTree->par.nucleus_n; // GeV
	
      
      Tenergy = NuWroTTree->in.back().t;
      Tcosx   = NuWroTTree->in.back().x;
      Tcosy   = NuWroTTree->in.back().y;
      Tcosz   = NuWroTTree->in.back().z;
      Tmass = std::sqrt(std::abs(Tenergy*Tenergy - Tcosx*Tcosx 
		       - Tcosy*Tcosy - Tcosz*Tcosz))/1000.;

      Tenergy = Tmass-0.1; // force this negative, cuz kinetic energy>eps
      // seems to make G4 hang.
      Tcosx =0.; Tcosy = 0.; Tcosz = 0.;


      simb::MCParticle mcpartT(trackid,
			      Tpdg,
			      primary,		    
			      FirstMother,
			      Tmass,
			      Tstatus
			      );
      

      //     Target = Hadron4mom - (Neutrino - Lepton); // commenting this out as target momentum no more is calculated by 4-momentum conservation 
    
      TLorentzVector Tmom;
      Tmom.SetPxPyPzE(Tcosx/1000., Tcosy/1000., Tcosz/1000., Tenergy); // this makes literally |P| = 0 or 1 GeV/c for target, this affects only the Kinematic variables; X and Y
      // target |p| = 0 GeV/c if interaction is
      // DIS, Coh, nu-e Elastic Scattering, nu-e inverse mu decay (this comes from Nuwro), 
      // else target |P| = 1 GeV/c
      Target = Tmom;
    

      mcpartT.AddTrajectoryPoint(Tpos,Tmom); 
      // for now, do(n't) target onto stack. EC, 20-Jul-2012.
      truth.Add(mcpartT);
    
      q = Neutrino - Lepton;
      Q2 = -(q * q);
      //    double W2 = Hadron4mom * Hadron4mom;
      //    double InvariantMass = std::sqrt(W2);
    
      double x = Q2/((2*Target*q));
      double y = (Target*q)/(Neutrino*Target);
    
      truth.SetOrigin(simb::kBeamNeutrino);
      int channel(-999);
      targetnucleusPdg = NuWroTTree->in[NuWroTTree->in.size()-1].pdg;
      truth.SetNeutrino(ccnc, mode, channel,
			targetnucleusPdg, 
			Tpdg, 
			hitquarkPdg,
			//InvariantMass, x, y, Q2
			0, x, y, Q2
			);
    
      std::cout << truth.GetNeutrino() << std::endl;
    
      FillHistograms(truth);  

      truthcol->push_back(truth);
      evt.put(std::move(truthcol));
    
    return;
  }
  
//   //......................................................................
  std::string NuWroGen::ParticleStatus(int StatusCode)
  {
    int code = StatusCode;
    std::string ParticleStatusName;

    switch(code)
      {
      case 0:
	ParticleStatusName = "kIStInitialState";
	break;
      case 1:
	ParticleStatusName = "kIStFinalState";
	break;
      case 11:
	ParticleStatusName = "kIStNucleonTarget";
	break;
      default:
	ParticleStatusName = "Status Unknown";
      }
    return ParticleStatusName;
  }


//   //......................................................................
  std::string NuWroGen::ReactionChannel(int ccnc,int mode)
  {
    std::string ReactionChannelName=" ";

    if(ccnc==0)
      ReactionChannelName = "kCC";
    else if(ccnc==1)
      ReactionChannelName = "kNC";
    else std::cout<<"Current mode unknown!! "<<std::endl;

    if(mode==0)
      ReactionChannelName += "_kQE";
    else if(mode==1)
      ReactionChannelName += "_kRes";
    else if(mode==2)
      ReactionChannelName += "_kDIS";
    else if(mode==3)
      ReactionChannelName += "_kCoh";
    else if(mode==4)
      ReactionChannelName += "_kNuElectronElastic";
    else if(mode==5)
      ReactionChannelName += "_kInverseMuDecay";
    else std::cout<<"interaction mode unknown!! "<<std::endl;

    return ReactionChannelName;
  }

//   //......................................................................
  void NuWroGen::FillHistograms(const simb::MCTruth &mc)
  {
    // Decide which histograms to put the spectrum in
    int id = -1;
    if (mc.GetNeutrino().CCNC()==simb::kCC) {
      fCCMode->Fill(mc.GetNeutrino().Mode());
      if      (mc.GetNeutrino().Nu().PdgCode() ==  12) id = 0;
      else if (mc.GetNeutrino().Nu().PdgCode() == -12) id = 1;
      else if (mc.GetNeutrino().Nu().PdgCode() ==  14) id = 2;
      else if (mc.GetNeutrino().Nu().PdgCode() == -14) id = 3;
      else return;
    }
    else {
      fNCMode->Fill(mc.GetNeutrino().Mode());
      if (mc.GetNeutrino().Nu().PdgCode() > 0) id = 4;
      else                                     id = 5;
    }
    if (id==-1) abort();
    
    // Fill the specta histograms
    fGenerated[id]->Fill(mc.GetNeutrino().Nu().E() );
      
    //< fill the vertex histograms from the neutrino - that is always 
    //< particle 0 in the list
    simb::MCNeutrino       mcnu = mc.GetNeutrino();
    const simb::MCParticle nu   = mcnu.Nu();
    
    fVertexX->Fill(nu.Vx());
    fVertexY->Fill(nu.Vy());
    fVertexZ->Fill(nu.Vz());
    
    fVertexXY->Fill(nu.Vx(), nu.Vy());
    fVertexXZ->Fill(nu.Vz(), nu.Vx());
    fVertexYZ->Fill(nu.Vz(), nu.Vy());
    
    double mom = nu.P();
    if(std::abs(mom) > 0.){
      fDCosX->Fill(nu.Px()/mom);
      fDCosY->Fill(nu.Py()/mom);
      fDCosZ->Fill(nu.Pz()/mom);
    }


//     LOG_DEBUG("GENIEInteractionInformation") 
//       << std::endl
//       << "REACTION:  " << ReactionChannel(mc.GetNeutrino().CCNC(),mc.GetNeutrino().Mode()) 
//       << std::endl
//       << "-----------> Particles in the Stack = " << mc.NParticles() << std::endl
//       << std::setiosflags(std::ios::left) 
//       << std::setw(20) << "PARTICLE"
//       << std::setiosflags(std::ios::left) 
//       << std::setw(32) << "STATUS"
//       << std::setw(18) << "E (GeV)"
//       << std::setw(18) << "m (GeV/c2)"
//       << std::setw(18) << "Ek (GeV)"
//       << std::endl << std::endl;

//     const TDatabasePDG* databasePDG = TDatabasePDG::Instance();

//     // Loop over the particle stack for this event 
//     for(int i = 0; i < mc.NParticles(); ++i){
//       simb::MCParticle part(mc.GetParticle(i));
//       std::string name = databasePDG->GetParticle(part.PdgCode())->GetName();
//       int code = part.StatusCode();
//       std::string status = ParticleStatus(code);
//       double mass = part.Mass();
//       double energy = part.E(); 
//       double Ek = (energy-mass); // Kinetic Energy (GeV)
//       if(status=="kIStFinalStB4Interactions")
// 	LOG_DEBUG("GENIEFinalState")
// 	  << std::setiosflags(std::ios::left) << std::setw(20) << name
// 	  << std::setiosflags(std::ios::left) << std::setw(32) <<status
// 	  << std::setw(18)<< energy
// 	  << std::setw(18)<< mass
// 	  << std::setw(18)<< Ek <<std::endl;
//       else 
// 	LOG_DEBUG("GENIEFinalState") 
// 	  << std::setiosflags(std::ios::left) << std::setw(20) << name
// 	  << std::setiosflags(std::ios::left) << std::setw(32) << status
// 	  << std::setw(18) << energy
// 	  << std::setw(18) << mass <<std::endl; 

    std::cout << "REACTION:  " << ReactionChannel(mc.GetNeutrino().CCNC(),mc.GetNeutrino().Mode()) << std::endl;
    std::cout << "-----------> Particles in the Stack = " << mc.NParticles() << std::endl;
    std::cout << std::setiosflags(std::ios::left) 
	      << std::setw(20) << "PARTICLE"
	      << std::setiosflags(std::ios::left) 
	      << std::setw(32) << "STATUS"
	      << std::setw(18) << "E (GeV)"
	      << std::setw(18) << "m (GeV/c2)"
	      << std::setw(18) << "Ek (GeV)"
	      << std::endl << std::endl;
    
    const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
    
    // Loop over the particle stack for this event 
    for(int i = 0; i < mc.NParticles(); ++i){
      simb::MCParticle part(mc.GetParticle(i));
      std::string name;
      if (part.PdgCode() == 18040)
	name = "Ar40 18040";
      else if (part.PdgCode() != -99999 )
	{
	  name = databasePDG->GetParticle(part.PdgCode())->GetName(); 
	}
      
      int code = part.StatusCode();
      std::string status = ParticleStatus(code);
      double mass = part.Mass();
      double energy = part.E(); 
      double Ek = (energy-mass); // Kinetic Energy (GeV)
      
      std::cout << std::setiosflags(std::ios::left) << std::setw(20) << name
		<< std::setiosflags(std::ios::left) << std::setw(32) <<status
		<< std::setw(18)<< energy
		<< std::setw(18)<< mass
		<< std::setw(18)<< Ek <<std::endl;  
    }

    if(mc.GetNeutrino().CCNC() == simb::kCC){
  
      ///look for the outgoing lepton in the particle stack
      ///just interested in the first one
      for(int i = 0; i < mc.NParticles(); ++i){
	simb::MCParticle part(mc.GetParticle(i));
	if(std::abs(part.PdgCode()) == 11){
	  fEMomentum->Fill(part.P());
	  fEDCosX->Fill(part.Px()/part.P());
	  fEDCosY->Fill(part.Py()/part.P());
	  fEDCosZ->Fill(part.Pz()/part.P());
	  fECons->Fill(nu.E() - part.E());
	  break;
	}
	else if(std::abs(part.PdgCode()) == 13){
	  fMuMomentum->Fill(part.P());
	  fMuDCosX->Fill(part.Px()/part.P());
	  fMuDCosY->Fill(part.Py()/part.P());
	  fMuDCosZ->Fill(part.Pz()/part.P());
	  fECons->Fill(nu.E() - part.E());
	  break;
	}
      }// end loop over particles
    }//end if CC interaction

    // fill fDyn
    double bin(0.0);
    double binNew(0.0);
    double binNewThresh(0.0);
    if      (NuWroTTree->flag.qel && NuWroTTree->flag.cc) bin = 1.;
    else if (NuWroTTree->flag.qel && NuWroTTree->flag.nc) bin = 2.;
    else if (NuWroTTree->flag.dis && NuWroTTree->flag.cc) bin = 10.;
    else if (NuWroTTree->flag.dis && NuWroTTree->flag.nc) bin = 11.;
    else if (NuWroTTree->flag.coh && NuWroTTree->flag.cc) bin = 12.;
    else if (NuWroTTree->flag.coh && NuWroTTree->flag.nc) bin = 13.;
    else if (NuWroTTree->flag.res)
      {
	unsigned int ii(0);
	unsigned int p(0), n(0), pip(0), pim(0), pi0(0);
	while(ii<NuWroTTree->post.size())
	  {
	    if  (NuWroTTree->out[ii].pdg==211) pip++;
	    if  (NuWroTTree->out[ii].pdg==-211) pim++;
	    if  (NuWroTTree->out[ii].pdg==111) pi0++;
	    if  (NuWroTTree->out[ii].pdg==2112) n++;
	    if  (NuWroTTree->out[ii].pdg==2212) p++;
	      ii++;
	  }
	if      (NuWroTTree->flag.cc &&  pip &&  p) bin = 3.;
	else if (NuWroTTree->flag.cc &&  pi0 &&  p) bin = 4.;
	else if (NuWroTTree->flag.cc &&  pip &&  n) bin = 5.;
	else if (NuWroTTree->flag.nc &&  pi0 &&  p) bin = 6.;
	else if (NuWroTTree->flag.nc &&  pip &&  n) bin = 7.;
	else if (NuWroTTree->flag.nc &&  pi0 &&  n) bin = 8.;
	else if (NuWroTTree->flag.nc &&  pim &&  p) bin = 9.;
	else      
	  std::cout << "NuWroGen: bin=0 events, cc?" << NuWroTTree->flag.cc << " nc?: " << NuWroTTree->flag.nc << " Num protons: " << p << " Num pi+s: " << pip << " Num pi-s: " << pim << " Num pi0s: " << pi0 << " Num neutrons: " << n  << std::endl;
      }

    unsigned int ii(0);
    unsigned int p(0), n(0), pip(0), pim(0), pi0(0), pThresh(0.);
    while(ii<NuWroTTree->post.size())
      {
	if  (NuWroTTree->out[ii].pdg==211) pip++;
	if  (NuWroTTree->out[ii].pdg==-211) pim++;
	if  (NuWroTTree->out[ii].pdg==111) pi0++;
	if  (NuWroTTree->out[ii].pdg==2112) n++;
	if  (NuWroTTree->out[ii].pdg==2212) p++;
	if  (NuWroTTree->out[ii].pdg==2212 && 
	     (NuWroTTree->out[ii].t/1000.-0.939)>0.050) pThresh++;
	ii++;
      }

    if      (NuWroTTree->flag.cc && p==0 && (pip+pim)==0) binNew = 1;
    else if (NuWroTTree->flag.cc && p==1 && (pip+pim)==0) binNew = 2;
    else if (NuWroTTree->flag.cc && p>=2 && (pip+pim)==0) binNew = 3;
    else if (NuWroTTree->flag.cc && p==0 && (pip+pim)==1) binNew = 4;
    else if (NuWroTTree->flag.cc && p==1 && (pip+pim)==1) binNew = 5;
    else if (NuWroTTree->flag.cc && p>=2 && (pip+pim)==1) binNew = 6;
    else if (NuWroTTree->flag.cc && p==0 && (pip+pim)>=2) binNew = 7;
    else if (NuWroTTree->flag.cc && p>=1 && (pip+pim)>=2) binNew = 8;
    else if (NuWroTTree->flag.nc)                         binNew = 9;
    else binNew = 10;

    if      (NuWroTTree->flag.cc && pThresh==0 && (pip+pim)==0) binNewThresh = 1;
    else if (NuWroTTree->flag.cc && pThresh==1 && (pip+pim)==0) binNewThresh = 2;
    else if (NuWroTTree->flag.cc && pThresh>=2 && (pip+pim)==0) binNewThresh = 3;
    else if (NuWroTTree->flag.cc && pThresh==0 && (pip+pim)==1) binNewThresh = 4;
    else if (NuWroTTree->flag.cc && pThresh==1 && (pip+pim)==1) binNewThresh = 5;
    else if (NuWroTTree->flag.cc && pThresh>=2 && (pip+pim)==1) binNewThresh = 6;
    else if (NuWroTTree->flag.cc && pThresh==0 && (pip+pim)>=2) binNewThresh = 7;
    else if (NuWroTTree->flag.cc && pThresh>=1 && (pip+pim)>=2) binNewThresh = 8;
    else if (NuWroTTree->flag.nc                            ) binNewThresh = 9;
    else binNewThresh = 10;

    // Needs to be normalized to 70 tonnes, 6e20pot
    double N_ArAtoms(70.*1000*1000/40*6.022e23);
    // arXiv:pdf/0806.1449v2.pdf adjusted to uBooNE distance
    double FluxNorm(5.19e-10 * (540./460.)*(540./460.));
    double nucleiPerAtom((double)(NuWroTTree->par.nucleus_p+NuWroTTree->par.nucleus_n));
    double NumEvtsRunThisJob(10000.);
    // The weight coming out of NuWro is proportional to the xsection for that process. Instead use the flux-wtd avg's for each dyn mechanism as given in the ouutput.root.txt file.
    //    double wt = NuWroTTree->weight * N_ArAtoms * nucleiPerAtom * 6.e20 * FluxNorm / NumEvtsRunThisJob; 

    //    double wt = fxsecFluxWtd.at(NuWroTTree->dyn) * N_ArAtoms * nucleiPerAtom * 6.e20 * FluxNorm / NumEvtsRunThisJob; 
    double wt = fxsecTotal * N_ArAtoms * nucleiPerAtom * 6.e20 * FluxNorm / NumEvtsRunThisJob; 
    fDyn->Fill(bin-0.5,wt);
    fDynNew->Fill(binNew-0.5,wt);
    fDynNewThresh->Fill(binNewThresh-0.5,wt);
    f2DynNew->Fill(bin-0.5,binNew-0.5,wt);
    f2DynNewThresh->Fill(bin-0.5,binNewThresh-0.5,wt);

    fWeight->Fill(wt);
    fWeightNW->Fill(NuWroTTree->weight);
    return;
  }

}


namespace evgen{

  DEFINE_ART_MODULE(NuWroGen)

}

#endif // EVGEN_NUWROGEN_H
////////////////////////////////////////////////////////////////////////


