// GenieWeightCalc.cxx
//
// Handles event weights for GENIE systematics studies
//
// Updated by Marco Del Tutto on Feb 18 2017
// Ported from uboonecode to larsim on Feb 14 2017 
//   by Marco Del Tutto <marco.deltutto@physics.ox.ac.uk>

#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nutools/NuReweight/art/NuReweight.h" //GENIEReweight.h"
#include "nutools/NuReweight/GENIEReweight.h"
#include "nutools/NuReweight/ReweightLabels.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

namespace evwgh {
  class GenieWeightCalc : public WeightCalc
  {
  public:
    GenieWeightCalc();
    void Configure(fhicl::ParameterSet const& pset,
                   CLHEP::HepRandomEngine& engine) override;
    std::vector<std::vector<double> > GetWeight(art::Event & e) override;
    
  private:
    // The reweighting utility class:
    std::vector<rwgt::NuReweight> reweightVector;

    std::string fGenieModuleLabel;

    // What follows is the list of sereighting parameters present in LArSoft.
    // Parameters with a (*) contains more that one reweighing parameter at the same time. 
    // They can only be modified changing the relative method in $NUTOOLS_DIR/source/NuReweight/GENIEReweight.cxx
      
    enum EReweight {kNCELaxial,        // Axial mass for NC elastic
                    kNCELeta,          // Strange axial form factor for NC elastic
                    kQEMA,             // Axial mass for CC quasi-elastic
                    kQEVec,            // Choice of CCQE vector form factor (sigma = 0 => BBA05; sigma = 1 => Dipole)
                    kCCResAxial,       // Axial mass for CC resonance neutrino production
                    kCCResVector,      // Vector mass for CC resonance neutrino production
                    kResGanged,        // CC Res && NC Res (NOT ACTIVE)
                    kNCResAxial,       // Axial mass for NC resonance neutrino production
                    kNCResVector,      // Vector mass for NC resonance neutrino production
		    kCohMA,            // Axial mass for CC and NC coherent pion production
                    kCohR0,            // Nuclear size param. controlling pi absorption in Rein-Sehgal model
                    kNonResRvp1pi,     // v+p and vbar + n (1 pi) type interactions (*)
                    kNonResRvbarp1pi,  // v+n and vbar + p (1 pi) type interactions (*)
                    kNonResRvp2pi,     // v+p and vbar + n (2 pi) type interactions (*)
		    kNonResRvbarp2pi,  // v+n and vbar + p (2 pi) type interactions (*)
                    kResDecayGamma,    // BR for radiative resonance decay
                    kResDecayEta,      // BR for single-eta resonance decay
                    kResDecayTheta,    // Pion angular distibution in Delta -> pi N (sigma = 0 => isotropic; sigma = 1 => RS)
                    kNC,
                    kDISAth,           // Ath higher twist param in BY model scaling variable xi_w
                    kDISBth,           // Bth higher twist param in BY model scaling variable xi_w
                    kDISCv1u,          // Cv1u u valence GRV98 PDF correction param in BY model
                    kDISCv2u,          // Cv2u u valence GRV98 PDF correction param in BY model
                    kDISnucl,          // NOT IMPLEMENTED IN GENIE
		    kAGKYxF,           // Pion Feynman x for Npi states in AGKY
                    kAGKYpT,           // Pion transverse momentum for Npi states in AGKY
                    kFormZone,         // Hadron Formation Zone
                    kFermiGasModelKf,  // CCQE Pauli Suppression via changes in Fermi level kF
                    kFermiGasModelSf,  // Choice of model (sigma = 0 => FermiGas; sigma = 1 => SF (spectral function))
                    kIntraNukeNmfp,    // Nucleon mean free path (total rescattering probability)
                    kIntraNukeNcex,    // Nucleon charge exchange probability
                    kIntraNukeNel,     // Nucleon elastic reaction probability
                    kIntraNukeNinel,   // Nucleon inelastic reaction probability
                    kIntraNukeNabs,    // Nucleon absorption probability
                    kIntraNukeNpi,     // Nucleon pi-production probability
                    kIntraNukePImfp,   // Pi mean free path (total rescattering probability)
                    kIntraNukePIcex,   // Pi charge exchange probability
                    kIntraNukePIel,    // Pi elastic reaction probability
                    kIntraNukePIinel,  // Pi inelastic reaction probability
                    kIntraNukePIabs,   // Pi absorption probability
                    kIntraNukePIpi,    // Pi pi-production probability
                    kNReWeights};      // ?
    
    DECLARE_WEIGHTCALC(GenieWeightCalc)
  };

  GenieWeightCalc::GenieWeightCalc()
  {}

  void GenieWeightCalc::Configure(fhicl::ParameterSet const& p,
                                  CLHEP::HepRandomEngine& engine)
  {
    // Global Config
    fGenieModuleLabel = p.get<std::string> ("genie_module_label");
    fhicl::ParameterSet const &pset = p.get<fhicl::ParameterSet> (GetName());

    // Calc Config
    auto const pars = pset.get<std::vector<std::string>>("parameter_list");
    auto const parsigmas = pset.get<std::vector<float>>("parameter_sigma");
    auto const mode = pset.get<std::string>("mode");

    if (pars.size() != parsigmas.size())
      throw cet::exception(__PRETTY_FUNCTION__) << GetName()
            << "::Bad fcl configuration. parameter_list and parameter_sigma need to have same number of parameters." << std::endl;

    auto number_of_multisims = pset.get<int>("number_of_multisims");
      
    std::vector<EReweight> erwgh;
    for (auto const& s : pars) {
      if      (s == "NCELaxial") erwgh.push_back(kNCELaxial);
      else if (s == "NCELeta") erwgh.push_back(kNCELeta);
      else if (s == "QEMA") erwgh.push_back(kQEMA);
      else if (s == "QEVec") erwgh.push_back(kQEVec);
      else if (s == "CCResAxial") erwgh.push_back(kCCResAxial);
      else if (s == "CCResVector") erwgh.push_back(kCCResVector);
      else if (s == "ResGanged") erwgh.push_back(kResGanged);
      else if (s == "NCResAxial") erwgh.push_back(kNCResAxial);
      else if (s == "NCResVector") erwgh.push_back(kNCResVector);
      else if (s == "CohMA") erwgh.push_back(kCohMA);
      else if (s == "CohR0") erwgh.push_back(kCohR0);
      else if (s == "NonResRvp1pi") erwgh.push_back(kNonResRvp1pi);
      else if (s == "NonResRvbarp1pi") erwgh.push_back(kNonResRvbarp1pi);
      else if (s == "NonResRvp2pi") erwgh.push_back(kNonResRvp2pi);
      else if (s == "NonResRvbarp2pi") erwgh.push_back(kNonResRvbarp2pi);
      else if (s == "ResDecayGamma") erwgh.push_back(kResDecayGamma);
      else if (s == "ResDecayEta") erwgh.push_back(kResDecayEta);
      else if (s == "ResDecayTheta") erwgh.push_back(kResDecayTheta);
      else if (s == "NC") erwgh.push_back(kNC);
      else if (s == "DISAth") erwgh.push_back(kDISAth);
      else if (s == "DISBth") erwgh.push_back(kDISBth);
      else if (s == "DISCv1u") erwgh.push_back(kDISCv1u);
      else if (s == "DISCv2u") erwgh.push_back(kDISCv2u);
      else if (s == "DISnucl") erwgh.push_back(kDISnucl);
      else if (s == "AGKYxF") erwgh.push_back(kAGKYxF);
      else if (s == "AGKYpT") erwgh.push_back(kAGKYpT);
      else if (s == "FormZone") erwgh.push_back(kFormZone);
      else if (s == "FermiGasModelKf") erwgh.push_back(kFermiGasModelKf);
      else if (s == "FermiGasModelSf") erwgh.push_back(kFermiGasModelSf);
      else if (s == "IntraNukeNmfp") erwgh.push_back(kIntraNukeNmfp);
      else if (s == "IntraNukeNcex") erwgh.push_back(kIntraNukeNcex);
      else if (s == "IntraNukeNel") erwgh.push_back(kIntraNukeNel);
      else if (s == "IntraNukeNinel") erwgh.push_back(kIntraNukeNinel);
      else if (s == "IntraNukeNabs") erwgh.push_back(kIntraNukeNabs);
      else if (s == "IntraNukeNpi") erwgh.push_back(kIntraNukeNpi);
      else if (s == "IntraNukePImfp") erwgh.push_back(kIntraNukePImfp);
      else if (s == "IntraNukePIcex") erwgh.push_back(kIntraNukePIcex);
      else if (s == "IntraNukePIel") erwgh.push_back(kIntraNukePIel);
      else if (s == "IntraNukePIinel") erwgh.push_back(kIntraNukePIinel);
      else if (s == "IntraNukePIabs") erwgh.push_back(kIntraNukePIabs);
      else if (s == "IntraNukePIpi") erwgh.push_back(kIntraNukePIpi);
      else {
	throw cet::exception(__PRETTY_FUNCTION__) << GetName()
              << "::Physical process " << s << " you requested is not available to reweight." << std::endl;
      }
    }
      
    //Prepare sigmas
    std::vector<std::vector<float>> reweightingSigmas(erwgh.size());

    if (mode.find("pm1sigma") != std::string::npos ) { 
      number_of_multisims = 2; // only +-1 sigma if pm1sigma is specified
    }
    for (unsigned int i = 0; i < reweightingSigmas.size(); ++i) {
      reweightingSigmas[i].resize(number_of_multisims);
      for (int j = 0; j < number_of_multisims; j ++) {
	if (mode.find("multisim") != std::string::npos )
          reweightingSigmas[i][j] = parsigmas[i]*CLHEP::RandGaussQ::shoot(&engine, 0., 1.);
        else if (mode.find("pm1sigma") != std::string::npos )
          reweightingSigmas[i][j] = (j == 0 ? 1.: -1.); // j==0 => 1; j==1 => -1 if pm1sigma is specified
	else
	  reweightingSigmas[i][j] = parsigmas[i];
      }
    }

    reweightVector.resize(number_of_multisims);
    
    for (unsigned weight_point = 0, e = reweightVector.size(); weight_point != e; ++weight_point) {
      auto& driver = reweightVector[weight_point];
      for (unsigned int i_reweightingKnob = 0; i_reweightingKnob<erwgh.size(); ++i_reweightingKnob) {
        std::cout << GetName()
                  << "::Setting up rwgh " << weight_point << "\t" << i_reweightingKnob << "\t" << erwgh[i_reweightingKnob] << std::endl;

	switch (erwgh[i_reweightingKnob]){

	case kNCELaxial:
          driver.ReweightNCEL(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kNCELeta:
          driver.ReweightNCEL(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
            
	case kQEMA:
          driver.ReweightQEMA(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kQEVec:
          driver.ReweightQEVec(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kResGanged:
	  break;
        
	case kCCResAxial:
          driver.ReweightCCRes(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kCCResVector:
          driver.ReweightCCRes(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;

	case kNCResAxial:
          driver.ReweightNCRes(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kNCResVector:
          driver.ReweightNCRes(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kCohMA:
          driver.ReweightCoh(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kCohR0:
          driver.ReweightCoh(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
        
	case kNonResRvp1pi:
          driver.ReweightNonResRvp1pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvbarp1pi:
          driver.ReweightNonResRvbarp1pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvp2pi:
          driver.ReweightNonResRvp2pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvbarp2pi:
          driver.ReweightNonResRvbarp2pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kResDecayGamma:
          driver.ReweightResDecay(reweightingSigmas[i_reweightingKnob][weight_point], 0., 0.);
	  break;
        case kResDecayEta:
          driver.ReweightResDecay(0., reweightingSigmas[i_reweightingKnob][weight_point], 0.);
          break;
        case kResDecayTheta:
          driver.ReweightResDecay(0., 0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kNC:
          driver.ReweightNC(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kDISAth:
          driver.ReweightDIS(reweightingSigmas[i_reweightingKnob][weight_point], 0., 0., 0.);
	  break;
        case kDISBth:
          driver.ReweightDIS(0., reweightingSigmas[i_reweightingKnob][weight_point], 0., 0.);
          break;
        case kDISCv1u:
          driver.ReweightDIS(0., 0., reweightingSigmas[i_reweightingKnob][weight_point], 0.);
          break;
        case kDISCv2u:
          driver.ReweightDIS(0., 0., 0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kDISnucl:
          driver.ReweightDISnucl(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
        
	case kAGKYxF:
          driver.ReweightAGKY(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
	  break;
        case kAGKYpT:
          driver.ReweightAGKY(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
      
        case kFormZone:
          driver.ReweightFormZone(reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
        case kFermiGasModelKf:
          driver.ReweightFGM(reweightingSigmas[i_reweightingKnob][weight_point], 0.);
          break;
        case kFermiGasModelSf:
          driver.ReweightFGM(0., reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
        case kIntraNukeNmfp:
          driver.ReweightIntraNuke(rwgt::fReweightMFP_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNcex:
          driver.ReweightIntraNuke(rwgt::fReweightFrCEx_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNel:
          driver.ReweightIntraNuke(rwgt::fReweightFrElas_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNinel:
          driver.ReweightIntraNuke(rwgt::fReweightFrInel_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNabs:
          driver.ReweightIntraNuke(rwgt::fReweightFrAbs_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukeNpi:
          driver.ReweightIntraNuke(rwgt::fReweightFrPiProd_N, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePImfp:
          driver.ReweightIntraNuke(rwgt::fReweightMFP_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIcex:
          driver.ReweightIntraNuke(rwgt::fReweightFrCEx_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIel:
          driver.ReweightIntraNuke(rwgt::fReweightFrElas_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIinel:
          driver.ReweightIntraNuke(rwgt::fReweightFrInel_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIabs:
          driver.ReweightIntraNuke(rwgt::fReweightFrAbs_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        case kIntraNukePIpi:
          driver.ReweightIntraNuke(rwgt::fReweightFrPiProd_pi, reweightingSigmas[i_reweightingKnob][weight_point]);
          break;
        
	case kNReWeights:
	  break;
	}
      }
    } // loop over nWeights

    // Tell all of the reweight drivers to configure themselves:
    std::cout<< GetName()<<"::Setting up "<<reweightVector.size()<<" reweightcalcs"<<std::endl;
    for(auto & driver : reweightVector){
      driver.Configure();
    }

  }

  std::vector<std::vector<double> > GenieWeightCalc::GetWeight(art::Event & e)
  { 
    // Returns a vector of weights for each neutrino interaction in the event

    // Get the MC generator information out of the event       
    // These are all handles to mc information.
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;  
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
    art::Handle< std::vector<simb::GTruth> > gTruthHandle;

    // Actually go and get the stuff
    e.getByLabel(fGenieModuleLabel,mcTruthHandle);
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);
    e.getByLabel(fGenieModuleLabel,gTruthHandle);

    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mcTruthHandle);

    std::vector<art::Ptr<simb::GTruth > > glist;
    art::fill_ptr_vector(glist, gTruthHandle);

    // Calculate weight(s) here 
    std::vector<std::vector<double> >weight(mclist.size());
    for ( unsigned int inu=0; inu<mclist.size();inu++) {
      weight[inu].resize(reweightVector.size());    
      for (unsigned int i_weight = 0; 
	   i_weight < reweightVector.size(); 
	   i_weight ++){
        weight[inu][i_weight]= reweightVector[i_weight].CalcWeight(*mclist[inu],*glist[inu]);
      }
    }
    return weight;

  }
  REGISTER_WEIGHTCALC(GenieWeightCalc)
}
