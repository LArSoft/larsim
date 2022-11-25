
#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "MakeReweight.h"
#include "dk2nu/tree/dk2nu.h"	
#include "dk2nu/tree/dkmeta.h"
#include "nutools/EventGeneratorBase/GENIE/MCTruthAndFriendsItr.h"

#include "TSystem.h"

namespace evwgh {
  class PPFXThinPionWeightCalc : public WeightCalc
  {
     public:
       PPFXThinPionWeightCalc();
       void Configure(fhicl::ParameterSet const& p, 
                   CLHEP::HepRandomEngine& engine) override;
       std::vector<std::vector<double> > GetWeight(art::Event & e) override;
     private:
       std::string fGenieModuleLabel;
    
       std::vector<std::string>  fInputLabels;
       std::string fPPFXMode;
       std::string fMode;
       std::string fHorn;
       std::string fTarget;
       int fSeed;
       int fVerbose;
       NeutrinoFluxReweight::MakeReweight* fPPFXrw;

     DECLARE_WEIGHTCALC(PPFXThinPionWeightCalc)
  };
  
  PPFXThinPionWeightCalc::PPFXThinPionWeightCalc()
  {
  }

  void PPFXThinPionWeightCalc::Configure(fhicl::ParameterSet const& p, 
                                  CLHEP::HepRandomEngine& engine)
  {
    //get configuration for this function
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    fGenieModuleLabel = p.get<std::string> ("genie_module_label");

    //ppfx setup
    fInputLabels = pset.get<std::vector<std::string>>("input_labels");
    fPPFXMode    = pset.get<std::string>("ppfx_mode");
    fVerbose     = pset.get<int>("verbose");
    fMode        = pset.get<std::string>("mode");  
    fHorn        = pset.get<std::string>("horn_curr");
    fTarget      = pset.get<std::string>("target_config");
    fSeed        = pset.get<int>("random_seed", -1);

    gSystem->Setenv("MODE", fPPFXMode.c_str());

    fPPFXrw = NeutrinoFluxReweight::MakeReweight::getInstance();
    std::cout<<"PPFX instance "<<fPPFXrw<<std::endl;
    std::string inputOptions  =std::string(getenv("PPFX_DIR"))+"/xml/inputs_"+fPPFXMode+".xml";
    if (fSeed != -1) fPPFXrw->setBaseSeed(fSeed); // Set the random seed
    std::cout << "is PPFX setup : " << fPPFXrw->AlreadyInitialized() << std::endl;  
    std::cout << "Setting PPFX, inputs: " << inputOptions << std::endl;
    std::cout << "Setting Horn Current Configuration to: " << fHorn << std::endl;
    std::cout << "Setting Target Configuration to: " << fTarget << std::endl;
    if(!(fPPFXrw->AlreadyInitialized())){
      fPPFXrw->SetOptions(inputOptions);	
    }
    std::cout << "PPFX just set with mode: " << fPPFXMode << std::endl;
  }

  std::vector<std::vector<double> > PPFXThinPionWeightCalc::GetWeight(art::Event & e)
  {
    std::vector<std::vector<double> > weight;
    evgb::MCTruthAndFriendsItr mcitr(e,fInputLabels);
    //calculate weight(s) here 

    int nmctruth=0, nmatched=0;
    bool flag = true;
    int  ievt = -1;
    std::vector<art::Handle<std::vector<bsim::Dk2Nu>>> dk2nus2;
    e.getManyByType(dk2nus2);	
    for (size_t dk2=0; dk2 < dk2nus2.size(); ++dk2) {
      art::Handle<std::vector<bsim::Dk2Nu>> dk2nus = dk2nus2[dk2];
    }
    while ( ( flag = mcitr.Next() ) ) {
      std::string label                  = mcitr.GetLabel();
      const simb::MCTruth*     pmctruth  = mcitr.GetMCTruth();
      // const simb::GTruth*  pgtruth  = mcitr.GetGTruth();
      //not-used//const simb::MCFlux*      pmcflux   = mcitr.GetMCFlux();
      const bsim::Dk2Nu*       pdk2nu    = mcitr.GetDk2Nu();
      //not-used//const bsim::NuChoice*    pnuchoice = mcitr.GetNuChoice();
      // art::Ptr<simb::MCTruth>  mctruthptr = mcitr.GetMCTruthPtr();

      ++ievt;
      ++nmctruth;
      if ( fVerbose > 0 ) {
        std::cout << "FluxWeightCalculator [" << std::setw(4) << ievt << "] "
                  << " label \"" << label << "\" MCTruth@ " << pmctruth 
                  << " Dk2Nu@ " << pdk2nu << std::endl;
      }
      if ( ! pdk2nu ) continue;
      ++nmatched;

      //RWH//bsim::Dk2Nu*  tmp_dk2nu  = fTmpDK2NUConversorAlg.construct_toy_dk2nu( &mctruth,&mcflux);
      //RWH//bsim::DkMeta* tmp_dkmeta = fTmpDK2NUConversorAlg.construct_toy_dkmeta(&mctruth,&mcflux);
      //RWH// those appear to have been memory leaks
      //RWH// herein begins the replacment for the "construct_toy_dkmeta"

      static bsim::DkMeta dkmeta_obj;        //RWH// create this on stack (destroyed when out-of-scope)  ... or static
      dkmeta_obj.tgtcfg  = fTarget;
      dkmeta_obj.horncfg = fHorn; 
      (dkmeta_obj.vintnames).push_back("Index_Tar_In_Ancestry");
      (dkmeta_obj.vintnames).push_back("Playlist");
      bsim::DkMeta* tmp_dkmeta = &dkmeta_obj;
      //RWH// herein ends this block

      // sigh ....
      //RWH// this is the signature in NeutrinoFluxReweight::MakeReweight :
      //      //! create an interaction chain from the new dk2nu(dkmeta) format
      //      void calculateWeights(bsim::Dk2Nu* nu, bsim::DkMeta* meta);
      //RWH// these _should_ be "const <class>*" because we don't need to change them
      //RWH// and the pointers we get out of the ART record are going to be const.
      bsim::Dk2Nu* tmp_dk2nu = const_cast<bsim::Dk2Nu*>(pdk2nu);  // remove const-ness
      try {
	fPPFXrw->calculateWeights(tmp_dk2nu,tmp_dkmeta);
      } catch (...) {
	std::cerr<<"Failed to calcualte weight"<<std::endl;
	continue;	
      }
      //Get weights:
      if (fMode=="reweight") {
	double ppfx_cv_wgt = fPPFXrw->GetCVWeight();
	std::vector<double> wvec = {ppfx_cv_wgt};
	weight.push_back(wvec);
      } else {
	std::vector<double> vttpcpion    = fPPFXrw->GetWeights("ThinTargetpCPion");
	
	std::vector<double> tmp_vhptot;
	for(unsigned int iuniv=0;iuniv<vttpcpion.size();iuniv++){
	  tmp_vhptot.push_back(float(vttpcpion[iuniv]));
	}
	weight.push_back(tmp_vhptot);
      }
    }
    return weight;
  }
  REGISTER_WEIGHTCALC(PPFXThinPionWeightCalc)
}
