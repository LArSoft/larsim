////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parmeters.cxx
/// \brief Store parameters for running LArG4
///
/// \version $Id: LArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// This service exists to pass parameters to various different
// classes in LArG4, which are not necessary directly called by
// the LArG4_module class.
//
// Ben Jones, MIT, March 2010


// Framework includes

#include "Simulation/LArG4Parameters.h"

namespace sim {

  //--------------------------------------------------------------------------
  LArG4Parameters::LArG4Parameters(fhicl::ParameterSet const& pset, art::ActivityRegistry& /* reg */)
  {
    this->reconfigure(pset);
  }

  //--------------------------------------------------------------------------
  void LArG4Parameters::reconfigure(fhicl::ParameterSet const& pset)
  {

    fOpVerbosity            = pset.get< int                      >("OpticalSimVerbosity"     );     
    fParticleKineticECut    = pset.get< double 			 >("ParticleKineticEnergyCut");
    fStoreTrajectories      = pset.get< bool   			 >("StoreTrajectories"       );	     
    fDrawNeutrals           = pset.get< bool   			 >("VisualizeNeutrals"       );	     
    fVisualizationEnergyCut = pset.get< double 			 >("VisualizationEnergyCut"  );  
    fUseCustomPhysics       = pset.get< bool   			 >("UseCustomPhysics"        );	     
    fKeepEMShowerDaughters  = pset.get< bool                     >("KeepEMShowerDaughters"   );
    fLongitudinalDiffusion  = pset.get< double 			 >("LongitudinalDiffusion"   );   
    fTransverseDiffusion    = pset.get< double 			 >("TransverseDiffusion"     );     
    fElectronClusterSize    = pset.get< double 			 >("ElectronClusterSize"     );     
    fEnabledPhysics         = pset.get< std::vector<std::string> >("EnabledPhysics"          );
    fK0Bias                 = pset.get< int                      >("CosmogenicK0Bias"        );	     
    fXBias                  = pset.get< int    			 >("CosmogenicXSMNBiasOn"    );    
    fXSBias                 = pset.get< int     		 >("CosmogenicXSMNBiasFactor");
    // First of last 3 flags above turns on secondary particle bias for 
    // K0s,Lambdas,neutrons in MuNuclear. 
    // The second turns on cross-section bias in MuNuclear.
    // The 3rd is the enhancement factor for XS bias in MuNuclear. Keep it 
    // <=100.

    fIonAndScintCalculator   = pset.get< std::string              >("IonAndScintCalculator", "Separate");
    fDisableWireplanes       = pset.get< bool                     >("DisableWireplanes"       );
    fUseModBoxRecomb         = pset.get< bool                     >("UseModBoxRecomb"         );
    fOpticalParamVolumes     = pset.get< std::vector<std::string> >("OpticalParamVolumes"     );
    fOpticalParamModels      = pset.get< std::vector<std::string> >("OpticalParamModels"      );
    fOpticalParamOrientations= pset.get< std::vector<int>         >("OpticalParamOrientations");
    fOpticalParamParameters  = pset.get< std::vector<std::vector<std::vector<double> > > >("OpticalParamParameters");
    fLitePhotons             = pset.get< bool                     >("UseLitePhotons"       );
    fEnableSCE               = pset.get< bool                     >("EnableSCE"            );

    return;
  }

}


namespace sim {

  DEFINE_ART_SERVICE(LArG4Parameters)

}
