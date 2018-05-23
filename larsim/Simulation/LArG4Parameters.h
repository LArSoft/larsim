////////////////////////////////////////////////////////////////////////
/// \file  LArG4Parameters.h
/// \brief Store parameters for running LArG4
///
/// \author bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
//
// This service exists to pass parameters to various different
// classes in LArG4, which are not necessary directly called by
// the LArG4_module class.
//
// Ben Jones, MIT, March 2010


#include <string>

#ifndef LArG4Parameters_h
#define LArG4Parameters_h 1

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"

namespace sim {

  class LArG4Parameters {
  public:
    LArG4Parameters(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~LArG4Parameters() {}
    
    void reconfigure(fhicl::ParameterSet const& pset);
    
    int    OpVerbosity()                                      const { return fOpVerbosity;            }
    double ParticleKineticEnergyCut()                         const { return fParticleKineticECut;    }
    bool   StoreTrajectories()                                const { return fStoreTrajectories;      }
    bool   DrawNeutrals()                                     const { return fDrawNeutrals;           }
    double VisualizationEnergyCut()                           const { return fVisualizationEnergyCut; }
    bool   UseCustomPhysics()                                 const { return fUseCustomPhysics;       }
    double RecombA()                                          const { return util::kRecombA;          }
    double Recombk()                                          const { return util::kRecombk;          }
    double ModBoxA()                                          const { return util::kModBoxA;          }
    double ModBoxB()                                          const { return util::kModBoxB;          }
    bool   UseModBoxRecomb()                                  const { return fUseModBoxRecomb;        }
    double GeVToElectrons()                                   const { return util::kGeVToElectrons;   }
    double LongitudinalDiffusion()                            const { return fLongitudinalDiffusion;  }
    double TransverseDiffusion()                              const { return fTransverseDiffusion;    }
    double ElectronClusterSize()                              const { return fElectronClusterSize;    }
    int    MinNumberOfElCluster()			      const { return fMinNumberOfElCluster;   }
    const std::vector<std::string>& EnabledPhysics()          const { return fEnabledPhysics;         }
    int    K0Bias()                                           const { return fK0Bias;                 }
    int    MNXBias()                                          const { return fXBias;                  }
    int    MNXSBias()                                         const { return fXSBias;                 }
    bool   KeepEMShowerDaughters()                            const { return fKeepEMShowerDaughters;  }
    bool   DisableWireplanes()                                const { return fDisableWireplanes;      }
    const std::vector<unsigned short int> SkipWireSignalInTPCs() const { return fSkipWireSignalInTPCs;}
    const std::string IonAndScintCalculator()                 const { return fIonAndScintCalculator;  }
    const std::vector<std::string> OpticalParamVolumes()      const { return fOpticalParamVolumes;    }
    const std::vector<std::string> OpticalParamModels()       const { return fOpticalParamModels;     }
    const std::vector<int>         OpticalParamOrientations() const { return fOpticalParamOrientations;}
    const std::vector<std::vector<std::vector<double>>> OpticalParamParameters() const{return fOpticalParamParameters;  }
    bool UseLitePhotons()                                     const { return fLitePhotons;            }

    bool   FillSimEnergyDeposits()                            const { return fFillSimEnergyDeposits;  }
    bool   NoElectronPropagation()                            const { return fNoElectronPropagation;  }
    bool   NoPhotonPropagation()                              const { return fNoPhotonPropagation;  }

  private:
    int                      fOpVerbosity;           ///< Verbosity of optical simulation - soon to be depricated
    double                   fParticleKineticECut;   ///< Minimum energy a particle needs before asking Geant4 
                                                     ///< to track it, GeV
    bool                     fStoreTrajectories;     ///< Whether to store full trajectories for every particle 
                                                     ///< simulated by Geant4
    bool                     fDrawNeutrals;          ///< depricated
    double                   fVisualizationEnergyCut;///< depricated, GeV
    bool                     fUseCustomPhysics;      ///< Whether to use a custom list of physics processes 
                                                     ///< or the default
    double                   fLongitudinalDiffusion; ///< Amount of diffusion in the longitudinal direction, cm^2/ns
    double                   fTransverseDiffusion;   ///< Amount of diffusion in the transverse direction, cm^2/ns
    double                   fElectronClusterSize;   ///< Number of ionization electrons in a given cluster 
                                                     ///< to be simulated in the readout simulation
    int 		     fMinNumberOfElCluster;   ///< Minimum number of electron clusters
    std::vector<std::string> fEnabledPhysics;        ///< List of enabled physics processes if using Custom physics
    int                      fK0Bias;                ///< Turns on secondary particle bias for K0, Lambda, 
                                                     ///< neutrons in MuNuclear
    int                      fXSBias;                ///< Turns on cross-section bian in MuNuclear
    int                      fXBias;                 ///< Enhancement factor for cross-section bian in MuNuclear, 
                                                     ///< should be <= 100
    bool                     fKeepEMShowerDaughters; ///< Whether to keep the secondary, tertiary, etc. 
                                                     ///< particles from an EM shower in the output
    bool                     fDisableWireplanes;     ///< Turn of LAr sensitivity and remove charge 
                                                     ///< drift simulation - use for running pure optical sims 
    std::vector<unsigned short int> fSkipWireSignalInTPCs;  ///< selective disabling of drift simulation 
    bool                     fUseModBoxRecomb;       ///< Use Modified Box model recombination instead of Birks
    std::string              fIonAndScintCalculator; ///< Name of algorithm to use to calculate the number of 
                                                     ///< ionization electrons and scintillation photons
                                                     ///< for each G4 step, used by 
                                                     ///< LArG4/IonizationAndScintillation.cxx
    std::vector<std::string> fOpticalParamVolumes;   ///< List of volume names which have parameterized 
                                                     ///< optical models
    std::vector<std::string> fOpticalParamModels;    ///< List of names of those models
    std::vector<int>         fOpticalParamOrientations; ///< List of orientations of (eg wireplane) in each 
                                                        ///< param volume
    std::vector<std::vector<std::vector<double> > > fOpticalParamParameters; ///< Model dependent list of 
                                                                             ///< parameters for optically 
                                                                             ///< paramaterized volumes
 
    bool fLitePhotons;

    bool   fFillSimEnergyDeposits;          ///< handle to fill SimEdeps or not
    //size_t fInitialSimEnergyDepositSize;    ///< reserve size for the edep collection in LArG4
    bool   fNoElectronPropagation;          ///< specifically prevents electron propagation
    bool   fNoPhotonPropagation;          ///< specifically prevents photon propagation in opfast
  };
}


DECLARE_ART_SERVICE(sim::LArG4Parameters, LEGACY)
#endif
