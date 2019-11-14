//=============================================================================
// NueAr40CCGenerator.h
//
// Gleb Sinev, Duke, 2015
// Supernova neutrino generator that simulates nu_e-Ar40 CC interactions
// and produces a text file with kinematics information,
// which can be used as input for LArSoft
// Based on a generator written by AJ Roeth
//=============================================================================

// C++ includes
#include <vector>
#include <map>
#include <string>

namespace CLHEP { class HepRandomEngine; }

// Framework includes
namespace fhicl { class ParameterSet; }

// nusimdata includes
namespace simb { class MCTruth; }

namespace evgen {

  class NueAr40CCGenerator {

    public:

      // Constructor
      NueAr40CCGenerator(fhicl::ParameterSet const& parameterSet);

      // Simulate interactions, produce plots
      std::vector<simb::MCTruth> Generate(CLHEP::HepRandomEngine& engine);

    private:

      // Generate a direction vector isotropically
      std::vector< double > GetIsotropicDirection
                                    (CLHEP::HepRandomEngine& engine) const;

      // Generate a position distributed uniformly inside the active volume
      std::vector< double > GetUniformPosition
                                    (CLHEP::HepRandomEngine& engine) const;

      // Decide how many neutrinos to generate
      int GetNumberOfNeutrinos(CLHEP::HepRandomEngine& engine) const;

      // Generate a neutrino time distributed uniformly
      // from fNeutrinoTimeBegin to fNeutrinoTimeEnd
      double GetNeutrinoTime(CLHEP::HepRandomEngine& engine) const;

      // Use fEnergyProbabilityMap to sample one energy value
      // or return a constant value
      double GetNeutrinoEnergy(CLHEP::HepRandomEngine& engine) const;

      // Read a ROOT file with a TGraph  and fill fEnergyProbabilityMap
      // Assume that the first point in TGraph is { 0, 0 }
      void ReadNeutrinoSpectrum();

      // Initialize vectors with branching ratios, probabilities,
      // and other quantities required to calculate
      // number and energy of gammas produced
      void InitializeVectors();

      // Simulate particles
      void CreateKinematicsVector(simb::MCTruth& truth,
                                  CLHEP::HepRandomEngine& engine) const;

      bool ProcessOneNeutrino
                 (simb::MCTruth& truth, double neutrinoEnergy,
                  double neutrinoTime, CLHEP::HepRandomEngine& engine) const;
      std::vector< double > CalculateCrossSections
                            (double neutrinoEnergy, int& highestLevel) const;

      // Assume that the first entry is { 0, 0 }
      std::map< double, double > fEnergyProbabilityMap;

      int fNumberOfLevels;
      int fNumberOfStartLevels;

      // Vectors containing information about levels
      std::vector< std::vector< double > > fBranchingRatios;
      std::vector< std::vector< int    > > fDecayTo;
      std::vector< double > fStartEnergyLevels;
      std::vector< double > fB;
      std::vector< double > fEnergyLevels;

      // Generate monoenergetic neutrinos if this variable is set to true
      bool fMonoenergeticNeutrinos;
      // Energy of monoenergetic neutrinos
      double fNeutrinoEnergy;

      // Otherwise sample the spectrum
      std::string fEnergySpectrumFileName;

      // Number of neutrinos is distributed according
      // to Poisson distribution if this is true
      bool fUsePoissonDistribution;

      // Allow zero neutrinos to be created if that is what the randonmly
      // generated number-of-nu ends up being.
      bool fAllowZeroNeutrinos;

      // Average or exact number of neutrinos generated
      int fNumberOfNeutrinos;

      // Generate neutrinos uniformly in a time interval
      double fNeutrinoTimeBegin;
      double fNeutrinoTimeEnd;

      // Two opposite vertices { { x0, y0, z0 }, { x1, y1, z1 } }
      // of the active volume box, such that x0 < x1, y0 < y1, z0 < z1
      std::vector< std::vector < double > > fActiveVolume;

  };

}
