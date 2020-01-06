//=============================================================================
// NueAr40CCGenerator.cxx
//
// Gleb Sinev, Duke, 2015
//=============================================================================

#include "NueAr40CCGenerator.h"

// Framework includes
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

// nusimdata includes
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLorentzVector.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

// C++ includes
#include <cmath>

namespace evgen {

  //----------------------------------------------------------------------------
  // Constructor
  NueAr40CCGenerator::NueAr40CCGenerator(fhicl::ParameterSet
                                                           const& parameterSet)
              : fNumberOfLevels        ( 73                )
              , fNumberOfStartLevels   ( 21                )
              , fBranchingRatios       (fNumberOfLevels    )
              , fDecayTo               (fNumberOfLevels    )
              , fMonoenergeticNeutrinos
                         (parameterSet.get< bool   >("MonoenergeticNeutrinos"))
              , fNeutrinoEnergy
                         (parameterSet.get< double >("NeutrinoEnergy"        ))
              , fEnergySpectrumFileName
                    (parameterSet.get< std::string >("EnergySpectrumFileName"))
              , fUsePoissonDistribution
                         (parameterSet.get< bool   >("UsePoissonDistribution"))
              , fAllowZeroNeutrinos
                         (parameterSet.get< bool   >("AllowZeroNeutrinos"    ))
              , fNumberOfNeutrinos
                         (parameterSet.get< int    >("NumberOfNeutrinos"     ))
              , fNeutrinoTimeBegin
                         (parameterSet.get< double >("NeutrinoTimeBegin"     ))
              , fNeutrinoTimeEnd
                         (parameterSet.get< double >("NeutrinoTimeEnd"       ))
  {

    fActiveVolume.push_back
      (parameterSet.get< std::vector< double > >("ActiveVolume0"));
    fActiveVolume.push_back
      (parameterSet.get< std::vector< double > >("ActiveVolume1"));

    if (!fMonoenergeticNeutrinos) ReadNeutrinoSpectrum();

    InitializeVectors();

  }

  //----------------------------------------------------------------------------
  // Main routine
  std::vector<simb::MCTruth> NueAr40CCGenerator::Generate(CLHEP::HepRandomEngine& engine)
  {

    std::vector<simb::MCTruth> truths;

    int NumberOfNu = this->GetNumberOfNeutrinos(engine);
    for(int i = 0; i < NumberOfNu; ++i) {

      simb::MCTruth truth;
      truth.SetOrigin(simb::kSuperNovaNeutrino);

      // Loop until at least one neutrino is simulated
      while (!truth.NParticles()) {
	CreateKinematicsVector(truth, engine);
      }

      truths.push_back(truth);

    }

    return truths;

  }

  //----------------------------------------------------------------------------
  // Return normalized direction cosines of isotropic vector
  std::vector< double > NueAr40CCGenerator::GetIsotropicDirection
                                         (CLHEP::HepRandomEngine& engine) const
  {

    CLHEP::RandFlat randFlat(engine);

    std::vector< double > isotropicDirection;

    double phi      = 2*TMath::Pi()*randFlat.fire();
    double cosTheta = 2*randFlat.fire() - 1;
    double theta    = TMath::ACos(cosTheta);

    // x, y, z
    isotropicDirection.push_back(cos(phi)*sin(theta));
    isotropicDirection.push_back(sin(phi)*sin(theta));
    isotropicDirection.push_back(cosTheta);

    return isotropicDirection;

  }

  //----------------------------------------------------------------------------
  // Return a random vector with 3D coordinates inside the active LAr volume
  std::vector< double > NueAr40CCGenerator::GetUniformPosition
                                        (CLHEP::HepRandomEngine& engine) const
  {

    CLHEP::RandFlat randFlat(engine);

    std::vector< double > position;

    position.push_back(randFlat.
      fire(fActiveVolume.at(0).at(0), fActiveVolume.at(1).at(0)));
    position.push_back(randFlat.
      fire(fActiveVolume.at(0).at(1), fActiveVolume.at(1).at(1)));
    position.push_back(randFlat.
      fire(fActiveVolume.at(0).at(2), fActiveVolume.at(1).at(2)));

    return position;

  }

  //----------------------------------------------------------------------------
  // Get number of neutrinos to generate
  int NueAr40CCGenerator::GetNumberOfNeutrinos
                                         (CLHEP::HepRandomEngine& engine) const
  {

    if (fUsePoissonDistribution)
    {
      CLHEP::RandPoisson randPoisson(engine);
      int N = randPoisson.fire(fNumberOfNeutrinos);
      if(N == 0 && !fAllowZeroNeutrinos) N = 1;
      return N;
    }

    return fNumberOfNeutrinos;

  }

  //----------------------------------------------------------------------------
  // Sample uniform distribution to get a neutrino interaction time
  double NueAr40CCGenerator::GetNeutrinoTime
                                         (CLHEP::HepRandomEngine& engine) const
  {

    CLHEP::RandFlat randFlat(engine);

    return randFlat.fire(fNeutrinoTimeBegin, fNeutrinoTimeEnd);

  }

  //----------------------------------------------------------------------------
  // Sample energy spectrum from fEnergyProbabilityMap
  // or return a constant value
  double NueAr40CCGenerator::GetNeutrinoEnergy
                                         (CLHEP::HepRandomEngine& engine) const
  {

    if (fMonoenergeticNeutrinos) return fNeutrinoEnergy;

    CLHEP::RandFlat randFlat(engine);

    double neutrinoEnergy = 0.0;

    double randomNumber   = randFlat.fire();

    // We need this to get a previous entry in the map
    std::pair< double, double > previousPair;

    for (std::map< double, double >::const_iterator energyProbability =
          fEnergyProbabilityMap.begin(); energyProbability !=
          fEnergyProbabilityMap.end(); ++energyProbability)
    {
      if (randomNumber < energyProbability->second)
      {
        if (energyProbability != fEnergyProbabilityMap.begin())
        {
          neutrinoEnergy = energyProbability->first -
            (energyProbability->second - randomNumber)*
            (energyProbability->first - previousPair.first)/
            (energyProbability->second - previousPair.second);
          break;
        }
        else
        {
          neutrinoEnergy = energyProbability->first;
          break;
        }
      }
      previousPair = *energyProbability;
    }

    return neutrinoEnergy;

  }

  //----------------------------------------------------------------------------
  // Read a neutrino spectrum from a ROOT file
  void NueAr40CCGenerator::ReadNeutrinoSpectrum()
  {

    cet::search_path searchPath("FW_SEARCH_PATH");
    std::string directoryName = "SupernovaNeutrinos/" +
                                                      fEnergySpectrumFileName;

    std::string fullName;
    searchPath.find_file(directoryName, fullName);

    if (fullName.empty())
      throw cet::exception("NueAr40CCGenerator")
        << "Cannot find file with neutrino energy spectrum "
        << fullName << '\n';

    TFile energySpectrumFile(fullName.c_str(), "READ");

    std::string energySpectrumGraphName = "NueSpectrum";
    TGraph *energySpectrumGraph =
      dynamic_cast< TGraph* >(energySpectrumFile.Get
                                     (energySpectrumGraphName.c_str()));

    double  integral       = 0.0;
    int     numberOfPoints = energySpectrumGraph->GetN();
    double *energyValues   = energySpectrumGraph->GetX();
    double *fluxValues     = energySpectrumGraph->GetY();
    for (int point = 0; point < numberOfPoints; ++point)
      integral += fluxValues[point];

    double  probability    = 0.0;
    for (int point = 0; point < numberOfPoints; ++point)
    {
      probability += fluxValues[point]/integral;
      fEnergyProbabilityMap.insert(std::make_pair(energyValues[point],
                                                          probability));
    }

  }

  //----------------------------------------------------------------------------
  // Fill vectors with quantities necessary to simulate gammas and electrons
  void NueAr40CCGenerator::InitializeVectors()
  {

    // The level data -- there's probably a better way to initialize
    fBranchingRatios.at(0).push_back(1.00);
    fDecayTo.at(0).push_back(1);

    fBranchingRatios.at(1).push_back(1.00);
    fDecayTo.at(1).push_back(4);

    fBranchingRatios.at(2).push_back(0.133891);
    fBranchingRatios.at(2).push_back(0.41841);
    fBranchingRatios.at(2).push_back(0.351464);
    fBranchingRatios.at(2).push_back(0.0962343);
    fDecayTo.at(2).push_back(38);
    fDecayTo.at(2).push_back(43);
    fDecayTo.at(2).push_back(59);
    fDecayTo.at(2).push_back(72);

    fBranchingRatios.at(3).push_back(0.391705);
    fBranchingRatios.at(3).push_back(0.460829);
    fBranchingRatios.at(3).push_back(0.147465);
    fDecayTo.at(3).push_back(61);
    fDecayTo.at(3).push_back(65);
    fDecayTo.at(3).push_back(71);

    fBranchingRatios.at(4).push_back(0.358974);
    fBranchingRatios.at(4).push_back(0.641026);
    fDecayTo.at(4).push_back(13);
    fDecayTo.at(4).push_back(58);

    fBranchingRatios.at(5).push_back(0.247808);
    fBranchingRatios.at(5).push_back(0.0468929);
    fBranchingRatios.at(5).push_back(0.213496);
    fBranchingRatios.at(5).push_back(0.11056);
    fBranchingRatios.at(5).push_back(0.381243);
    fDecayTo.at(5).push_back(40);
    fDecayTo.at(5).push_back(52);
    fDecayTo.at(5).push_back(67);
    fDecayTo.at(5).push_back(71);
    fDecayTo.at(5).push_back(72);

    fBranchingRatios.at(6).push_back(0.361025);
    fBranchingRatios.at(6).push_back(0.056677);
    fBranchingRatios.at(6).push_back(0.194099);
    fBranchingRatios.at(6).push_back(0.388199);
    fDecayTo.at(6).push_back(52);
    fDecayTo.at(6).push_back(55);
    fDecayTo.at(6).push_back(63);
    fDecayTo.at(6).push_back(68);

    fBranchingRatios.at(7).push_back(0.0300963);
    fBranchingRatios.at(7).push_back(0.0613965);
    fBranchingRatios.at(7).push_back(0.0152488);
    fBranchingRatios.at(7).push_back(0.0321027);
    fBranchingRatios.at(7).push_back(0.0513644);
    fBranchingRatios.at(7).push_back(0.0682183);
    fBranchingRatios.at(7).push_back(0.073435);
    fBranchingRatios.at(7).push_back(0.0100321);
    fBranchingRatios.at(7).push_back(0.0120385);
    fBranchingRatios.at(7).push_back(0.0922953);
    fBranchingRatios.at(7).push_back(0.152488);
    fBranchingRatios.at(7).push_back(0.401284);
    fDecayTo.at(7).push_back(17);
    fDecayTo.at(7).push_back(25);
    fDecayTo.at(7).push_back(27);
    fDecayTo.at(7).push_back(32);
    fDecayTo.at(7).push_back(49);
    fDecayTo.at(7).push_back(54);
    fDecayTo.at(7).push_back(56);
    fDecayTo.at(7).push_back(62);
    fDecayTo.at(7).push_back(63);
    fDecayTo.at(7).push_back(67);
    fDecayTo.at(7).push_back(68);
    fDecayTo.at(7).push_back(70);

    fBranchingRatios.at(8).push_back(0.0089877);
    fBranchingRatios.at(8).push_back(0.06386);
    fBranchingRatios.at(8).push_back(0.13245);
    fBranchingRatios.at(8).push_back(0.473037);
    fBranchingRatios.at(8).push_back(0.321665);
    fDecayTo.at(8).push_back(43);
    fDecayTo.at(8).push_back(56);
    fDecayTo.at(8).push_back(67);
    fDecayTo.at(8).push_back(70);
    fDecayTo.at(8).push_back(71);

    fBranchingRatios.at(9).push_back(0.16129);
    fBranchingRatios.at(9).push_back(0.193548);
    fBranchingRatios.at(9).push_back(0.645161);
    fDecayTo.at(9).push_back(37);
    fDecayTo.at(9).push_back(65);
    fDecayTo.at(9).push_back(72);

    fBranchingRatios.at(10).push_back(0.245826);
    fBranchingRatios.at(10).push_back(0.0816327);
    fBranchingRatios.at(10).push_back(0.20872);
    fBranchingRatios.at(10).push_back(0.463822);
    fDecayTo.at(10).push_back(40);
    fDecayTo.at(10).push_back(56);
    fDecayTo.at(10).push_back(60);
    fDecayTo.at(10).push_back(71);

    fBranchingRatios.at(11).push_back(0.170984);
    fBranchingRatios.at(11).push_back(0.310881);
    fBranchingRatios.at(11).push_back(0.518135);
    fDecayTo.at(11).push_back(29);
    fDecayTo.at(11).push_back(54);
    fDecayTo.at(11).push_back(66);

    fBranchingRatios.at(12).push_back(0.242424);
    fBranchingRatios.at(12).push_back(0.757576);
    fDecayTo.at(12).push_back(54);
    fDecayTo.at(12).push_back(62);

    fBranchingRatios.at(13).push_back(0.159664);
    fBranchingRatios.at(13).push_back(0.840336);
    fDecayTo.at(13).push_back(48);
    fDecayTo.at(13).push_back(58);

    fBranchingRatios.at(14).push_back(0.288136);
    fBranchingRatios.at(14).push_back(0.145763);
    fBranchingRatios.at(14).push_back(0.118644);
    fBranchingRatios.at(14).push_back(0.108475);
    fBranchingRatios.at(14).push_back(0.338983);
    fDecayTo.at(14).push_back(56);
    fDecayTo.at(14).push_back(66);
    fDecayTo.at(14).push_back(70);
    fDecayTo.at(14).push_back(71);
    fDecayTo.at(14).push_back(72);

    fBranchingRatios.at(15).push_back(0.00869263);
    fBranchingRatios.at(15).push_back(0.087274);
    fBranchingRatios.at(15).push_back(0.0973574);
    fBranchingRatios.at(15).push_back(0.288595);
    fBranchingRatios.at(15).push_back(0.347705);
    fBranchingRatios.at(15).push_back(0.170376);
    fDecayTo.at(15).push_back(40);
    fDecayTo.at(15).push_back(64);
    fDecayTo.at(15).push_back(65);
    fDecayTo.at(15).push_back(68);
    fDecayTo.at(15).push_back(70);
    fDecayTo.at(15).push_back(71);

    fBranchingRatios.at(16).push_back(1);
    fDecayTo.at(16).push_back(65);

    fBranchingRatios.at(17).push_back(0.341763);
    fBranchingRatios.at(17).push_back(0.0567327);
    fBranchingRatios.at(17).push_back(0.164046);
    fBranchingRatios.at(17).push_back(0.239234);
    fBranchingRatios.at(17).push_back(0.198223);
    fDecayTo.at(17).push_back(35);
    fDecayTo.at(17).push_back(39);
    fDecayTo.at(17).push_back(51);
    fDecayTo.at(17).push_back(67);
    fDecayTo.at(17).push_back(69);

    fBranchingRatios.at(18).push_back(0.106406);
    fBranchingRatios.at(18).push_back(0.254103);
    fBranchingRatios.at(18).push_back(0.0465855);
    fBranchingRatios.at(18).push_back(0.529381);
    fBranchingRatios.at(18).push_back(0.0635257);
    fDecayTo.at(18).push_back(60);
    fDecayTo.at(18).push_back(61);
    fDecayTo.at(18).push_back(63);
    fDecayTo.at(18).push_back(70);
    fDecayTo.at(18).push_back(72);

    fBranchingRatios.at(19).push_back(0.0905797);
    fBranchingRatios.at(19).push_back(0.228261);
    fBranchingRatios.at(19).push_back(0.181159);
    fBranchingRatios.at(19).push_back(0.137681);
    fBranchingRatios.at(19).push_back(0.362319);
    fDecayTo.at(19).push_back(43);
    fDecayTo.at(19).push_back(45);
    fDecayTo.at(19).push_back(52);
    fDecayTo.at(19).push_back(70);
    fDecayTo.at(19).push_back(71);

    fBranchingRatios.at(20).push_back(0.0316414);
    fBranchingRatios.at(20).push_back(0.0415293);
    fBranchingRatios.at(20).push_back(0.0362558);
    fBranchingRatios.at(20).push_back(0.0481213);
    fBranchingRatios.at(20).push_back(0.0903098);
    fBranchingRatios.at(20).push_back(0.0929466);
    fBranchingRatios.at(20).push_back(0.659196);
    fDecayTo.at(20).push_back(30);
    fDecayTo.at(20).push_back(32);
    fDecayTo.at(20).push_back(46);
    fDecayTo.at(20).push_back(61);
    fDecayTo.at(20).push_back(64);
    fDecayTo.at(20).push_back(66);
    fDecayTo.at(20).push_back(70);

    fBranchingRatios.at(21).push_back(0.310757);
    fBranchingRatios.at(21).push_back(0.398406);
    fBranchingRatios.at(21).push_back(0.290837);
    fDecayTo.at(21).push_back(64);
    fDecayTo.at(21).push_back(66);
    fDecayTo.at(21).push_back(70);

    fBranchingRatios.at(22).push_back(0.152542);
    fBranchingRatios.at(22).push_back(0.847458);
    fDecayTo.at(22).push_back(67);
    fDecayTo.at(22).push_back(71);

    fBranchingRatios.at(23).push_back(0.168367);
    fBranchingRatios.at(23).push_back(0.321429);
    fBranchingRatios.at(23).push_back(0.510204);
    fDecayTo.at(23).push_back(52);
    fDecayTo.at(23).push_back(70);
    fDecayTo.at(23).push_back(71);

    fBranchingRatios.at(24).push_back(0.0276437);
    fBranchingRatios.at(24).push_back(0.078982);
    fBranchingRatios.at(24).push_back(0.0245722);
    fBranchingRatios.at(24).push_back(0.162352);
    fBranchingRatios.at(24).push_back(0.184291);
    fBranchingRatios.at(24).push_back(0.438789);
    fBranchingRatios.at(24).push_back(0.0833699);
    fDecayTo.at(24).push_back(36);
    fDecayTo.at(24).push_back(53);
    fDecayTo.at(24).push_back(62);
    fDecayTo.at(24).push_back(64);
    fDecayTo.at(24).push_back(70);
    fDecayTo.at(24).push_back(71);
    fDecayTo.at(24).push_back(72);

    fBranchingRatios.at(25).push_back(0.0163934);
    fBranchingRatios.at(25).push_back(0.336276);
    fBranchingRatios.at(25).push_back(0.226986);
    fBranchingRatios.at(25).push_back(0.420345);
    fDecayTo.at(25).push_back(43);
    fDecayTo.at(25).push_back(67);
    fDecayTo.at(25).push_back(68);
    fDecayTo.at(25).push_back(70);

    fBranchingRatios.at(26).push_back(0.184332);
    fBranchingRatios.at(26).push_back(0.0460829);
    fBranchingRatios.at(26).push_back(0.460829);
    fBranchingRatios.at(26).push_back(0.078341);
    fBranchingRatios.at(26).push_back(0.230415);
    fDecayTo.at(26).push_back(53);
    fDecayTo.at(26).push_back(54);
    fDecayTo.at(26).push_back(60);
    fDecayTo.at(26).push_back(61);
    fDecayTo.at(26).push_back(71);

    fBranchingRatios.at(27).push_back(0.0147145);
    fBranchingRatios.at(27).push_back(0.0170689);
    fBranchingRatios.at(27).push_back(0.0500294);
    fBranchingRatios.at(27).push_back(0.329606);
    fBranchingRatios.at(27).push_back(0.588582);
    fDecayTo.at(27).push_back(36);
    fDecayTo.at(27).push_back(46);
    fDecayTo.at(27).push_back(56);
    fDecayTo.at(27).push_back(67);
    fDecayTo.at(27).push_back(68);

    fBranchingRatios.at(28).push_back(1);
    fDecayTo.at(28).push_back(70);

    fBranchingRatios.at(29).push_back(0.280702);
    fBranchingRatios.at(29).push_back(0.0935673);
    fBranchingRatios.at(29).push_back(0.0409357);
    fBranchingRatios.at(29).push_back(0.584795);
    fDecayTo.at(29).push_back(63);
    fDecayTo.at(29).push_back(66);
    fDecayTo.at(29).push_back(68);
    fDecayTo.at(29).push_back(70);

    fBranchingRatios.at(30).push_back(0.393701);
    fBranchingRatios.at(30).push_back(0.283465);
    fBranchingRatios.at(30).push_back(0.188976);
    fBranchingRatios.at(30).push_back(0.133858);
    fDecayTo.at(30).push_back(61);
    fDecayTo.at(30).push_back(67);
    fDecayTo.at(30).push_back(71);
    fDecayTo.at(30).push_back(72);

    fBranchingRatios.at(31).push_back(0.0477737);
    fBranchingRatios.at(31).push_back(0.194805);
    fBranchingRatios.at(31).push_back(0.0245826);
    fBranchingRatios.at(31).push_back(0.269017);
    fBranchingRatios.at(31).push_back(0.463822);
    fDecayTo.at(31).push_back(45);
    fDecayTo.at(31).push_back(60);
    fDecayTo.at(31).push_back(65);
    fDecayTo.at(31).push_back(71);
    fDecayTo.at(31).push_back(72);

    fBranchingRatios.at(32).push_back(0.105769);
    fBranchingRatios.at(32).push_back(0.129808);
    fBranchingRatios.at(32).push_back(0.0528846);
    fBranchingRatios.at(32).push_back(0.480769);
    fBranchingRatios.at(32).push_back(0.230769);
    fDecayTo.at(32).push_back(46);
    fDecayTo.at(32).push_back(56);
    fDecayTo.at(32).push_back(66);
    fDecayTo.at(32).push_back(70);
    fDecayTo.at(32).push_back(71);

    fBranchingRatios.at(33).push_back(0.21875);
    fBranchingRatios.at(33).push_back(0.78125);
    fDecayTo.at(33).push_back(67);
    fDecayTo.at(33).push_back(71);

    fBranchingRatios.at(34).push_back(0.102011);
    fBranchingRatios.at(34).push_back(0.179598);
    fBranchingRatios.at(34).push_back(0.718391);
    fDecayTo.at(34).push_back(43);
    fDecayTo.at(34).push_back(61);
    fDecayTo.at(34).push_back(66);

    fBranchingRatios.at(35).push_back(0.00826446);
    fBranchingRatios.at(35).push_back(0.393546);
    fBranchingRatios.at(35).push_back(0.334514);
    fBranchingRatios.at(35).push_back(0.263676);
    fDecayTo.at(35).push_back(64);
    fDecayTo.at(35).push_back(67);
    fDecayTo.at(35).push_back(68);
    fDecayTo.at(35).push_back(70);

    fBranchingRatios.at(36).push_back(0.056338);
    fBranchingRatios.at(36).push_back(0.704225);
    fBranchingRatios.at(36).push_back(0.239437);
    fDecayTo.at(36).push_back(51);
    fDecayTo.at(36).push_back(70);
    fDecayTo.at(36).push_back(71);

    fBranchingRatios.at(37).push_back(0.21875);
    fBranchingRatios.at(37).push_back(0.78125);
    fDecayTo.at(37).push_back(67);
    fDecayTo.at(37).push_back(70);

    fBranchingRatios.at(38).push_back(0.181818);
    fBranchingRatios.at(38).push_back(0.757576);
    fBranchingRatios.at(38).push_back(0.0606061);
    fDecayTo.at(38).push_back(66);
    fDecayTo.at(38).push_back(71);
    fDecayTo.at(38).push_back(72);

    fBranchingRatios.at(39).push_back(0.157258);
    fBranchingRatios.at(39).push_back(0.403226);
    fBranchingRatios.at(39).push_back(0.237903);
    fBranchingRatios.at(39).push_back(0.201613);
    fDecayTo.at(39).push_back(62);
    fDecayTo.at(39).push_back(70);
    fDecayTo.at(39).push_back(71);
    fDecayTo.at(39).push_back(72);

    fBranchingRatios.at(40).push_back(0.0740741);
    fBranchingRatios.at(40).push_back(0.925926);
    fDecayTo.at(40).push_back(52);
    fDecayTo.at(40).push_back(72);

    fBranchingRatios.at(41).push_back(0.0535714);
    fBranchingRatios.at(41).push_back(0.35119);
    fBranchingRatios.at(41).push_back(0.595238);
    fDecayTo.at(41).push_back(67);
    fDecayTo.at(41).push_back(68);
    fDecayTo.at(41).push_back(70);

    fBranchingRatios.at(42).push_back(0.00816803);
    fBranchingRatios.at(42).push_back(0.0583431);
    fBranchingRatios.at(42).push_back(0.350058);
    fBranchingRatios.at(42).push_back(0.583431);
    fDecayTo.at(42).push_back(49);
    fDecayTo.at(42).push_back(62);
    fDecayTo.at(42).push_back(71);
    fDecayTo.at(42).push_back(72);

    fBranchingRatios.at(43).push_back(0.0961538);
    fBranchingRatios.at(43).push_back(0.423077);
    fBranchingRatios.at(43).push_back(0.480769);
    fDecayTo.at(43).push_back(66);
    fDecayTo.at(43).push_back(67);
    fDecayTo.at(43).push_back(68);

    fBranchingRatios.at(44).push_back(0.450549);
    fBranchingRatios.at(44).push_back(0.549451);
    fDecayTo.at(44).push_back(69);
    fDecayTo.at(44).push_back(72);

    fBranchingRatios.at(45).push_back(0.469925);
    fBranchingRatios.at(45).push_back(0.0836466);
    fBranchingRatios.at(45).push_back(0.446429);
    fDecayTo.at(45).push_back(61);
    fDecayTo.at(45).push_back(65);
    fDecayTo.at(45).push_back(72);

    fBranchingRatios.at(46).push_back(0.0408163);
    fBranchingRatios.at(46).push_back(0.510204);
    fBranchingRatios.at(46).push_back(0.44898);
    fDecayTo.at(46).push_back(67);
    fDecayTo.at(46).push_back(70);
    fDecayTo.at(46).push_back(71);

    fBranchingRatios.at(47).push_back(1);
    fDecayTo.at(47).push_back(72);

    fBranchingRatios.at(48).push_back(0.641026);
    fBranchingRatios.at(48).push_back(0.358974);
    fDecayTo.at(48).push_back(58);
    fDecayTo.at(48).push_back(69);

    fBranchingRatios.at(49).push_back(0.0458015);
    fBranchingRatios.at(49).push_back(0.954198);
    fDecayTo.at(49).push_back(66);
    fDecayTo.at(49).push_back(70);

    fBranchingRatios.at(50).push_back(0.401639);
    fBranchingRatios.at(50).push_back(0.188525);
    fBranchingRatios.at(50).push_back(0.409836);
    fDecayTo.at(50).push_back(61);
    fDecayTo.at(50).push_back(69);
    fDecayTo.at(50).push_back(72);

    fBranchingRatios.at(51).push_back(0.0188679);
    fBranchingRatios.at(51).push_back(0.173585);
    fBranchingRatios.at(51).push_back(0.754717);
    fBranchingRatios.at(51).push_back(0.0528302);
    fDecayTo.at(51).push_back(61);
    fDecayTo.at(51).push_back(67);
    fDecayTo.at(51).push_back(71);
    fDecayTo.at(51).push_back(72);

    fBranchingRatios.at(52).push_back(0.0100849);
    fBranchingRatios.at(52).push_back(0.00796178);
    fBranchingRatios.at(52).push_back(0.530786);
    fBranchingRatios.at(52).push_back(0.451168);
    fDecayTo.at(52).push_back(59);
    fDecayTo.at(52).push_back(68);
    fDecayTo.at(52).push_back(70);
    fDecayTo.at(52).push_back(71);

    fBranchingRatios.at(53).push_back(0.0379902);
    fBranchingRatios.at(53).push_back(0.0490196);
    fBranchingRatios.at(53).push_back(0.612745);
    fBranchingRatios.at(53).push_back(0.300245);
    fDecayTo.at(53).push_back(67);
    fDecayTo.at(53).push_back(70);
    fDecayTo.at(53).push_back(71);
    fDecayTo.at(53).push_back(72);

    fBranchingRatios.at(54).push_back(0.106557);
    fBranchingRatios.at(54).push_back(0.819672);
    fBranchingRatios.at(54).push_back(0.0737705);
    fDecayTo.at(54).push_back(59);
    fDecayTo.at(54).push_back(68);
    fDecayTo.at(54).push_back(70);

    fBranchingRatios.at(55).push_back(0.699301);
    fBranchingRatios.at(55).push_back(0.300699);
    fDecayTo.at(55).push_back(64);
    fDecayTo.at(55).push_back(70);

    fBranchingRatios.at(56).push_back(1);
    fDecayTo.at(56).push_back(71);

    fBranchingRatios.at(57).push_back(1);
    fDecayTo.at(57).push_back(72);

    fBranchingRatios.at(58).push_back(0.888099);
    fBranchingRatios.at(58).push_back(0.111901);
    fDecayTo.at(58).push_back(69);
    fDecayTo.at(58).push_back(72);

    fBranchingRatios.at(59).push_back(0.00647298);
    fBranchingRatios.at(59).push_back(0.752672);
    fBranchingRatios.at(59).push_back(0.165588);
    fBranchingRatios.at(59).push_back(0.0752672);
    fDecayTo.at(59).push_back(65);
    fDecayTo.at(59).push_back(70);
    fDecayTo.at(59).push_back(71);
    fDecayTo.at(59).push_back(72);

    fBranchingRatios.at(60).push_back(0.0708556);
    fBranchingRatios.at(60).push_back(0.668449);
    fBranchingRatios.at(60).push_back(0.260695);
    fDecayTo.at(60).push_back(65);
    fDecayTo.at(60).push_back(71);
    fDecayTo.at(60).push_back(72);

    fBranchingRatios.at(61).push_back(0.166667);
    fBranchingRatios.at(61).push_back(0.833333);
    fDecayTo.at(61).push_back(69);
    fDecayTo.at(61).push_back(72);

    fBranchingRatios.at(62).push_back(0.0898551);
    fBranchingRatios.at(62).push_back(0.57971);
    fBranchingRatios.at(62).push_back(0.330435);
    fDecayTo.at(62).push_back(67);
    fDecayTo.at(62).push_back(68);
    fDecayTo.at(62).push_back(70);

    fBranchingRatios.at(63).push_back(0.813008);
    fBranchingRatios.at(63).push_back(0.186992);
    fDecayTo.at(63).push_back(71);
    fDecayTo.at(63).push_back(72);

    fBranchingRatios.at(64).push_back(0.29078);
    fBranchingRatios.at(64).push_back(0.70922);
    fDecayTo.at(64).push_back(70);
    fDecayTo.at(64).push_back(71);

    fBranchingRatios.at(65).push_back(0.05);
    fBranchingRatios.at(65).push_back(0.08);
    fBranchingRatios.at(65).push_back(0.5);
    fBranchingRatios.at(65).push_back(0.37);
    fDecayTo.at(65).push_back(69);
    fDecayTo.at(65).push_back(70);
    fDecayTo.at(65).push_back(71);
    fDecayTo.at(65).push_back(72);

    fBranchingRatios.at(66).push_back(0.398406);
    fBranchingRatios.at(66).push_back(0.310757);
    fBranchingRatios.at(66).push_back(0.290837);
    fDecayTo.at(66).push_back(70);
    fDecayTo.at(66).push_back(71);
    fDecayTo.at(66).push_back(72);

    fBranchingRatios.at(67).push_back(0.819672);
    fBranchingRatios.at(67).push_back(0.180328);
    fDecayTo.at(67).push_back(70);
    fDecayTo.at(67).push_back(71);

    fBranchingRatios.at(68).push_back(0.186992);
    fBranchingRatios.at(68).push_back(0.813008);
    fDecayTo.at(68).push_back(70);
    fDecayTo.at(68).push_back(71);

    fBranchingRatios.at(69).push_back(1);
    fDecayTo.at(69).push_back(72);

    fBranchingRatios.at(70).push_back(1);
    fDecayTo.at(70).push_back(71);

    fBranchingRatios.at(71).push_back(1);
    fDecayTo.at(71).push_back(72);

    // Info from table VII (http://prc.aps.org/pdf/PRC/v58/i6/p3677_1)
    // in MeV
    double startEnergyLevels[] = {              2.281, 2.752, 2.937,
      3.143, 3.334, 3.569, 3.652, 3.786, 3.861, 4.067, 4.111, 4.267,
      4.364, 4.522, 4.655, 4.825, 5.017, 5.080, 5.223, 5.696, 6.006 };

    fStartEnergyLevels.insert(fStartEnergyLevels.end(), &startEnergyLevels[0],
                                     &startEnergyLevels[fNumberOfStartLevels]);

    double b[] = { 0.9,  1.5,  0.11, 0.06, 0.04, 0.01,
                   0.16, 0.26, 0.01, 0.05, 0.11, 0.29,
                   3.84, 0.31, 0.38, 0.47, 0.36, 0.23,
                                     0.03, 0.11, 0.13 };

    fB.insert(fB.end(), &b[0], &b[fNumberOfStartLevels]);

    double energyLevels[] = {       7.4724,    6.2270,    5.06347,
      4.99294,  4.8756,   4.87255,  4.78865,   4.744093,  4.53706,
      4.47299,  4.41936,  4.39588,  4.3840,    4.3656,    4.28052,
      4.25362,  4.21307,  4.18003,  4.14901,   4.11084,   4.10446,
      4.02035,  3.92390,  3.88792,  3.86866,   3.840228,  3.82143,
      3.79757,  3.76779,  3.73848,  3.663739,  3.62995,   3.59924,
      3.55697,  3.48621,  3.439144, 3.41434,   3.39363,   3.36803,
      3.22867,  3.15381,  3.14644,  3.12836,   3.109721,  3.1002,
      3.02795,  2.98587,  2.9508,   2.87901,   2.80788,   2.7874,
      2.786644, 2.75672,  2.74691,  2.730372,  2.625990,  2.57593,
      2.558,    2.54277,  2.419171, 2.397165,  2.290493,  2.289871,
      2.26040,  2.103668, 2.069809, 2.047354,  1.959068,  1.643639,
                          0.891398, 0.8001427, 0.0298299, 0.0      };

    fEnergyLevels.insert(fEnergyLevels.end(), &energyLevels[0],
                                &energyLevels[fNumberOfLevels]);

    // I have to put checks that the arrays really
    // have as many elements as they should

  }

  //----------------------------------------------------------------------------
  // Simulate particles and fill truth
  void NueAr40CCGenerator::CreateKinematicsVector(simb::MCTruth& truth,
                                        CLHEP::HepRandomEngine& engine) const
  {

    bool success = false;
    while (!success) {
      double neutrinoEnergy = GetNeutrinoEnergy(engine);
      double neutrinoTime   = GetNeutrinoTime  (engine);
      success = ProcessOneNeutrino(truth, neutrinoEnergy, neutrinoTime, engine);
    }

  }

  //----------------------------------------------------------------------------
  // Simulate particles and fill truth
  // Return true if one neutrino is processed, false otherwise
  bool NueAr40CCGenerator::ProcessOneNeutrino(simb::MCTruth& truth,
                        double neutrinoEnergy, double neutrinoTime,
                                    CLHEP::HepRandomEngine& engine) const
  {

    CLHEP::RandFlat randFlat(engine);

    int highestLevel = 0;
    std::vector< double > levelCrossSections =
                           CalculateCrossSections(neutrinoEnergy, highestLevel);

    double totalCrossSection = 0;
    // Calculating total cross section
    for (std::vector< double >::iterator crossSection =
                                                     levelCrossSections.begin();
                       crossSection != levelCrossSections.end(); ++crossSection)
      totalCrossSection += *crossSection;

    if (totalCrossSection == 0)
      return false;

    std::vector< double > startLevelProbabilities;
    // Calculating each level's probability
    for (std::vector< double >::iterator crossSection =
                                                     levelCrossSections.begin();
                       crossSection != levelCrossSections.end(); ++crossSection)
      startLevelProbabilities.push_back((*crossSection)/totalCrossSection);

    double randomNumber     = randFlat.fire();
    double tprob            =  0;
    int    chosenStartLevel = -1;
    // Picking a starting level
    for (int level = 0; level < highestLevel; ++level)
    {
      if (randomNumber < (startLevelProbabilities.at(level) + tprob))
      {
        chosenStartLevel = level;
        break;
      }
      tprob += startLevelProbabilities.at(level);
    }

    int lastLevel = -1;
    int level     = -1;

    // Time delays
    // Seems like this is not implemented yet
    std::vector< double > levelDelay;
    for (int n = 0; n < fNumberOfLevels; ++n) levelDelay.push_back(0.0);

    // The highest n for which fStartEnergyLevels.at(chosenStartLevel)
    // is higher than fEnergyLevels.at(n)
    int highestHigher = 0;
    // The lowest n for which fStartEnergyLevels.at(chosenStartLevel)
    // is lower than fEnergyLevels.at(n)
    int lowestLower   = 0;

    // Finding lowestLower and highestHigher
    for (int n = 0; n < fNumberOfLevels; ++n)
    {
      if (fStartEnergyLevels.at(chosenStartLevel) < fEnergyLevels.at(n))
        lowestLower = n;
      if (fStartEnergyLevels.at(chosenStartLevel) > fEnergyLevels.at(n))
      {
        highestHigher = n;
        break;
      }
    }

    if (std::abs(fStartEnergyLevels.at(chosenStartLevel) -
                                           fEnergyLevels.at(lowestLower))
          < std::abs(fStartEnergyLevels.at(chosenStartLevel)
                                       - fEnergyLevels.at(highestHigher)))
    {
      lastLevel = lowestLower;
      level     = lowestLower;
    }
    // If the chosen start level energy is closest to the lowest level energy
    // that it's lower than than the highest level energy
    // that it's higher than, it starts at the level of
    // the lowest level energy that it's lower than

    if (std::abs(fStartEnergyLevels.at(chosenStartLevel)
                                     - fEnergyLevels.at(highestHigher))
          < std::abs(fStartEnergyLevels.at(chosenStartLevel)
                                       - fEnergyLevels.at(lowestLower)))
    {
      lastLevel = highestHigher;
      level     = highestHigher;
    }
    // If the chosen start level energy is closest to the highest level energy
    // that it's higher than than the lowest level energy that it's lower than,
    // it starts at the level of the highest level energy that it's higher than

    std::vector< double > vertex = GetUniformPosition(engine);

    std::vector< double > neutrinoDirection = GetIsotropicDirection(engine);



    //
    // Add the first particle (the neutrino) to the truth list...
    //

    // NOTE: all of the values below calculated for the neutrino should be checked!

    // Primary particles have negative IDs
    int neutrinoTrackID = -1*(truth.NParticles() + 1);
    std::string primary("primary");

    int    nuePDG     = 12;
    double neutrinoP  = neutrinoEnergy/1000.0; // use GeV...
    double neutrinoPx = neutrinoDirection.at(0)*neutrinoP;
    double neutrinoPy = neutrinoDirection.at(1)*neutrinoP;
    double neutrinoPz = neutrinoDirection.at(2)*neutrinoP;



    // Create the neutrino:
    //   * set the mother to -1 since it is primary
    //   * set the mass to something < 0 so that the constructor looks up the mass from the pdg tables
    //   * set the status bit to 0 so that geant doesn't waste any CPU tracking the neutrino
    simb::MCParticle neutrino(neutrinoTrackID, nuePDG, primary, -1, -1, 0);

    TLorentzVector neutrinoPosition(vertex.at(0), vertex.at(1),
                                    vertex.at(2), neutrinoTime);
    TLorentzVector neutrinoMomentum(neutrinoPx, neutrinoPy,
                                    neutrinoPz, neutrinoEnergy/1000.0);
    neutrino.AddTrajectoryPoint(neutrinoPosition, neutrinoMomentum);

    truth.Add(neutrino);



    //
    // Adding the electron to truth
    //

    // In MeV
    double electronEnergy    = neutrinoEnergy -
                               (fEnergyLevels.at(lastLevel) + 1.5);
    double electronEnergyGeV = electronEnergy/1000.0;

    double electronM  = 0.000511;
    double electronP  = std::sqrt(std::pow(electronEnergyGeV,2)
                                         - std::pow(electronM,2));

    // For the moment, choose an isotropic direction for the electron.
    std::vector< double > electronDirection = GetIsotropicDirection(engine);
    double electronPx = electronDirection.at(0)*electronP;
    double electronPy = electronDirection.at(1)*electronP;
    double electronPz = electronDirection.at(2)*electronP;

    // Primary particles have negative IDs
    int trackID = -1*(truth.NParticles() + 1);
    int electronPDG = 11;
    simb::MCParticle electron(trackID, electronPDG, primary);

    TLorentzVector electronPosition(vertex.at(0), vertex.at(1),
                                    vertex.at(2), neutrinoTime);
    TLorentzVector electronMomentum(electronPx, electronPy,
                                    electronPz, electronEnergyGeV);
    electron.AddTrajectoryPoint(electronPosition, electronMomentum);

    truth.Add(electron);

    double ttime       = neutrinoTime;
    int    noMoreDecay = 0;

    // Level loop
    int groundLevel = fNumberOfLevels - 1;
    while (level != groundLevel)
    {

      double rl = randFlat.fire();

      int decayNum = 0;

      tprob = 0; // Used this variable above for cross section

      // Decay loop
      for (unsigned int iLevel = 0;
             iLevel < fBranchingRatios.at(level).size(); ++iLevel)
      {

        if (rl < (fBranchingRatios.at(level).at(iLevel) + tprob))
        {

          // We have a decay

          level = fDecayTo.at(level).at(decayNum);

          // Really should be using a map
          double gammaEnergy    = fEnergyLevels.at(lastLevel) -
                                                 fEnergyLevels.at(level);
          double gammaEnergyGeV = gammaEnergy/1000;

          std::vector< double > gammaDirection = GetIsotropicDirection(engine);
          double gammaPx = gammaDirection.at(0)*gammaEnergyGeV;
          double gammaPy = gammaDirection.at(1)*gammaEnergyGeV;
          double gammaPz = gammaDirection.at(2)*gammaEnergyGeV;

          //double gammaM  = 0.0; // unused

          double gammaTime = (-TMath::Log(randFlat.fire())/
                              (1/(levelDelay.at(lastLevel)))) + ttime;

          // Adding the gamma to truth

          // Primary particles have negative IDs
          trackID = -1*(truth.NParticles() + 1);
          int gammaPDG = 22;
          simb::MCParticle gamma(trackID, gammaPDG, primary); // NOTE: should the gammas be primary?

          TLorentzVector gammaPosition(vertex.at(0), vertex.at(1),
                                       vertex.at(2), neutrinoTime);
          TLorentzVector gammaMomentum(gammaPx, gammaPy,
                                       gammaPz, gammaEnergyGeV);
          gamma.AddTrajectoryPoint(gammaPosition, gammaMomentum);

          truth.Add(gamma);

          lastLevel = level;
          ttime     = gammaTime;

          break;

        }

        if ((tprob + fBranchingRatios.at(level).at(iLevel)) > 1) {
          std::cout << "(tprob + *iLevel) > 1" << std::endl;
          noMoreDecay = 1; // If it doesn't do any more gamma decay
          break;
        }

        ++decayNum;
        tprob += fBranchingRatios.at(level).at(iLevel);

      } // End of decay loop

      if (noMoreDecay == 1) break;

    } // End of level loop

    if (level != 72)
    {
      std::cout << "level != 72" << std::endl;
      std::cout << "level  = "   << level << std::endl;
    }

    // Set the neutrino for the MCTruth object:
    // NOTE: currently these parameters are all pretty much a guess...
    truth.SetNeutrino(simb::kCC,
		      simb::kQE, // not sure what the mode should be here, assumed that these are all QE???
		      simb::kCCQE, // not sure what the int_type should be either...
		      0, // target is AR40? not sure how to specify that...
		      0, // nucleon pdg ???
		      0, // quark pdg ???
		      -1.0,  // W ??? - not sure how to calculate this from the above
		      -1.0,  // X ??? - not sure how to calculate this from the above
		      -1.0,  // Y ??? - not sure how to calculate this from the above
		      -1.0); // Qsqr ??? - not sure how to calculate this from the above

    return true;

  }

  //----------------------------------------------------------------------------
  // Calculate cross sections for neutrino with neutrinoEnergy to excite
  // the final nucleus with highestLevel being the highest possible level.
  // highestLevel is output, so, whatever integer was used as the argument,
  // it may have a different value after this function is executed
  std::vector< double > NueAr40CCGenerator::CalculateCrossSections
                               (double neutrinoEnergy, int& highestLevel) const
  {

    highestLevel = 0;
    std::vector< double > levelCrossSections;

    // Loop through energy levels, if neutrino has enough energy,
    // calculate cross section
    for (int level = 0; level < fNumberOfStartLevels; ++level)
    {
      // Electron energy in keV
      double w = (neutrinoEnergy - (fStartEnergyLevels.at(level) + 1.5))*1000;

      double sigma = 0.0;
      if (neutrinoEnergy > (fStartEnergyLevels.at(level) + 1.5) && w >= 511.)
      {
        ++highestLevel;
        for (int n = 0; n <= level; ++n)
        {
          // Electron energy for level n to use in cross section
          w = (neutrinoEnergy - (fStartEnergyLevels.at(n) + 1.5))*1000;
          // Electron momentum in keV/c
          double p = std::sqrt(pow(w, 2) - pow(511.0, 2));
          // Fermi function approximation
          double f = std::sqrt(3.0634 + (0.6814/(w - 1)));
          // In cm^2*10^-42
          sigma += 1.6e-8*(p*w*f*fB.at(n));
        }
      }
      levelCrossSections.push_back(sigma);
    }

    return levelCrossSections;

  }

}
