/**
 * @author      : Daniele Guffanti (University and INFN Milano-Bicocca)
 * @file        : G4GPSAngDist.h
 * @created     : Friday Oct 16, 2026 15:34:16 CEST
 */

#ifndef LARSIM_G4GPSANGDIST_H

#define LARSIM_G4GPSANGDIST_H

#include "DirectionGenerator.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"
#include "canvas/Utilities/Exception.h"

#include "Geant4/G4SPSAngDistribution.hh"
#include "Geant4/G4SPSRandomGenerator.hh"
#include "Geant4/G4SPSPosDistribution.hh"
#include "Geant4/G4RandomTools.hh"

namespace evgen {
  namespace {
    class G4EngineGuard {
    public:
      explicit G4EngineGuard(CLHEP::HepRandomEngine& engine)
        : fPrevious{G4Random::getTheEngine()}
      {
        G4Random::setTheEngine(&engine);
      }
      ~G4EngineGuard() { G4Random::setTheEngine(fPrevious); }
      G4EngineGuard(G4EngineGuard const&) = delete;
      G4EngineGuard& operator=(G4EngineGuard const&) = delete;
    private:
      CLHEP::HepRandomEngine* fPrevious = {};
    };
  }

  class G4GPSAngDist : public IDirectionGenerator {
    public:
      struct Config {
        fhicl::Atom<std::string> type{fhicl::Name{"distType"}, 
          fhicl::Comment{
            "Type of angular distribution to use in G4GPSAngDistribution\n"
            "Options: iso [default], cos, planar, beam1d, beam2d, focused, user"
          }, "iso"};
        fhicl::Sequence<double, 3> rot1{fhicl::Name{"rot1"}, 
          fhicl::Comment{
            "Defines the first (x' direction) rotation vector AR1 [default (1,0,0)]\n"
            "for the angular distribution and is not necessarily a unit vector.\n"
            "Used with /gps/ang/rot2 to compute the angular distribution rotation matrix."
          }, {1,0,0}};
        fhicl::Sequence<double, 3> rot2{fhicl::Name{"rot2"}, 
          fhicl::Comment{
            "Defines the second rotation vector AR2 [default (0,1,0)]\n"
            "for the angular distribution and is not necessarily a unit vector.\n"
            "Used with /gps/ang/rot1 to compute the angular distribution rotation matrix."
          }, {0,1,0}};
        fhicl::Atom<double> minTheta{fhicl::Name{"minTheta"},
          fhicl::Comment{
            "Sets a minimum value [default 0 rad] for the θ distribution. Units are rad."
          }, 0.};
        fhicl::Atom<double> maxTheta{fhicl::Name{"maxTheta"},
          fhicl::Comment{
            "Sets a maximum value [default π rad] for the θ distribution. Units are rad."
          }, CLHEP::pi};
        fhicl::Atom<double> minPhi{fhicl::Name{"minPhi"},
          fhicl::Comment{
            "Sets a minimum value [default 0 rad] for the φ distribution. Units are rad."
          }, 0.};
        fhicl::Atom<double> maxPhi{fhicl::Name{"maxPhi"},
          fhicl::Comment{
            "Sets a maximum value [default 2π rad] for the φ distribution. Units are rad."
          }, 2*CLHEP::pi};
        fhicl::Atom<double> sigma_r{fhicl::Name{"sigma_r"},
          fhicl::Comment{
            "Sets the standard deviation [default 0 rad] of beam directional profile in radial. Units are rad."
          }, 0.};
        fhicl::Atom<double> sigma_x{fhicl::Name{"sigma_x"},
          fhicl::Comment{
            "Sets the standard deviation [default 0 rad] of beam directional profile in x. Units are rad."
          }, 0.};
        fhicl::Atom<double> sigma_y{fhicl::Name{"sigma_y"},
          fhicl::Comment{
            "Sets the standard deviation [default 0 rad] of beam directional profile in y. Units are rad."
          }, 0.};
        fhicl::Sequence<double,3> focuspoint{fhicl::Name{"focuspoint"},
          fhicl::Comment{
            "Sets the focus point [default 0,0,0] of the beam. Units are mm."
          }, {0.,0.,0.}};
        fhicl::Atom<bool> surfnorm{fhicl::Name{"surfnorm"},
          fhicl::Comment{
            "Allows user to choose whether angular distributions are with respect to\n"
            "the coordinate system (false, default) or surface normals (true)\n"
            "for user-defined distributions."
          }, false};
        fhicl::Atom<int> verbosity{fhicl::Name{"verbosity"},
          fhicl::Comment{
            "Sets the verbosity level [default 0] for the G4GPSAngDist tool."
          }, 0};
      };

      using Parameters = art::ToolConfigTable<Config>;

      explicit G4GPSAngDist(Parameters const& config);

      geo::Vector_t generate(CLHEP::HepRandomEngine& engine) override;

    protected:
      void dumpSettingsImpl(std::ostream& out) const override;


    private:
      G4ThreeVector fAngRef1 = {};
      G4ThreeVector fAngRef2 = {};

      std::unique_ptr<G4SPSAngDistribution> fAngDist;
      std::unique_ptr<G4SPSRandomGenerator> fRandGen;
      std::unique_ptr<G4SPSPosDistribution> fPosDist;
  };

    G4GPSAngDist::G4GPSAngDist(Parameters const& config)
    : fAngDist{std::make_unique<G4SPSAngDistribution>()},
      fRandGen{std::make_unique<G4SPSRandomGenerator>()},
      fPosDist{std::make_unique<G4SPSPosDistribution>()}
  {
    // Set verbosity
    fAngDist->SetVerbosity(config().verbosity());
    
    // Validate type
    static const std::array<std::string, 7> validTypes{
      "iso", "cos", "planar", "beam1d", "beam2d", "focused", "user"};
    const std::string type = config().type();
    if (std::find(validTypes.begin(), validTypes.end(), type) == validTypes.end()) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: invalid 'type' \"" << type << "\".\n"
        << "Valid options: iso, cos, planar, beam1d, beam2d, focused, user.";
    }
    fAngDist->SetAngDistType( type );

    // Validate angular ranges
    if (config().minTheta() < 0. || config().maxTheta() > CLHEP::pi
        || config().minTheta() >= config().maxTheta()) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: invalid theta range ["
        << config().minTheta() << ", " << config().maxTheta()
        << "] rad. Must satisfy 0 <= minTheta < maxTheta <= pi.";
    }
    if (config().minPhi() < 0. || config().maxPhi() > 2*CLHEP::pi
        || config().minPhi() >= config().maxPhi()) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: invalid phi range ["
        << config().minPhi() << ", " << config().maxPhi()
        << "] rad. Must satisfy 0 <= minPhi < maxPhi <= 2pi.";
    }
    fAngDist->SetMinTheta(config().minTheta());
    fAngDist->SetMaxTheta(config().maxTheta());
    fAngDist->SetMinPhi(config().minPhi());
    fAngDist->SetMaxPhi(config().maxPhi());

    // Rotation vectors must be non-zero
    auto isZero = [](auto const& v) { return std::abs(v[0]) < 1e-8 && std::abs(v[1]) < 1e-8 && std::abs(v[2]) < 1e-8; };
    if (isZero(config().rot1()) || isZero(config().rot2())) {
      throw art::Exception(art::errors::Configuration)
        << "G4GPSAngDist: 'rot1' and 'rot2' must be non-zero vectors.";
    }
    fAngRef1.set(config().rot1()[0], config().rot1()[1], config().rot1()[2]);
    fAngDist->DefineAngRefAxes("angref1", 
        G4ThreeVector(config().rot1()[0], config().rot1()[1], config().rot1()[2]));
    fAngRef2.set(config().rot2()[0], config().rot2()[1], config().rot2()[2]);
    fAngDist->DefineAngRefAxes("angref2",
        G4ThreeVector(config().rot2()[0], config().rot2()[1], config().rot2()[2]));

    fAngDist->SetBeamSigmaInAngR(config().sigma_r());
    fAngDist->SetBeamSigmaInAngX(config().sigma_x());
    fAngDist->SetBeamSigmaInAngY(config().sigma_y());
    fAngDist->SetFocusPoint(
        G4ThreeVector(config().focuspoint()[0], config().focuspoint()[1], config().focuspoint()[2]));
    fAngDist->SetUserWRTSurface(config().surfnorm());

    // Set the random generator for pos and ang distributions to use the same engine
    fPosDist->SetBiasRndm( fRandGen.get() );
    fAngDist->SetBiasRndm( fRandGen.get() );
    fAngDist->SetPosDistribution( fPosDist.get() );
  }

  geo::Vector_t G4GPSAngDist::generate(CLHEP::HepRandomEngine& engine)
  {
    G4EngineGuard guard{engine};
    G4ParticleMomentum dir = fAngDist->GenerateOne();
    return {dir.x(), dir.y(), dir.z()};
  }

  void G4GPSAngDist::dumpSettingsImpl(std::ostream& out) const
  {
    out << "G4GPSAngDist settings:\n"
      << "  type: " << fAngDist->GetDistType() << "\n"
      << "  minTheta: " << fAngDist->GetMinTheta() << " rad\n"
      << "  maxTheta: " << fAngDist->GetMaxTheta() << " rad\n"
      << "  minPhi: " << fAngDist->GetMinPhi() << " rad\n"
      << "  maxPhi: " << fAngDist->GetMaxPhi() << " rad\n"
      << "  rot1: (" << fAngRef1.x() << ", " << fAngRef1.y() << ", " << fAngRef1.z() << ")\n"
      << "  rot2: (" << fAngRef2.x() << ", " << fAngRef2.y() << ", " << fAngRef2.z() << ")\n";
  }
} // namespace evgen

#endif /* end of include guard LARSIM_G4GPSANGDIST_H */

