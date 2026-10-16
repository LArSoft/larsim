/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : G4GPSAngDist_tool.cc
 * @created     : Friday Oct 16, 2026 15:34:16 CEST
 */

#ifndef LARSIM_G4GPSANGDIST_H

#define LARSIM_G4GPSANGDIST_H

#include "IAngDist.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/ToolConfigTable.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Geant4/G4SPSAngDistribution.hh"
#include "Geant4/G4SPSRandomGenerator.hh"

namespace evgen {
  class G4GPSAngDist : public IAngDist {
    public:
      struct Config {
        fhicl::Atom<std::string> type{fhicl::Name{"type"}, 
          fhicl::Comment{
            "Type of angular distribution to use in G4GPSAngDistribution\n"
            "Options: iso [default], cos, planar, beam1d, beam2d, focused, user"
          }, "iso"};
        fhicl::Sequence<int, 3> rot1{fhicl::Name{"rot1"}, 
          fhicl::Comment{
            "Defines the first (x' direction) rotation vector AR1 [default (1,0,0)]\n"
            "for the angular distribution and is not necessarily a unit vector.\n"
            "Used with /gps/ang/rot2 to compute the angular distribution rotation matrix."
          }, {1,0,0}};
        fhicl::Sequence<int, 3> rot2{fhicl::Name{"rot2"}, 
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
        fhicl::Atom<bool> user_coord{fhicl::Name{"user_coord"},
          fhicl::Comment{
            "Calculate the angular distribution with respect to the user defined co-ordinate system (true),\n"
            "or with respect to the global co-ordinate system (false, default)."
          }, false};
        fhicl::Atom<bool> surfnorm{fhicl::Name{"surfnorm"},
          fhicl::Comment{
            "Allows user to choose whether angular distributions are with respect to\n"
            "the co-ordinate system (false, default) or surface normals (true)\n"
            "for user-defined distributions."
          }, false};
      };

      using Parameters = art::ToolConfigTable<Config>;

      explicit G4GPSAngDist(Parameters const& config);
      geo::Vector_t generate(CLHEP::HepRandomEngine& engine) override;

    private:
      std::unique_ptr<G4SPSAngDistribution> fAngDist;
      std::unique_ptr<G4SPSRandomGenerator> fRandGen;
  };

    G4GPSAngDist::G4GPSAngDist(Parameters const& config)
    : fAngDist{std::make_unique<G4SPSAngDistribution>()},
      fRandGen{std::make_unique<G4SPSRandomGenerator>()}
  {
    fAngDist->SetAngDistType(config().type());
    fAngDist->SetMinTheta(config().minTheta());
    fAngDist->SetMaxTheta(config().maxTheta());
  }

  geo::Vector_t G4GPSAngDist::generate(CLHEP::HepRandomEngine& engine)
  {
    // TODO: point G4's global engine at `engine` (see earlier note on seeding).
    G4ParticleMomentum dir = fAngDist->GenerateOne();
    return {dir.x(), dir.y(), dir.z()};
  }
} // namespace evgen

DEFINE_ART_CLASS_TOOL(evgen::G4GPSAngDist)

#endif /* end of include guard LARSIM_G4GPSANGDIST_H */

