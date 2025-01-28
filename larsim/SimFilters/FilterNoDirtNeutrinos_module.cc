////////////////////////////////////////////////////////////////////////
/// \file  FilterNoDirtNeutrinos_module.cc
/// \brief Simple EDFilter to require neutrino interaction in TPC
///
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// C++ Includes
#include <cstring>
#include <iostream>

///Geant4 interface
namespace simfilter {

  class FilterNoDirtNeutrinos : public art::EDFilter {
  public:
    explicit FilterNoDirtNeutrinos(fhicl::ParameterSet const& pset);

  private:
    bool filter(art::Event&) override;

    std::string fLArG4ModuleLabel;
    std::string fGenModuleLabel;
    bool fKeepCryostatNeutrinos;
  };

} // namespace simfilter

namespace simfilter {

  //-----------------------------------------------------------------------
  // Constructor
  FilterNoDirtNeutrinos::FilterNoDirtNeutrinos(fhicl::ParameterSet const& pset)
    : EDFilter{pset}
    , fLArG4ModuleLabel(pset.get<std::string>("LArG4ModuleLabel", "NoLabel"))
    , fGenModuleLabel(pset.get<std::string>("GenModuleLabel", "NoLabel"))
    , fKeepCryostatNeutrinos(pset.get<bool>("KeepCryostatNeutrinos", false))
  {}

  //-----------------------------------------------------------------------
  bool FilterNoDirtNeutrinos::filter(art::Event& evt)
  {
    bool interactionDesired(false);
    //get the list of particles from this event
    art::ServiceHandle<geo::Geometry const> geom;
    auto const& [cryostat, tpc] = std::make_tuple(geom->Cryostat(), geom->TPC());

    // * MC truth information
    art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth>> mclist;
    art::Handle<std::vector<simb::MCParticle>> mcpHandle;

    if (evt.getByLabel(fGenModuleLabel, mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);
    evt.getByLabel(fLArG4ModuleLabel, mcpHandle);

    std::set<art::Ptr<simb::MCTruth>> mctSetGENIE;
    for (size_t i = 0; i < mctruthListHandle->size(); ++i) {
      art::Ptr<simb::MCTruth> mct_ptr(mctruthListHandle, i);
      if (mctSetGENIE.find(mct_ptr) == mctSetGENIE.end()) mctSetGENIE.insert(mct_ptr);
    }

    // Get the MCTruths from associations to our particles
    art::FindOneP<simb::MCTruth> assMCT(mcpHandle, evt, "largeant");

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    if (fKeepCryostatNeutrinos) {
      // Get cryostat (box) volume boundary.
      xmin = tpc.HalfWidth() - cryostat.HalfWidth();
      xmax = tpc.HalfWidth() + cryostat.HalfWidth();
      ymin = -cryostat.HalfHeight();
      ymax = cryostat.HalfHeight();
      zmin = tpc.Length() / 2. - tpc.Length() / 2.;
      zmax = tpc.Length() / 2. + tpc.Length() / 2.;
    }
    else {
      // Get fiducial volume boundary.
      xmin = 0.;
      xmax = 2. * tpc.HalfWidth();
      ymin = -tpc.HalfHeight();
      ymax = tpc.HalfHeight();
      zmin = 0.;
      zmax = tpc.Length();
    }

    // Now let's loop over G4 MCParticle list and track back MCTruth
    bool inTPC(false);
    for (size_t i = 0; i < mcpHandle->size() && !inTPC; ++i) {
      const art::Ptr<simb::MCParticle> mcp_ptr(mcpHandle, i);
      const art::Ptr<simb::MCTruth>& mct = assMCT.at(i);
      if (mctSetGENIE.find(mct) == mctSetGENIE.end()) {
        // This is non-genie
        continue;
      }
      else {
        // This is genie

        const simb::MCParticle* part(&mcpHandle->at(i));

        // Now walk through trajectory and see if it enters the TPC
        int n = part->NumberTrajectoryPoints();
        for (int j = 0; j < n && !inTPC; ++j) {

          TVector3 pos = part->Position(j).Vect();
          if (pos.X() >= xmin && pos.X() <= xmax && pos.Y() >= ymin && pos.Y() <= ymax &&
              pos.Z() >= zmin && pos.Z() <= zmax) {
            interactionDesired = true;
            std::cout << "FilterNoDirtNeutrinos: Genie daughter found in TPC. G4Particle "
                      << std::endl;
            inTPC = true;
          }
        } // trajectory loop
      }   // end Genie particle
    }     // loop on MCPHandle

    return interactionDesired;

  } // end FilterNoDirtNeutrinos()function

} // namespace simfilter

DEFINE_ART_MODULE(simfilter::FilterNoDirtNeutrinos)
