////////////////////////////////////////////////////////////////////////
/// \file  PhotonBackTracker.h
/// \brief back track the reconstruction to the simulation
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  jstock@fnal.gov
//  \adapted from BackTracker.h by brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef CHEAT_PHOTONBACKTRACKERER_H
#define CHEAT_PHOTONBACKTRACKERER_H
#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#else 
#define DEPRECATED
#endif



#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // int
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/Simulation/LArVoxelList.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/OpHit.h"

/*namespace recob{
  class SpacePoint;
}*/

///code to link reconstructed objects back to the MC truth information
namespace cheat{

  class PhotonBackTracker
  {

  public:

    PhotonBackTracker(fhicl::ParameterSet const& pset,
    art::ActivityRegistry&     reg);
    ~PhotonBackTracker();

    void reconfigure(fhicl::ParameterSet const& pset);

    // The Rebuild function rebuilds the various maps we need to answer backtracking queries.
    // It is called automatically before each event is processed. For jobs involving
    // Monte Carlo generation, this is too soon. So, we'll call rebuild after those data
    // products are put into the event in LArG4.  This is the least bad way of ensuring the
    // PhotonBackTracker works in jobs that combine MC production and reconstruction analysis based
    // on MC truth.  Don't copy this design pattern without talking to brebel@fnal.gov first
    void Rebuild(const art::Event& evt);

    // Get a reference to the ParticleList
    const sim::ParticleList&          ParticleList()  const { return fParticleList; }

    // Set the EveIdCalculator for the owned ParticleList
    void  SetEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }

    // Return a pointer to the simb::MCParticle object corresponding to
    // the given TrackID
    const simb::MCParticle*              TrackIDToParticle(int const& id)       const;
    const simb::MCParticle*              TrackIDToMotherParticle(int const& id) const;

    std::vector<sim::SDP>                TrackIDToSimSDP(int const& id)         const;

    // Get art::Ptr<> to simb::MCTruth and related information
    const art::Ptr<simb::MCTruth>&       TrackIDToMCTruth(int const& id)                        const;
    const art::Ptr<simb::MCTruth>&       ParticleToMCTruth(const simb::MCParticle* p)           const;
    std::vector<const simb::MCParticle*> MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const;
    const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector() const { return fMCTruthList; }

    // this method will return the Geant4 track IDs of 
    // the particles contributing ionization electrons to the identified hit
    DEPRECATED std::vector<sim::TrackSDP> OpHitToTrackID(art::Ptr<recob::OpHit> const& hit)
      {return OpHitToTrackSDPs(hit);}
    std::vector<sim::TrackSDP> OpHitToTrackSDPs(art::Ptr<recob::OpHit> const& hit);
    
    // method to return a subset of allhits that are matched to a list of TrackIDs
    const std::vector<std::vector<art::Ptr<recob::OpHit>>> TrackIDsToOpHits(std::vector<art::Ptr<recob::OpHit>> const& allhits,
                  std::vector<int> const& tkIDs);
    
    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    std::vector<sim::TrackSDP> OpHitToEveSDPs(art::Ptr<recob::OpHit> const& hit);
    std::vector<sim::TrackSDP> OpHitToEveID(art::Ptr<recob::OpHit> const& hit);
    
    //@{
    // method to return sim::SDP objects associated with a given hit
    void                 OpHitToSDPs(recob::OpHit const& hit,
                                      std::vector<sim::SDP>&      ides) const;
    DEPRECATED void      OpHitToSimSDPs(recob::OpHit const& hit,
                                      std::vector<sim::SDP>&      ides) const
                                      { OpHitToSDPs( hit, ides); }
    void                 OpHitToSDPs(art::Ptr<recob::OpHit> const& hit,
                                      std::vector<sim::SDP>&      ides) const
                                      { OpHitToSDPs(*hit, ides); }
    DEPRECATED void      OpHitToSimSDPs(art::Ptr<recob::OpHit> const& hit,
                                      std::vector<sim::SDP>&      ides) const
                                      { OpHitToSDPs(*hit, ides); }
    //@}
    
    // method to return the XYZ position of the weighted average energy deposition for a given hit
    std::vector<double>  SimSDPsToXYZ(std::vector<sim::SDP> const& ides);
    
    // method to return the XYZ position of the weighted average energy deposition for a given hit
    std::vector<double>  OpHitToXYZ(art::Ptr<recob::OpHit> const& hit);          
    
    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
/*    std::vector<double> SpacePointToXYZ(art::Ptr<recob::SpacePoint> const& spt,
          art::Event                  const& evt,
          std::string                 const& label);*/

    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
/*    std::vector<double> SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits);*/
    
    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids 
    double              OpHitCollectionPurity(std::set<int>                              trackIDs, 
              std::vector< art::Ptr<recob::OpHit> > const& hits);
    
    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are 
    // represented in a collection of hits
    double              OpHitCollectionEfficiency(std::set<int>                              trackIDs, 
            std::vector< art::Ptr<recob::OpHit> > const& hits,
            std::vector< art::Ptr<recob::OpHit> > const& allhits);
    double              OpHitCollectionEfficiency(std::set<int>                              trackIDs, 
            std::vector< art::Ptr<recob::OpHit> > const& hits,
            std::vector< art::Ptr<recob::OpHit> > const& allhits,
            geo::View_t                         const& view);

    // method to return the fraction of charge in a collection that come from the specified Geant4 track ids 
    double              OpHitChargeCollectionPurity(std::set<int>                              trackIDs, 
              std::vector< art::Ptr<recob::OpHit> > const& hits);
    
    // method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are 
    // represented in a collection of hits
    double              OpHitChargeCollectionEfficiency(std::set<int>                              trackIDs, 
                  std::vector< art::Ptr<recob::OpHit> > const& hits,
                  std::vector< art::Ptr<recob::OpHit> > const& allhits);
    double              OpHitChargeCollectionEfficiency(std::set<int>                              trackIDs, 
                  std::vector< art::Ptr<recob::OpHit> > const& hits,
                  std::vector< art::Ptr<recob::OpHit> > const& allhits,
                  geo::View_t                         const& view);
  
    // method to return all EveIDs corresponding to the current sim::ParticleList
    std::set<int>       GetSetOfEveIDs();

    // method to return all TrackIDs corresponding to the current sim::ParticleList
    std::set<int>       GetSetOfTrackIDs();

    // method to return all EveIDs corresponding to the given list of hits
    std::set<int>       GetSetOfEveIDs(std::vector< art::Ptr<recob::OpHit> > const& hits);

    // method to return all TrackIDs corresponding to the given list of hits
    std::set<int>       GetSetOfTrackIDs(std::vector< art::Ptr<recob::OpHit> > const& hits);

    const std::vector< art::Ptr< sim::OpDetBacktrackerRecord >>& OpDetBacktrackerRecords() const { return cOpDetBacktrackerRecords; } 

    void ChannelToTrackSDPs(std::vector<sim::TrackSDP>& trackSDPs,
        int channel,
        const double hit_start_time,
        const double hit_end_time);
    
  private:

    geo::GeometryCore const* geom = lar::providerFrom<geo::Geometry>();

    const art::Ptr< sim::OpDetBacktrackerRecord > FindOpDetBacktrackerRecord(int channel) const;

    const void shouldThisFail() const;

    bool have_complained;

    sim::ParticleList                      fParticleList;          ///< ParticleList to map track ID to sim::Particle
    sim::LArVoxelList                      fVoxelList;             ///< List to map the position of energy depostions
                                                                   ///< in voxels to the particles depositing the 
                                                                   ///< energy
    std::vector< art::Ptr<simb::MCTruth> > fMCTruthList;           ///< all the MCTruths for the event
    std::vector< art::Ptr< sim::OpDetBacktrackerRecord >>  cOpDetBacktrackerRecords;         ///< all the OpDetBacktrackerRecords for the event
    std::map<int, int>                     fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry
    std::string                            fG4ModuleLabel;         ///< label for geant4 module
    double                                 fMinOpHitEnergyFraction;  ///< minimum fraction of energy a track id has to 
    double                                 fDelay;                 //Shift in time
                                                                   ///< contribute to a hit to be counted in
                                                                   ///< purity and efficiency calculations 
                                                                   ///< based on hit collections
  };
} // namespace

DECLARE_ART_SERVICE(cheat::PhotonBackTracker, LEGACY)
#endif // CHEAT_PHOTONBACKTRACKER_H
