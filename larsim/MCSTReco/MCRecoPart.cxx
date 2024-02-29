////////////////////////////////////////////////////////////////////////
//
//  MCRecoPart source
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "MCRecoPart.h"

namespace sim {

  //--------------------------------------------------------------------------------------------
  MCRecoPart::MCRecoPart(fhicl::ParameterSet const& pset)
    : _x_max{std::numeric_limits<double>::min()}
    , _x_min{std::numeric_limits<double>::max()}
    , _y_max{std::numeric_limits<double>::min()}
    , _y_min{std::numeric_limits<double>::max()}
    , _z_max{std::numeric_limits<double>::min()}
    , _z_min{std::numeric_limits<double>::max()}
    , _trackIDOffsets{pset.get<std::vector<unsigned int>>("TrackIDOffsets", {0})}
  {
    this->clear();
    _track_index.clear();
    _pdg_list.clear();
    for (auto const& id : pset.get<std::vector<int>>("SavePathPDGList"))

      _pdg_list.insert(id);

    art::ServiceHandle<geo::Geometry const> geo;
    // Build "Fiducial" Volume Definition:
    //
    // Iterate over all TPC's to get bounding box that covers volumes of each individual TPC in the detector
    for (auto const& tpc : geo->Iterate<geo::TPCGeo>()) {
      _x_max = std::max(_x_max, tpc.BoundingBox().MaxX());
      _x_min = std::min(_x_min, tpc.BoundingBox().MinX());
      _y_max = std::max(_y_max, tpc.BoundingBox().MaxY());
      _y_min = std::min(_y_min, tpc.BoundingBox().MinY());
      _z_max = std::max(_z_max, tpc.BoundingBox().MaxZ());
      _z_min = std::min(_z_min, tpc.BoundingBox().MinZ());
    }
  }

  //--------------------------------------------------------------------------------------------
  unsigned int MCRecoPart::MotherTrackID(const unsigned int part_index) const
  //--------------------------------------------------------------------------------------------
  {
    if (this->size() <= part_index) return ::sim::kINVALID_UINT;

    unsigned int result = this->at(part_index)._mother;

    //Is a primary particle. Account for possible track ID offsets
    for (auto const& offset : _trackIDOffsets) {
      if (result == offset) return this->at(part_index)._track_id;
    }

    if (TrackToParticleIndex(result) != ::sim::kINVALID_UINT) return result;

    //std::cout<< "\033[95mWarning:\033[00m Mother particle not in the particle list!"<<std::endl;
    // Do brute search
    unsigned int daughter_id = this->at(part_index)._track_id;

    for (auto const& part : *this) {

      if (part.HasDaughter(daughter_id)) return part._track_id;
    }
    return result;
  }

  //--------------------------------------------------------------------------------------------
  unsigned int MCRecoPart::AncestorTrackID(const unsigned int part_index)
  //--------------------------------------------------------------------------------------------
  {
    if (part_index >= this->size()) return kINVALID_UINT;

    if ((*this)[part_index]._ancestor != kINVALID_UINT) return (*this)[part_index]._ancestor;

    auto result = MotherTrackID(part_index);

    if (result == this->at(part_index)._track_id) return result;

    //Is a primary particle. Account for possible track ID offsets
    for (auto const& offset : _trackIDOffsets) {
      if (result == offset) return this->at(part_index)._track_id;
    }

    auto mother_index = TrackToParticleIndex(result);

    while (1) {

      if (mother_index != kINVALID_UINT) {

        auto const new_result = MotherTrackID(mother_index);

        if (new_result == this->at(mother_index)._track_id) break;

        result = new_result;
      }
      else {

        // Look for a particle that has a daughter = this mother
        auto const old_result = result;
        for (auto const& p : *this) {

          if (p.HasDaughter(result)) {
            result = p._track_id;
            break;
          }
        }
        if (result == old_result) break;
      }

      mother_index = TrackToParticleIndex(result);
    }

    (*this)[part_index]._ancestor = result;
    return result;
  }

  //--------------------------------------------------------------------------------------------
  bool MCRecoPart::InDetector(const double& x, const double& y, const double& z) const
  //--------------------------------------------------------------------------------------------
  {
    return !(x > _x_max || x < _x_min || z > _z_max || z < _z_min || y > _y_max || y < _y_min);
  }

  //--------------------------------------------------------------------------------------------
  void MCRecoPart::AddParticles(const std::vector<simb::MCParticle>& mcp_v,
                                const std::vector<simb::Origin_t>& orig_v,
                                const std::vector<sim::MCParticleLite>& mcmp_v)
  //--------------------------------------------------------------------------------------------
  {
    if (orig_v.size() != mcp_v.size())
      throw cet::exception(__FUNCTION__) << "MCParticle and Origin_t vector size not same!";

    this->clear();
    _track_index.clear();

    for (size_t i = 0; i < mcp_v.size(); ++i) {

      auto const& mcp = mcp_v[i];

      //std::cout<<" Track ID : "<<mcp.TrackId()<<" ... Index : " <<this->size()<<std::endl;

      _track_index.insert(std::make_pair((size_t)(mcp.TrackId()), (size_t)(this->size())));

      // Change units to LArSoft (MeV, cm, us)
      // (done inside constructor of MCMiniPart)
      this->push_back(MCMiniPart(mcp));

      auto& mini_mcp = (*this->rbegin());

      for (size_t i = 0; i < (size_t)(mcp.NumberDaughters()); ++i) {
        mini_mcp.AddDaughter(mcp.Daughter(i));
      }
      mini_mcp._origin = orig_v[i];

      if (_pdg_list.find(mcp.PdgCode()) != _pdg_list.end()) {

        std::set<size_t> det_path_index;

        for (size_t i = 0; i < mcp.NumberTrajectoryPoints(); ++i) {

          if (InDetector(mcp.Vx(i), mcp.Vy(i), mcp.Vz(i))) det_path_index.insert(i);
        }

        if (det_path_index.size()) {
          if ((*det_path_index.begin())) det_path_index.insert((*det_path_index.begin()) - 1);
          if (det_path_index.size() > 1) {
            if (((*det_path_index.rbegin()) + 1) < mcp.NumberTrajectoryPoints())
              det_path_index.insert((*det_path_index.rbegin()) + 1);
          }
          std::vector<std::pair<TLorentzVector, TLorentzVector>> det_path;
          det_path.reserve(det_path_index.size());
          for (auto const& index : det_path_index) {

            TLorentzVector vec(mcp.Momentum(index));
            for (size_t i = 0; i < 4; ++i)
              vec[i] *= 1.e3;

            det_path.emplace_back(mcp.Position(index), vec);
          }
          mini_mcp._det_path = std::move(det_path);
        }
      } // end if in _pdg_list
    }   // end for loop over mcp_v

    // Now loop over dropped particles
    for (auto const& mcmp : mcmp_v) {

      _track_index.try_emplace(mcmp.TrackID(), this->size());

      this->push_back(sim::MCMiniPart(mcmp));

    } // end for loop over mcmp_v
  }   // end AddParticles
}
