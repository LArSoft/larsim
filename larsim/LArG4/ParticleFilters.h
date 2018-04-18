/**
 * @file    ParticleFilters.h
 * @brief   Defines classes to filter particles based on their trajectory.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 24, 2016
 *
 *
 */

#ifndef LARSIM_LARG4_PARTICLEFILTERS_H
#define LARSIM_LARG4_PARTICLEFILTERS_H

// LArSoft libraries

// ROOT libraries
#include "TGeoVolume.h"
#include "TGeoMatrix.h" // TGeoCombiTrans
#include "TLorentzVector.h"
#include "TVector3.h"

// C/C++ standard libraries
#include <array>
#include <vector>
#include <utility> // std::move()


namespace larg4 {

  /**
   * @brief Tag for filters
   * 
   * A filter of this type acts on trajectory coordinates, and if it ever
   * returns true, then the particle is safe.
   */
  struct KeepByPositionFilterTag {};


  /** **************************************************************************
   * @brief Use to keep particles with at least part of trajectory in a volume
   * @author Matt Bass, Gianluca Petrillo
   * 
   * The class stores a list of pointers to volumes to consider.
   * If a point specified in the mustKeep() call is within one of the blessed
   * volumes, the whole track it belongs to must to be kept.
   * 
   * No condition for prompt rejection is provided.
   */
  class PositionInVolumeFilter: public KeepByPositionFilterTag {
      public:
    
    using Point_t = std::array<double, 3>;
    
    /// Due to the structure of ROOT geometry, volumes and their transformations
    /// are not living in the same place; so we need to keep both.
    struct VolumeInfo_t {
      VolumeInfo_t(TGeoVolume const* new_vol, TGeoCombiTrans const* new_trans)
        : vol(new_vol), trans(new_trans) {}
      
      TGeoVolume const*     vol;   ///< ROOT volume
      TGeoCombiTrans const* trans; ///< volume transformation (has both ways)
    }; // VolumeInfo_t
    
    using AllVolumeInfo_t = std::vector<VolumeInfo_t>;
    
    /// @{
    /// @brief Constructors: read the volumes from the specified list
    /// @param volumes list of interesting volumes
    PositionInVolumeFilter(std::vector<VolumeInfo_t> const& volumes)
      : volumeInfo(volumes)
      {}
    PositionInVolumeFilter(std::vector<VolumeInfo_t>&& volumes)
      : volumeInfo(std::move(volumes))
      {}
    /// @}
    
    /**
     * @brief Returns whether a track along the specified point must be kept
     * @param pos point on the track, a [x,y,z] array in "Geant4 coordinates"
     * @brief whether a track along the specified point must be kept
     * 
     * If the return value is true, the track must be kept.
     * If the return value is false, no decision can be made yet.
     * If no decision is made after the track is over, the track should be
     * dropped.
     */
    bool mustKeep(Point_t const& pos) const
      {
        // if no volume is specified, it means we don't filter
        if (volumeInfo.empty()) return true;
        double local[3];
        for(auto const& info: volumeInfo) {
          // transform the point to relative to the volume
          info.trans->MasterToLocal(pos.data(), local);
          // containment check
          if (info.vol->Contains(local)) return true;
        } // for volumes
        return false;
      } // mustKeep()
    
    bool mustKeep(TVector3 const& pos) const
      { return mustKeep(Point_t{{ pos.X(), pos.Y(), pos.Z() }}); }

    bool mustKeep(TLorentzVector const& pos) const
      { return mustKeep(Point_t{{ pos.X(), pos.Y(), pos.Z() }}); }

    
      protected:
    std::vector<VolumeInfo_t> volumeInfo; ///< all good volumes
    
  }; // PositionInVolumeFilter
  
} // namespace larg4

#endif // LARSIM_LARG4_PARTICLEFILTERS_H
