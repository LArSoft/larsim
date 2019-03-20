/**
 * @file  larsim/Simulation/PhotonVoxels.h
 * @brief Definitions of voxel data structures.
 */

#ifndef LARSIM_SIMULATION_PHOTONVOXELS_H
#define LARSIM_SIMULATION_PHOTONVOXELS_H

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C/C++ standard libraries
#include <array>
#include <optional>


namespace sim {


  /// Representation of a single small volume (voxel).
  class PhotonVoxel{
  public:
    PhotonVoxel() = default;
    PhotonVoxel(geo::Point_t const& min, geo::Point_t const& max)
      : fVoxelMin(min), fVoxelMax(max) {}
    PhotonVoxel(double xMin, 
                double xMax, 
                double yMin, 
                double yMax, 
                double zMin, 
                double zMax)
      : PhotonVoxel({ xMin, yMin, zMin }, { xMax, yMax, zMax }) {}

  private:
    geo::Point_t fVoxelMin;
    geo::Point_t fVoxelMax;
    
  public:
    using DefaultPoint = geo::Point_t;
    
    /// @{
    // the choice of `decltype(auto)` is because in case `geo::Point_t` is the
    // requested `Point` type, a reference to the data member is returned
    // instead of a copy.
    
    /// Returns the voxel vertex (type `Point`) with the lowest coordinates.
    template <typename Point = DefaultPoint>
    decltype(auto) GetLowerCorner() const;
    
    /// Returns the voxel vertex (type `Point`) with the highest coordinates.
    template <typename Point = DefaultPoint>
    decltype(auto) GetUpperCorner() const;
    
    /// Returns the center of the voxel (type `Point`).
    template <typename Point = DefaultPoint>
    Point GetCenter() const;

    /// @}
    
  }; // class PhotonVoxel

  
  /// Representation of a region of space diced into voxels.
  class PhotonVoxelDef
  {
    using DefaultPoint = geo::Point_t;
    using DefaultVector = geo::Vector_t;
    
    geo::Point_t fLowerCorner;
    geo::Point_t fUpperCorner;
    unsigned int fxSteps = 1U;
    unsigned int fySteps = 1U;
    unsigned int fzSteps = 1U;
    
  public:
    PhotonVoxelDef() = default;
    PhotonVoxelDef(double xMin, 
                   double xMax, 
                   int xN, 
                   double yMin, 
                   double yMax, 
                   int yN, 
                   double zMin, 
                   double zMax, 
                   int z);
    
    /// Returns the volume vertex (type `Point`) with the lowest coordinates.
    template <typename Point = DefaultPoint>
    decltype(auto) GetRegionLowerCorner() const;
    
    /// Returns the volume vertex (type `Point`) with the highest coordinates.
    template <typename Point = DefaultPoint>
    decltype(auto) GetRegionUpperCorner() const;
    
    /// Returns the number of voxels along each of the three dimensions.
    std::array<unsigned int, 3U> GetSteps() const;

    
    template <typename Vector = DefaultVector>
    Vector GetVoxelSize() const;

    unsigned int GetNVoxels() const;
    
    /// Returns the ID of the voxel containing `p`, or `-1` if none.
    template <typename Point>
    int GetVoxelID(Point const& p) const;
    
    int GetVoxelID(double const*)  const;
    bool IsLegalVoxelID(int) const;

    struct NeiInfo
    {
      NeiInfo() = default;
      NeiInfo(int i, double w) : id(i), weight(w) {}
      int id = -1;
      double weight = 0.0;
    };
    
    
    /**
     * @brief Returns IDs of the eight neighboring voxels around `v`.
     * @param v location within the mapped volume
     * @return an optional collection of eight neighboring voxels
     * 
     * If `v` is not inside the mapped volume, no list is returned (the optional
     * return value evaluates to `false`).
     * Otherwise, each of the eight voxels with the center closest to `v` are
     * returned, each with a weight proportional to the distance of `v` from
     * that center.
     */
    template <typename Point>
    std::optional<std::array<NeiInfo, 8U>> GetNeighboringVoxelIDs(Point const& v) const;

    PhotonVoxel      GetPhotonVoxel(int ID) const;
    std::array<int, 3U> GetVoxelCoords(int ID) const;
    
    /// Returns whether point `p` is inside the region (upper border excluded).
    bool isInside(geo::Point_t const& p) const
      { return isInsideImpl(p); }
    

    bool operator==(const PhotonVoxelDef &rhs) const;
    bool operator!=(const PhotonVoxelDef &rhs) const 
      { return ! ((*this)==rhs); }
    
  private:
    
    int GetVoxelIDImpl(geo::Point_t const& p) const;
    
    std::optional<std::array<NeiInfo, 8U>> GetNeighboringVoxelIDsImpl
      (geo::Point_t const& v) const;
    
    /// Returns the coordinates of the cvoxel containing `p` in step units.
    std::array<double, 3U> GetVoxelStepCoordsUnchecked
      (geo::Point_t const& p) const;
      
    /// Returns whether the specified point is within the volume.
    bool isInsideImpl(geo::Point_t const& point) const
      { return isInsideVolume(point, fLowerCorner, fUpperCorner); }
    
    static bool isInsideVolume(
      geo::Point_t const& point,
      geo::Point_t const& lower, geo::Point_t const& upper
      );
    static bool isInsideRange(double value, double lower, double upper);
    
  }; // class PhotonVoxelDef
  
} // namespace sim


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
//--- sim::PhotonVoxel
//------------------------------------------------------------------------------
template <typename Point /* = DefaultPoint */>
decltype(auto) sim::PhotonVoxel::GetLowerCorner() const
  { return geo::vect::convertTo<Point>(fVoxelMin); }

template <typename Point /* = DefaultPoint */>
decltype(auto) sim::PhotonVoxel::GetUpperCorner() const
  { return geo::vect::convertTo<Point>(fVoxelMax); }
  
template <typename Point /* = DefaultPoint */>
Point sim::PhotonVoxel::GetCenter() const
  { return geo::vect::convertTo<Point>(geo::vect::middlePoint({ fVoxelMin, fVoxelMax })); }


//------------------------------------------------------------------------------
//--- sim::PhotonVoxelDef
//------------------------------------------------------------------------------
template <typename Point /* = DefaultPoint */>
decltype(auto) sim::PhotonVoxelDef::GetRegionLowerCorner() const
  { return geo::vect::convertTo<Point>(fLowerCorner); }
  
template <typename Point /* = DefaultPoint */>
decltype(auto) sim::PhotonVoxelDef::GetRegionUpperCorner() const
  { return geo::vect::convertTo<Point>(fUpperCorner); }

//------------------------------------------------------------------------------
template <typename Vector /* = DefaultVector */>
Vector sim::PhotonVoxelDef::GetVoxelSize() const {
  return {
    (fUpperCorner.X() - fLowerCorner.X()) / fxSteps,
    (fUpperCorner.Y() - fLowerCorner.Y()) / fySteps,
    (fUpperCorner.Z() - fLowerCorner.Z()) / fzSteps
    };
} // sim::PhotonVoxelDef::GetVoxelSize()

//------------------------------------------------------------------------------
template <typename Point>
int sim::PhotonVoxelDef::GetVoxelID(Point const& p) const
  { return GetVoxelIDImpl(geo::vect::toPoint(p)); }

//------------------------------------------------------------------------------
template <typename Point>
std::optional<std::array<sim::PhotonVoxelDef::NeiInfo, 8U>>
sim::PhotonVoxelDef::GetNeighboringVoxelIDs(Point const& v) const
  { return GetNeighboringVoxelIDsImpl(geo::vect::toPoint(v)); }

//------------------------------------------------------------------------------

#endif // LARSIM_SIMULATION_PHOTONVOXELS_H
