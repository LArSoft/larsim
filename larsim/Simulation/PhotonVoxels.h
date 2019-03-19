/**
 * @file  larsim/Simulation/PhotonVoxels.h
 * @brief Definitions of voxel data structures.
 */

#ifndef LARSIM_SIMULATION_PHOTONVOXELS_H
#define LARSIM_SIMULATION_PHOTONVOXELS_H

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "TVector3.h"

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
    using DefaultPoint = TVector3; // legacy; it should really be `geo::Point_t`
    
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
    using DefaultPoint = TVector3; // legacy; it should really be `geo::Point_t`
    using DefaultVector = TVector3; // legacy; it should really be `geo::Vector_t`
    
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
    
    TVector3 GetSteps() const;

    
    template <typename Vector = DefaultVector>
    Vector GetVoxelSize() const;

    int GetNVoxels() const;

    int GetVoxelID(const TVector3&) const;
    int GetVoxelID(double const*)  const;
    bool IsLegalVoxelID(int) const;

    struct NeiInfo
    {
      NeiInfo(int i, double w) : id(i), weight(w) {}
      int id;
      double weight;
    };

    // Out-param allows less allocation if caller re-uses a buffer
    void GetNeighboringVoxelIDs(const TVector3& v,
                                std::vector<NeiInfo>& ret) const;

    PhotonVoxel      GetPhotonVoxel(int ID) const;
    std::vector<int> GetVoxelCoords(int ID) const;
    PhotonVoxel      GetContainingVoxel(TVector3) const;

    bool operator==(const PhotonVoxelDef &rhs) const;
    bool operator!=(const PhotonVoxelDef &rhs) const 
      { return ! ((*this)==rhs); }

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

#endif // LARSIM_SIMULATION_PHOTONVOXELS_H
