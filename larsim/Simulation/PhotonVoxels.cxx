/**
 * @file  larsim/Simulation/PhotonVoxels.cxx
 * @brief Definitions of voxel data structures: implementation.
 * @see   larsim/Simulation/PhotonVoxels.h
 */

// library header
#include "larsim/Simulation/PhotonVoxels.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcorealg/Geometry/geo_vectors_utils_TVector.h"

// C++ standard libraries
#include <vector>
#include <string>
#include <algorithm> // std::min(), std::max()
#include <stdexcept> // std::runtime_error
#include <cmath> // std::abs(), std::floor()


namespace sim {


  //----------------------------------------------------------------------------
  // PhotonVoxelDef class
  //----------------------------------------------------------------------------
  PhotonVoxelDef::PhotonVoxelDef(double xMin, 
                                 double xMax, 
                                 int xN, 
                                 double yMin, 
                                 double yMax, 
                                 int yN, 
                                 double zMin, 
                                 double zMax, 
                                 int zN)
    : fLowerCorner(xMin, yMin, zMin)
    , fUpperCorner(xMax, yMax, zMax)
    , fxSteps(xN)
    , fySteps(yN)
    , fzSteps(zN)
    {}

  //----------------------------------------------------------------------------
  std::array<unsigned int, 3U> PhotonVoxelDef::GetSteps() const
  {
    return { fxSteps, fySteps, fzSteps };
  }

  //----------------------------------------------------------------------------
  bool PhotonVoxelDef::operator==(const PhotonVoxelDef & right) const
  {
    return ( ( GetRegionUpperCorner() == right.GetRegionUpperCorner() ) &&
             ( GetRegionLowerCorner() == right.GetRegionLowerCorner() ) &&
             ( GetSteps() == right.GetSteps()) );
  }

  //----------------------------------------------------------------------------
  unsigned int PhotonVoxelDef::GetNVoxels() const
  {
    return fxSteps * fySteps * fzSteps;
  }

  //----------------------------------------------------------------------------
  int PhotonVoxelDef::GetVoxelID(double const* Position) const {
    return GetVoxelIDImpl(geo::vect::makeFromCoords<geo::Point_t>(Position));
  }

  //----------------------------------------------------------------------------
  std::vector<PhotonVoxelDef::NeiInfo> PhotonVoxelDef::GetNeighboringVoxelIDsImpl
    (geo::Point_t const& v) const
  {
    if (!isInside(v)) return {};
    
    std::vector<NeiInfo> ret;
    ret.reserve(8);

    // Position in voxel coordinates including floating point part
    auto const rStepD = GetVoxelStepCoordsUnchecked(v);
    
    // The neighbours are the 8 corners of a cube around this point
    for(int dx = 0; dx <= 1; ++dx){
      for(int dy = 0; dy <= 1; ++dy){
        for(int dz = 0; dz <= 1; ++dz){
          // The full 3D step
          const int dr[3] = {dx, dy, dz};

          // The integer-only position of the current corner
          int rStepI[3];
          for(int d = 0; d < 3; ++d){
            // Round down to get the "lower left" corner
            rStepI[d] = int(rStepD[d]);
            // Ensure we'll stay in-bounds
            rStepI[d] = std::max(0, rStepI[d]);
            rStepI[d] = std::min(rStepI[d], int(GetSteps()[d])-2);
            // Adjust to the corner we're actually considering
            rStepI[d] += dr[d];
          }

          double w = 1;
          for(int d = 0; d < 3; ++d){
            // These expressions will interpolate when between the 8 corners,
            // and extrapolate in the half-voxel space around the edges.
            if(dr[d] == 0)
              w *= 1+rStepI[d]-rStepD[d];
            else
              w *= 1-rStepI[d]+rStepD[d];
          }

          const int id = (rStepI[0] +
                          rStepI[1] * (fxSteps) +
                          rStepI[2] * (fxSteps * fySteps));

          ret.emplace_back(id, w);
        }
      }
    }

    // Sanity check the weights sum to 1
    double wSum = 0;
    for(const NeiInfo& n: ret) wSum += n.weight;
    if(std::abs(wSum-1) > 1e-3){
      std::string msg
        = "PhotonVoxelDef::GetNeighboringVoxelIDs():"
          " Weights sum to " + std::to_string(wSum) + " (should be 1)."
          " Weights are:";
      for(const NeiInfo& n: ret) {
        msg += ' ';
        msg += std::to_string(n.weight);
      }
      throw std::runtime_error(msg);
    }
    return ret;
  }

  //----------------------------------------------------------------------------
  PhotonVoxel PhotonVoxelDef::GetPhotonVoxel(int ID) const
  {
    //    float TempID = (float) ID;

    // Decompose ID into steps in each direction
    int xStep =  ID % fxSteps ;
    int yStep =  ((ID - xStep ) / fxSteps) % fySteps ;
    int zStep =  ((ID - xStep - (yStep * fxSteps)) / (fySteps * fxSteps)) % fzSteps ;


    TVector3 VoxelSize = GetVoxelSize();

    double xMin = VoxelSize[0] * (xStep)   + fLowerCorner.X();
    double xMax = VoxelSize[0] * (xStep+1) + fLowerCorner.X();
    double yMin = VoxelSize[1] * (yStep)   + fLowerCorner.Y();
    double yMax = VoxelSize[1] * (yStep+1) + fLowerCorner.Y();
    double zMin = VoxelSize[2] * (zStep)   + fLowerCorner.Z();
    double zMax = VoxelSize[2] * (zStep+1) + fLowerCorner.Z();


   
    return PhotonVoxel(xMin, xMax, yMin, yMax, zMin, zMax);
  }

  //----------------------------------------------------------------------------
  bool PhotonVoxelDef::IsLegalVoxelID(int ID) const
  {
    return (( ID >= 0) && (static_cast<unsigned int>(ID) < GetNVoxels()));
  }

  std::vector<int> PhotonVoxelDef::GetVoxelCoords(int ID) const
  {
    std::vector<int> ReturnVector;
    ReturnVector.resize(3);
    ReturnVector.at(0) =  ID % fxSteps ;
    ReturnVector.at(1) =  ((ID - ReturnVector.at(0) ) / fxSteps) % fySteps ;
    ReturnVector.at(2) =  ((ID - ReturnVector.at(0) - (ReturnVector.at(1) * fxSteps)) / (fySteps * fxSteps)) % fzSteps ;
    return ReturnVector;
    
  }
  
  //----------------------------------------------------------------------------
  std::array<double, 3U> PhotonVoxelDef::GetVoxelStepCoordsUnchecked
    (geo::Point_t const& p) const
  {
    
    auto const span = fUpperCorner - fLowerCorner;
    auto const relPos = p - fLowerCorner;
    
    return {
      (relPos.X() / span.X()) * fxSteps,
      (relPos.Y() / span.Y()) * fySteps,
      (relPos.Z() / span.Z()) * fzSteps
      };
  } // PhotonVoxelDef::GetVoxelStepCoordsUnchecked()
  
  //----------------------------------------------------------------------------
  int PhotonVoxelDef::GetVoxelIDImpl(geo::Point_t const& p) const {
    if (!isInside(p)) return -1;
    
    auto const stepCoords = GetVoxelStepCoordsUnchecked(p);
    
    // figure out how many steps this point is in the x,y,z directions;
    // `p` is guaranteed to be in the mapped volume by the previous check
    int xStep = static_cast<int>(stepCoords[0]);
    int yStep = static_cast<int>(stepCoords[1]);
    int zStep = static_cast<int>(stepCoords[2]);

    // if within bounds, generate the voxel ID
    return (xStep
            + yStep * (fxSteps)
            + zStep * (fxSteps * fySteps));
  }
  
  //----------------------------------------------------------------------------
  bool PhotonVoxelDef::isInsideVolume(
    geo::Point_t const& point,
    geo::Point_t const& lower, geo::Point_t const& upper
    )
  {
    return
         isInsideRange(point.X(), lower.X(), upper.X())
      && isInsideRange(point.Y(), lower.Y(), upper.Y())
      && isInsideRange(point.Z(), lower.Z(), upper.Z())
      ;
  }
  
  bool PhotonVoxelDef::isInsideRange(double value, double lower, double upper) {
    
    return (value >= lower) && (value < upper);
    
  } // PhotonVoxelDef::isInsideRange()
  
  
  //----------------------------------------------------------------------------
  
  
} // namespace sim
