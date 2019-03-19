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
  TVector3 PhotonVoxelDef::GetSteps() const
  {
    TVector3 Steps = TVector3(fxSteps, fySteps, fzSteps);
    return Steps;
  }

  //----------------------------------------------------------------------------
  bool PhotonVoxelDef::operator==(const PhotonVoxelDef & right) const
  {
    return ( ( GetRegionUpperCorner() == right.GetRegionUpperCorner() ) &&
             ( GetRegionLowerCorner() == right.GetRegionLowerCorner() ) &&
             ( GetSteps() == right.GetSteps()) );
  }

  //----------------------------------------------------------------------------
  int PhotonVoxelDef::GetNVoxels() const
  {
    return fxSteps * fySteps * fzSteps;
  }

  //----------------------------------------------------------------------------
  int PhotonVoxelDef::GetVoxelID(const TVector3& p) const
  {
    const double xyz[3] = {p.X(), p.Y(), p.Z()};
    return GetVoxelID(xyz);
  }

  //----------------------------------------------------------------------------
  int PhotonVoxelDef::GetVoxelID(double const* Position) const
  {
    // figure out how many steps this point is in the x,y,z directions
    int xStep = int (std::floor((Position[0]-fLowerCorner.X()) / (fUpperCorner.X()-fLowerCorner.X()) * fxSteps ));
    int yStep = int (std::floor((Position[1]-fLowerCorner.Y()) / (fUpperCorner.Y()-fLowerCorner.Y()) * fySteps ));
    int zStep = int (std::floor((Position[2]-fLowerCorner.Z()) / (fUpperCorner.Z()-fLowerCorner.Z()) * fzSteps ));

    // check if point lies within the voxelized region
    if((0 <= xStep) && ((unsigned int)xStep < fxSteps) &&
       (0 <= yStep) && ((unsigned int)yStep < fySteps) &&
       (0 <= zStep) && ((unsigned int)zStep < fzSteps) ){
      // if within bounds, generate the voxel ID
      return (xStep
              + yStep * (fxSteps)
              + zStep * (fxSteps * fySteps));
      }
    else{
      // out of bounds
      return -1;
    }
  }

  //----------------------------------------------------------------------------
  void PhotonVoxelDef::
  GetNeighboringVoxelIDs(const TVector3& v, std::vector<NeiInfo>& ret) const
  {
    ret.clear();
    ret.reserve(8);

    // Position in voxel coordinates including floating point part
    double rStepD[3];
    for(int i = 0; i < 3; ++i){
      // If we're outside the cuboid we have values for, return empty vector,
      // ie failure.
      if(v[i] < GetRegionLowerCorner()[i] || v[i] > GetRegionUpperCorner()[i]) return;// {};
      // Figure out our position wrt to the centres of the voxels
      rStepD[i] = ((v[i]-GetRegionLowerCorner()[i]) / (GetRegionUpperCorner()[i]-GetRegionLowerCorner()[i]) * GetSteps()[i] ) - 0.5;
    }

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
  }

  //----------------------------------------------------------------------------
  TVector3 PhotonVoxelDef::GetVoxelSize() const
  {
    TVector3 TheSize = TVector3((GetRegionUpperCorner()[0]-GetRegionLowerCorner()[0]) / fxSteps,
                                (GetRegionUpperCorner()[1]-GetRegionLowerCorner()[1]) / fySteps,
                                (GetRegionUpperCorner()[2]-GetRegionLowerCorner()[2]) / fzSteps);
    return TheSize;
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
    return (( ID > -1) && (ID<GetNVoxels()));
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
  
} // namespace sim
