/**
 * @file  larsim/Simulation/PhotonVoxels.h
 * @brief Definitions of voxel data structures.
 */

#ifndef LARSIM_SIMULATION_PHOTONVOXELS_H
#define LARSIM_SIMULATION_PHOTONVOXELS_H

#include "TVector3.h"

namespace sim {


  class PhotonVoxel{
  public:
    PhotonVoxel() = default;
    PhotonVoxel(double xMin, 
                double xMax, 
                double yMin, 
                double yMax, 
                double zMin, 
                double zMax, 
                int N = 0) ;

  private:
    double xVoxelMin = 0.0;
    double xVoxelMax = 0.0;
    double yVoxelMin = 0.0;
    double yVoxelMax = 0.0;
    double zVoxelMin = 0.0;
    double zVoxelMax = 0.0;

    int NPhotons;

  public:

    TVector3 GetLowerCorner() const;
    TVector3 GetUpperCorner() const;
    TVector3 GetCenter()      const;

  }; // class PhotonVoxel


  class PhotonVoxelDef
  {
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
    
  private:
    TVector3 fLowerCorner;
    TVector3 fUpperCorner;
    int      fxSteps;
    int      fySteps;
    int      fzSteps;

  public:

    TVector3 GetRegionUpperCorner() const;
    TVector3 GetRegionLowerCorner() const;
    TVector3 GetSteps() const;


    TVector3 GetVoxelSize()const;

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

#endif // LARSIM_SIMULATION_PHOTONVOXELS_H
