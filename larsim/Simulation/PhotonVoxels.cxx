#include "larsim/Simulation/PhotonVoxels.h"

#include <iostream>

namespace sim {


  // PhotonVoxel class
  //----------------------------------------------------------------------------
  PhotonVoxel::PhotonVoxel(double xMin, 
			   double xMax, 
			   double yMin, 
			   double yMax, 
			   double zMin, 
			   double zMax, 
			   int N)
  {
    xVoxelMin = xMin;
    xVoxelMax = xMax;
    yVoxelMin = yMin;
    yVoxelMax = yMax;
    zVoxelMin = zMin;
    zVoxelMax = zMax;
    NPhotons = N;
  }

  //----------------------------------------------------------------------------
  PhotonVoxel::PhotonVoxel()
  {
  }
  
  //----------------------------------------------------------------------------
  TVector3 PhotonVoxel::GetLowerCorner() const
  {
    TVector3 LowerCorner = TVector3(xVoxelMin, yVoxelMin, zVoxelMin);
    return LowerCorner;
  }


  //----------------------------------------------------------------------------
  TVector3 PhotonVoxel::GetUpperCorner() const
  {
    TVector3 UpperCorner = TVector3(xVoxelMax, yVoxelMax, zVoxelMax);
    return UpperCorner;
  }


  //----------------------------------------------------------------------------
  TVector3 PhotonVoxel::GetCenter() const
  {
    TVector3 Center = TVector3((xVoxelMin+xVoxelMax)/2.0, (yVoxelMin+yVoxelMax)/2.0, (zVoxelMin+zVoxelMax)/2.0);
    return Center;
  }


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
  {
    fxSteps = xN;
    fySteps = yN;
    fzSteps = zN;
    
    
    fLowerCorner = TVector3(xMin,yMin,zMin);
    fUpperCorner = TVector3(xMax,yMax,zMax);
  }

  //----------------------------------------------------------------------------
  PhotonVoxelDef::PhotonVoxelDef()
  {
  }

  //----------------------------------------------------------------------------
  TVector3 PhotonVoxelDef::GetRegionLowerCorner() const
  {
    return fLowerCorner;
  }

  //----------------------------------------------------------------------------
  TVector3 PhotonVoxelDef::GetRegionUpperCorner() const
  {
    return fUpperCorner;
  }

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
    int xStep = int ((Position[0]-fLowerCorner[0]) / (fUpperCorner[0]-fLowerCorner[0]) * fxSteps );
    int yStep = int ((Position[1]-fLowerCorner[1]) / (fUpperCorner[1]-fLowerCorner[1]) * fySteps );
    int zStep = int ((Position[2]-fLowerCorner[2]) / (fUpperCorner[2]-fLowerCorner[2]) * fzSteps );

    // check if point lies within the voxelized region
    if((0 <= xStep) && (xStep < fxSteps) &&
       (0 <= yStep) && (yStep < fySteps) &&
       (0 <= zStep) && (zStep < fzSteps) ){
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
  std::vector<PhotonVoxelDef::NeiInfo> PhotonVoxelDef::
  GetNeighboringVoxelIDs(const TVector3& v) const
  {
    std::vector<NeiInfo> ret;

    // Deliberately keep the floating part
    double rStepD[3];
    for(int i = 0; i < 3; ++i){
      if(v[i] < fLowerCorner[i] || v[i] > fUpperCorner[i]) return {};
      rStepD[i] = ((v[i]-fLowerCorner[i]) / (fUpperCorner[i]-fLowerCorner[i]) * GetSteps()[i] );
    }

    for(int dx = 0; dx <= 1; ++dx){
      for(int dy = 0; dy <= 1; ++dy){
        for(int dz = 0; dz <= 1; ++dz){
          const int dr[3] = {dx, dy, dz};

          // TODO comment about how the maths works
          int rStepI[3];
          for(int d = 0; d < 3; ++d){
            rStepI[d] = int(rStepD[d]-.5);
            rStepI[d] = std::max(0, rStepI[d]);
            rStepI[d] = std::min(rStepI[d], int(GetSteps()[d])-2);
            rStepI[d] += dr[d];
          }

          double w = 1;
          for(int d = 0; d < 3; ++d) w *= 1-fabs(rStepD[d]-rStepI[d]);

          const int id = (rStepI[0] +
                          rStepI[1] * (fxSteps) +
                          rStepI[2] * (fxSteps * fySteps));

          ret.emplace_back(id, w);
        }
      }
    }

    // Check the weights sum to 1
    double wSum = 0;
    for(const NeiInfo& n: ret) wSum += n.weight;
    if(fabs(wSum-1) > 1e-3){
      std::cout << "PhotonVoxelDef::GetNeighboringVoxelIDs(): "
                << "But weights, sum to " << wSum << " (should be 1)."
                << "Weights are:";
      for(const NeiInfo& n: ret) std::cout << " " << n.weight;
      std::cout << " Aborting." << std::endl;
      abort();
    }

    return ret;
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

    double xMin = VoxelSize[0] * (xStep)   + fLowerCorner[0];
    double xMax = VoxelSize[0] * (xStep+1) + fLowerCorner[0];
    double yMin = VoxelSize[1] * (yStep)   + fLowerCorner[1];
    double yMax = VoxelSize[1] * (yStep+1) + fLowerCorner[1];
    double zMin = VoxelSize[2] * (zStep)   + fLowerCorner[2];
    double zMax = VoxelSize[2] * (zStep+1) + fLowerCorner[2];


   
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
}
