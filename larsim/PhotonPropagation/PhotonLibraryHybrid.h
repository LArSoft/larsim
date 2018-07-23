////# PhotonLibraryHybrid.h header file
////#
////# Chris Backhouse, UCL, 2017
#ifndef PHOTONLIBRARYHYBRID_H
#define PHOTONLIBRARYHYBRID_H

#include "larsim/PhotonPropagation/IPhotonLibrary.h"

#include <cmath>
#include <string>
#include <vector>

namespace sim{class PhotonVoxelDef;}

namespace phot
{
  /// TODO doc
  class PhotonLibraryHybrid: public IPhotonLibrary
  {
  public:
    PhotonLibraryHybrid(const std::string& fname,
                        const sim::PhotonVoxelDef& voxdef);
    virtual ~PhotonLibraryHybrid();

    virtual float GetCount(size_t Voxel, size_t OpChannel) const override;

    // This one is unimplemented
    virtual const float* GetCounts(size_t Voxel) const override;

    /// Don't implement reflected light
    virtual bool hasReflected() const override {return false;}
    virtual const float* GetReflCounts(size_t Voxel) const override {return 0;}
    virtual float GetReflCount(size_t Voxel, size_t OpChannel) const override {return 0;}

    /// Don't implement reflected light timing
    virtual bool hasReflectedT0() const override {return false;}
    virtual const float* GetReflT0s(size_t Voxel) const override {return 0;}
    virtual float GetReflT0(size_t Voxel, size_t OpChannel) const override {return 0;}

    virtual int NOpChannels() const override {return fRecords.size();}
    virtual int NVoxels() const override;

  protected:
    const sim::PhotonVoxelDef& fVoxDef;

    struct FitFunc
    {
      FitFunc() {}
      FitFunc(float n, float d) : norm(n), decay(d) {}
      double Eval(double x) const {return exp(norm+decay*x)/(x*x);}

      float norm, decay;
    };

    struct Exception
    {
      Exception(size_t _vox, float _vis) : vox(_vox), vis(_vis) {}

      bool operator<(const Exception& e) const {return vox < e.vox;}
      bool operator<(size_t v) const {return vox < v;}

      size_t vox;
      float vis;
    };

    struct OpDetRecord
    {
      FitFunc fit;
      // std::map is particularly space-inefficient, and sorted vector can also
      // be searched in log(N) time.
      std::vector<Exception> exceptions;
    };

    std::vector<OpDetRecord> fRecords;
  };
} // namespace

#endif
