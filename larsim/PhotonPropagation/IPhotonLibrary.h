////# IPhotonLibrary.h header file
////#
////# Chris Backhouse, UCL, 2017
#ifndef IPHOTONLIBRARY_H
#define IPHOTONLIBRARY_H

#include <vector>
#include <cstddef> // size_t

class TF1; // forward declaration

namespace phot
{
  /// Interface shared by all PhotonLibrary-like classes
  class IPhotonLibrary
  {
  public:
    /// Type for visibility count per optical channel.
    using Counts_t = const float*;
    
    /// Type for time of arrival per optical channel.
    using T0s_t = const float*;
    
    /// Type for function parameters
    /// (which is not part of this interface yet).
    using Params_t = std::vector<float> const*;
    
    /// Type for parametrization function
    /// (which is not part of this interface yet).
    using Functions_t = TF1*;
    
    
    virtual ~IPhotonLibrary() = default;

    virtual float GetCount(size_t Voxel, size_t OpChannel) const = 0;
    virtual float GetReflCount(size_t Voxel, size_t OpChannel) const = 0;
    virtual float GetReflT0(size_t Voxel, size_t OpChannel) const = 0;

    /// Returns a pointer to NOpChannels() visibility values, one per channel
    virtual Counts_t GetCounts(size_t Voxel) const = 0;
    virtual Counts_t GetReflCounts(size_t Voxel) const = 0;
    virtual T0s_t    GetReflT0s(size_t Voxel) const = 0;

    /// Returns whether the current library deals with reflected light count.
    virtual bool hasReflected() const = 0;

    /// Returns whether the current library deals with reflected light timing.
    virtual bool hasReflectedT0() const = 0;

    virtual int NOpChannels() const = 0;
    virtual int NVoxels() const = 0;
    
    virtual bool isVoxelValid(size_t Voxel) const
      { return Voxel < (std::size_t) NVoxels(); }

    /// Returns the number of elements in the library
    size_t LibrarySize() const { return NVoxels() * NOpChannels(); }
  };
} // namespace

#endif
