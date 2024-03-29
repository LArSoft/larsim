////# PhotonLibrary.h header file
////#
////# Ben Jones, MIT, 2012
#ifndef PHOTONLIBRARY_H
#define PHOTONLIBRARY_H

#include "larsim/PhotonPropagation/IPhotonLibrary.h"

#include "larsim/Simulation/PhotonVoxels.h"

#include "TF1.h"
class TTree;

#include "lardataobj/Utilities/LazyVector.h"

#include <limits> // std::numeric_limits
#include <optional>

namespace art {
  class TFileDirectory;
}

namespace phot {

  class PhotonLibrary : public IPhotonLibrary {
  public:
    /// If no valid `pDir` is provided, storage features will not be supported.
    PhotonLibrary(art::TFileDirectory* pDir = nullptr);

    TTree* ProduceTTree() const;

    virtual float GetCount(size_t Voxel, size_t OpChannel) const override;
    void SetCount(size_t Voxel, size_t OpChannel, float Count);

    float GetTimingPar(size_t Voxel, size_t OpChannel, size_t parnum) const;
    void SetTimingPar(size_t Voxel, size_t OpChannel, float Count, size_t parnum);

    //    TF1& GetTimingTF1(size_t Voxel, size_t OpChannel) const;
    void SetTimingTF1(size_t Voxel, size_t OpChannel, TF1 func);

    virtual float GetReflCount(size_t Voxel, size_t OpChannel) const override;
    void SetReflCount(size_t Voxel, size_t OpChannel, float Count);

    virtual float GetReflT0(size_t Voxel, size_t OpChannel) const override;
    void SetReflT0(size_t Voxel, size_t OpChannel, float reflT0);

    /// Returns a pointer to NOpChannels() visibility values, one per channel
    virtual float const* GetCounts(size_t Voxel) const override;
    const std::vector<float>* GetTimingPars(size_t Voxel) const;
    TF1* GetTimingTF1s(size_t Voxel) const;

    virtual float const* GetReflCounts(size_t Voxel) const override;
    virtual float const* GetReflT0s(size_t Voxel) const override;

    ///Returns whether the current library deals with time propagation distributions.
    bool hasTiming() const { return fHasTiming; }

    /// Returns whether the current library deals with reflected light count.
    virtual bool hasReflected() const override { return fHasReflected; }

    /// Returns whether the current library deals with reflected light timing.
    virtual bool hasReflectedT0() const override { return fHasReflectedT0; }

    void StoreLibraryToFile(std::string LibraryFile,
                            bool storeReflected = false,
                            bool storeReflT0 = false,
                            size_t storeTiming = 0) const;
    void LoadLibraryFromFile(std::string LibraryFile,
                             size_t NVoxels,
                             bool storeReflected = false,
                             bool storeReflT0 = false,
                             size_t storeTiming = 0,
                             int maxrange = 200);
    void CreateEmptyLibrary(size_t NVoxels,
                            size_t NChannels,
                            bool storeReflected = false,
                            bool storeReflT0 = false,
                            size_t storeTiming = 0);

    // --- BEGIN --- Metadata: voxel information -------------------------------
    /// @name Metadata: voxel information
    /// @{

    /// Returns whether voxel metadata is available.
    bool hasVoxelDef() const { return fVoxelDef.has_value(); }

    /// Returns the current voxel metadata (undefined behaviour if none).
    /// @see `hasVoxelDef()`
    sim::PhotonVoxelDef const& GetVoxelDef() const
    {
      assert(fVoxelDef);
      return *fVoxelDef;
    }

    /// Copies the specified voxel definition into our own
    /// (overwrites the existing metadata if any).
    void SetVoxelDef(sim::PhotonVoxelDef const& voxelDef) { fVoxelDef = voxelDef; }

    /// @}
    // --- END --- Metadata: voxel information ---------------------------------

  private:
    virtual int NOpChannels() const override { return fNOpChannels; }
    virtual int NVoxels() const override { return fNVoxels; }

    virtual bool isVoxelValid(size_t Voxel) const override { return isVoxelValidImpl(Voxel); }

    bool fHasReflected = false; ///< Whether the current library deals with reflected light counts.
    bool fHasReflectedT0 =
      false; ///< Whether the current library deals with reflected light timing.

    size_t fHasTiming =
      0; ///< Whether the current library deals with time propagation distribution.

    // fLookupTable[unchecked_index(Voxel, OpChannel)] = Count
    // for each voxel, all NChannels() channels are stored in sequence
    util::LazyVector<float> fLookupTable;
    util::LazyVector<float> fReflLookupTable;
    util::LazyVector<float> fReflTLookupTable;
    util::LazyVector<std::vector<float>> fTimingParLookupTable;
    util::LazyVector<TF1> fTimingParTF1LookupTable;
    std::string fTimingParFormula;
    size_t fTimingParNParameters;

    size_t fNOpChannels;
    size_t fNVoxels;

    /// Voxel definition loaded from library metadata.
    std::optional<sim::PhotonVoxelDef> fVoxelDef;

    /// ROOT directory where to write data.
    art::TFileDirectory* fDir = nullptr;

    bool isVoxelValidImpl(size_t Voxel) const { return Voxel < fNVoxels; }

    /// Returns the index of visibility of specified voxel and cell
    size_t uncheckedIndex(size_t Voxel, size_t OpChannel) const
    {
      return Voxel * fNOpChannels + OpChannel;
    }

    /// Unchecked access to a visibility datum
    float uncheckedAccess(size_t Voxel, size_t OpChannel) const
    {
      return fLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a visibility datum
    float& uncheckedAccess(size_t Voxel, size_t OpChannel)
    {
      return fLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a reflected visibility datum
    float uncheckedAccessRefl(size_t Voxel, size_t OpChannel) const
    {
      return fReflLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a reflected visibility datum
    float& uncheckedAccessRefl(size_t Voxel, size_t OpChannel)
    {
      return fReflLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a reflected T0 visibility datum
    float uncheckedAccessReflT(size_t Voxel, size_t OpChannel) const
    {
      return fReflTLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a reflected T0 visibility datum
    float& uncheckedAccessReflT(size_t Voxel, size_t OpChannel)
    {
      return fReflTLookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a parameter the time distribution
    float uncheckedAccessTimingPar(size_t Voxel, size_t OpChannel, size_t parnum) const
    {
      return fTimingParLookupTable[uncheckedIndex(Voxel, OpChannel)][parnum];
    }

    /// Unchecked access to a parameter of the time distribution
    float& uncheckedAccessTimingPar(size_t Voxel, size_t OpChannel, size_t parnum)
    {
      return fTimingParLookupTable[uncheckedIndex(Voxel, OpChannel)][parnum];
    }

    /// Unchecked access to a parameter of the time distribution
    TF1& uncheckedAccessTimingTF1(size_t Voxel, size_t OpChannel)
    {
      return fTimingParTF1LookupTable[uncheckedIndex(Voxel, OpChannel)];
    }

    /// Unchecked access to a parameter of the time distribution
    const TF1& uncheckedAccessTimingTF1(size_t Voxel, size_t OpChannel) const
    {
      // note that this will produce a segmentation fault if the formula is not there
      return *(fTimingParTF1LookupTable.data_address(uncheckedIndex(Voxel, OpChannel)));
    }

    /// Reads the metadata from specified ROOT directory and sets it as current.
    void LoadMetadata(TDirectory& srcDir);

    /// Writes the current metadata (if any) into the ROOT output file.
    void StoreMetadata() const;

    /// Name of the optical channel number in the input tree
    static std::string const OpChannelBranchName;

    /// Returns the number of optical channels in the specified tree
    static size_t ExtractNOpChannels(TTree* tree);

    /// Converts size_t into integer
    static int size_t2int(size_t val)
    {
      constexpr size_t int_max_as_size_t{std::numeric_limits<int>::max()};
      return val <= int_max_as_size_t ? static_cast<int>(val) : -1;
    }
  };

}

#endif
