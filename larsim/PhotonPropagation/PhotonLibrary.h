////# PhotonLibrary.h header file
////#
////# Ben Jones, MIT, 2012
#ifndef PHOTONLIBRARY_H
#define PHOTONLIBRARY_H

#include "larsim/PhotonPropagation/IPhotonLibrary.h"

#include "TTree.h"
#include "larsim/Simulation/PhotonVoxels.h"


namespace phot{
  
  class PhotonLibrary : public IPhotonLibrary
  {
  public:
    PhotonLibrary();
    ~PhotonLibrary();

    TTree * ProduceTTree();
    

    virtual float GetCount(size_t Voxel, size_t OpChannel) const override;
    void SetCount(size_t Voxel, size_t OpChannel, float Count);

    virtual float GetReflCount(size_t Voxel, size_t OpChannel) const override;
    void SetReflCount(size_t Voxel, size_t OpChannel, float Count);

    virtual float GetReflT0(size_t Voxel, size_t OpChannel) const override;
    void SetReflT0(size_t Voxel, size_t OpChannel, float reflT0);
    
    /// Returns a pointer to NOpChannels() visibility values, one per channel
    virtual float const* GetCounts(size_t Voxel) const override;
    virtual float const* GetReflCounts(size_t Voxel) const override;
    virtual float const* GetReflT0s(size_t Voxel) const override;

    /// Returns whether the current library deals with reflected light count.
    virtual bool hasReflected() const override { return fHasReflected; }
    
    /// Returns whether the current library deals with reflected light timing.
    virtual bool hasReflectedT0() const override { return fHasReflectedT0; }

    
    void StoreLibraryToFile(std::string LibraryFile, bool storeReflected=false, bool storeReflT0=false);
    void LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels, bool storeReflected=false, bool storeReflT0=false);
    void CreateEmptyLibrary(size_t NVoxels, size_t NChannels, bool storeReflected=false, bool storeReflT0=false);
    

    virtual int NOpChannels() const override { return fNOpChannels; }
    virtual int NVoxels() const override { return fNVoxels; }

  private:
    
    bool fHasReflected   = false; ///< Whether the current library deals with reflected light counts.
    bool fHasReflectedT0 = false; ///< Whether the current library deals with reflected light timing.
    
    // fLookupTable[unchecked_index(Voxel, OpChannel)] = Count
    // for each voxel, all NChannels() channels are stored in sequence
    std::vector<float> fLookupTable;
    std::vector<float> fReflLookupTable;
    std::vector<float> fReflTLookupTable;
    size_t fNOpChannels;
    size_t fNVoxels;
    
    /// Returns the index of visibility of specified voxel and cell
    size_t uncheckedIndex(size_t Voxel, size_t OpChannel) const
      { return Voxel * fNOpChannels + OpChannel; }
    
    /// Unchecked access to a visibility datum
    float const& uncheckedAccess (size_t Voxel, size_t OpChannel) const
      { return fLookupTable[uncheckedIndex(Voxel, OpChannel)]; }
    
    /// Unchecked access to a visibility datum
    float& uncheckedAccess(size_t Voxel, size_t OpChannel)
      { return fLookupTable[uncheckedIndex(Voxel, OpChannel)]; }

    /// Unchecked access to a reflected visibility datum
    float const& uncheckedAccessRefl (size_t Voxel, size_t OpChannel) const
    { return fReflLookupTable[uncheckedIndex(Voxel, OpChannel)]; }

    /// Unchecked access to a reflected visibility datum                                                                                         
    float& uncheckedAccessRefl(size_t Voxel, size_t OpChannel)
    { return fReflLookupTable[uncheckedIndex(Voxel, OpChannel)]; }
    
    /// Unchecked access to a reflected T0 visibility datum
    float const& uncheckedAccessReflT (size_t Voxel, size_t OpChannel) const
    { return fReflTLookupTable[uncheckedIndex(Voxel, OpChannel)]; }

    /// Unchecked access to a reflected T0 visibility datum                                                                                        
    float& uncheckedAccessReflT(size_t Voxel, size_t OpChannel)
    { return fReflTLookupTable[uncheckedIndex(Voxel, OpChannel)]; }

    /// Name of the optical channel number in the input tree
    static std::string const OpChannelBranchName;
    
    /// Returns the number of optical channels in the specified tree
    static size_t ExtractNOpChannels(TTree* tree);
    
  };

}

#endif
