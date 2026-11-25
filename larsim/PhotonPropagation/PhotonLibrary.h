////# PhotonLibrary.h header file
////#
////# Ben Jones, MIT, 2012
#ifndef PHOTONLIBRARY_H
#define PHOTONLIBRARY_H

#include "larsim/PhotonPropagation/IPhotonLibrary.h"

#include "larsim/Simulation/PhotonVoxels.h"

#include "TF1.h"
class TTree;

#include <cstdint>
#include <limits> // std::numeric_limits
#include <memory>
#include <optional>
#include <string>
#include <vector>

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

    /**
     * @brief Returns the time parametrization parameters of all the channels for `Voxel`.
     * @return pointer to `NOpChannels()` vectors of parameters, `nullptr` if no timing
     *
     * The vectors live in a single scratch row owned by the library and are refilled when a
     * different voxel is requested: the pointer is valid only until the next call to this
     * method (or to a timing setter).
     */
    const std::vector<float>* GetTimingPars(size_t Voxel) const;

    /**
     * @brief Returns the arrival time distributions of all the channels for `Voxel`.
     * @return pointer to `NOpChannels()` functions, `nullptr` if no timing
     *
     * The library keeps one function per channel (not per voxel), all from the formula of the
     * library, and sets their range and parameters for the requested voxel. The pointer is valid
     * only until the next call to this method (or to a timing setter); after that the functions
     * describe the newly requested voxel.
     * The functions are not constant because `TF1::GetRandom()` is not (see ROOT-9549).
     */
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

    /**
     * @brief Values of one property (e.g. the visibility) of the filled voxels.
     *
     * Each filled voxel has a row of `width()` values (e.g. one per channel), found through
     * `fVoxelTableIdx`. Rows are allocated only on request, zero-initialised, in chunks of
     * `RowsPerChunk` rows; chunks never move, so pointers to rows stay valid while rows are added.
     */
    class PropertyTable {
    public:
      /// Removes all rows and sets the number of values per row (0 disables the storage).
      void reset(size_t width)
      {
        fWidth = width;
        fChunks.clear();
      }

      bool enabled() const { return fWidth > 0; }
      size_t width() const { return fWidth; }

      /// Makes room for row `row`, which must be the one following the last row added.
      void addRow(size_t row)
      {
        if (enabled() && (row % RowsPerChunk == 0))
          fChunks.emplace_back(new float[RowsPerChunk * fWidth]());
      }

      float* row(size_t row) { return fChunks[row / RowsPerChunk].get() + offset(row); }
      float const* row(size_t row) const { return fChunks[row / RowsPerChunk].get() + offset(row); }

      /// Memory allocated for the values, in bytes.
      size_t allocatedBytes() const
      {
        return fChunks.size() * RowsPerChunk * fWidth * sizeof(float);
      }

    private:
      static constexpr size_t RowsPerChunk = 4096;

      size_t fWidth = 0;
      std::vector<std::unique_ptr<float[]>> fChunks;

      size_t offset(size_t row) const { return (row % RowsPerChunk) * fWidth; }
    }; // PropertyTable

    // Only the filled voxels (with at least one value) have a row in the property tables, at the
    // same index `fVoxelTableIdx[v]` in each of the enabled tables below (`NotFoundIdx` if voxel `v`
    // is not filled). Within a row, the channels are stored in sequence (for the timing,
    // `fTimingParNParameters` per channel). The voxels that are not filled read as zero (`fZeroRow`).
    static constexpr std::int32_t NotFoundIdx = -1;

    // int32 to match the upstream voxel ID limit. We may need more in the future...
    std::vector<std::int32_t> fVoxelTableIdx; ///< Index of each voxel in the property tables.
    size_t fNFilledVoxels = 0;                ///< Number of filled voxels (rows in each table).
    PropertyTable fCounts;                    ///< Direct visibilities.
    PropertyTable fReflCounts;                ///< Reflected visibilities.
    PropertyTable fReflT0s;                   ///< Reflected light first arrival times.
    PropertyTable fTimingPars;                ///< Time parametrization parameters.
    std::vector<float> fZeroRow;              ///< Row of zeros, as wide as the widest table.

    std::string fTimingParFormula;    ///< Formula of the time parametrization (from the library).
    size_t fTimingParNParameters = 0; ///< Number of parameters of the time parametrization.
    double fTimingMaxRange = 0.;      ///< Upper end of the range of the time distributions.

    // Scratch rows handed out by `GetTimingTF1s()` and `GetTimingPars()`, one entry per channel,
    // describing the voxel they were last filled for.
    //
    // TF1s are very large in memory compared to its underlying data, so we keep one per channel
    // and set its range and parameters for the requested voxel. We implictly require calls to
    // relinquish the TF1s after use.
    mutable std::vector<TF1> fTimingTF1Row;
    mutable std::optional<size_t> fTimingTF1RowVoxel;
    mutable std::vector<std::vector<float>> fTimingParRow;
    mutable std::optional<size_t> fTimingParRowVoxel;

    size_t fNOpChannels = 0;
    size_t fNVoxels = 0;

    /// Voxel definition loaded from library metadata.
    std::optional<sim::PhotonVoxelDef> fVoxelDef;

    /// ROOT directory where to write data.
    art::TFileDirectory* fDir = nullptr;

    bool isVoxelValidImpl(size_t Voxel) const { return Voxel < fNVoxels; }

    /// Sets up empty storage for `fNVoxels` voxels and `fNOpChannels` channels.
    void resetStorage(bool storeReflected, bool storeReflT0, size_t storeTiming);

    /// Returns the table index of `Voxel`, filling the voxel (with zeros) if it is not yet.
    size_t tableIdxForWriting(size_t Voxel);

    /// Returns the values of `Voxel` in `table` (a row of zeros if the voxel is not filled).
    float const* valuesOf(PropertyTable const& table, size_t Voxel) const
    {
      auto const tableIdx = fVoxelTableIdx[Voxel];
      return (tableIdx == NotFoundIdx) ? fZeroRow.data() : table.row(tableIdx);
    }

    /// Forgets the voxel the timing scratch rows were filled for.
    void invalidateTimingRows() const
    {
      fTimingTF1RowVoxel.reset();
      fTimingParRowVoxel.reset();
    }

    /// Unchecked access to a visibility datum
    float uncheckedAccess(size_t Voxel, size_t OpChannel) const
    {
      return valuesOf(fCounts, Voxel)[OpChannel];
    }

    /// Unchecked access to a visibility datum (adds the row of `Voxel` if needed)
    float& uncheckedAccess(size_t Voxel, size_t OpChannel)
    {
      return fCounts.row(tableIdxForWriting(Voxel))[OpChannel];
    }

    /// Unchecked access to a reflected visibility datum
    float uncheckedAccessRefl(size_t Voxel, size_t OpChannel) const
    {
      return valuesOf(fReflCounts, Voxel)[OpChannel];
    }

    /// Unchecked access to a reflected visibility datum (adds the row of `Voxel` if needed)
    float& uncheckedAccessRefl(size_t Voxel, size_t OpChannel)
    {
      return fReflCounts.row(tableIdxForWriting(Voxel))[OpChannel];
    }

    /// Unchecked access to a reflected T0 visibility datum
    float uncheckedAccessReflT(size_t Voxel, size_t OpChannel) const
    {
      return valuesOf(fReflT0s, Voxel)[OpChannel];
    }

    /// Unchecked access to a reflected T0 visibility datum (adds the row of `Voxel` if needed)
    float& uncheckedAccessReflT(size_t Voxel, size_t OpChannel)
    {
      return fReflT0s.row(tableIdxForWriting(Voxel))[OpChannel];
    }

    /// Unchecked access to a parameter the time distribution
    float uncheckedAccessTimingPar(size_t Voxel, size_t OpChannel, size_t parnum) const
    {
      return valuesOf(fTimingPars, Voxel)[OpChannel * fTimingParNParameters + parnum];
    }

    /// Unchecked access to a parameter of the time distribution (adds the row of `Voxel` if
    /// needed)
    float& uncheckedAccessTimingPar(size_t Voxel, size_t OpChannel, size_t parnum)
    {
      invalidateTimingRows();
      return fTimingPars.row(tableIdxForWriting(Voxel))[OpChannel * fTimingParNParameters + parnum];
    }

    /// Creates one timing function per channel from the formula of the library.
    void makeTimingTF1Row();

    /// Fills the timing scratch rows with the parameters of `Voxel`, unless they already are.
    void fillTimingTF1Row(size_t Voxel) const;
    void fillTimingParRow(size_t Voxel) const;

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
