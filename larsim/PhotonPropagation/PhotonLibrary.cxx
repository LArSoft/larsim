#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/Exception.h"

#include "cetlib_except/exception.h"
#include "larsim/PhotonPropagation/PhotonLibrary.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RtypesCore.h"
#include "TBranch.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TKey.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include "TTree.h"
#include "TVector.h"

#include <algorithm>
#include <cassert>
#include <limits>

namespace {

  /// Reads a scalar metadata value with the specified `key` from `srcDir`,
  /// supporting both the current format (`TParameter<T>`) and the legacy one
  /// (RooFit's `RooInt` and `RooDouble`). Missing keys end up in `missingKeys`.
  template <typename T>
  struct MetadataReader {
    TDirectory& srcDir;
    std::vector<std::string>& missingKeys;

    std::optional<T> operator()(std::string const& key) const
    {
      if (auto const* param = srcDir.Get<TParameter<T>>(key.c_str())) return param->GetVal();
      if (auto value = readLegacy(key)) return value;
      missingKeys.push_back(key);
      return std::nullopt;
    }

  private:
    /// Reads a value stored as RooFit `RooInt` or `RooDouble` without relying
    /// on RooFit itself: `RooInt` does not exist any more (since ROOT 6.30),
    /// and its objects can be read back only through the "emulated" class that
    /// ROOT builds from the streamer information stored in the file itself.
    std::optional<T> readLegacy(std::string const& key) const
    {
      TKey* dirKey = srcDir.GetKey(key.c_str());
      if (!dirKey) return std::nullopt;
      std::string const className = dirKey->GetClassName();
      bool const isInt = (className == "RooInt");
      if (!isInt && className != "RooDouble") return std::nullopt;
      TClass* actualClass = TClass::GetClass(className.c_str());
      if (!actualClass) return std::nullopt;
      // ReadObjectAny() (unlike ReadObj()) streams emulated classes correctly
      void* obj = dirKey->ReadObjectAny(nullptr);
      if (!obj) return std::nullopt;
      std::optional<T> value;
      if (Longptr_t const offset = actualClass->GetDataMemberOffset("_value"); offset > 0) {
        auto const* valuePtr = reinterpret_cast<char const*>(obj) + offset;
        value = isInt ? static_cast<T>(*reinterpret_cast<Int_t const*>(valuePtr)) :
                        static_cast<T>(*reinterpret_cast<Double_t const*>(valuePtr));
      }
      actualClass->Destructor(obj);
      return value;
    }
  }; // MetadataReader

  /// Stores `value` as a `TParameter<T>` with the specified `name` in `destDir`.
  /// `art::TFileDirectory::makeAndRegister()` can't be used since `TParameter`
  /// provides neither `SetName()` nor `SetTitle()`; instead, the parameter is
  /// registered on construction, in the directory `art::TFileDirectory::make()`
  /// has made current.
  template <typename T>
  void storeParameter(art::TFileDirectory& destDir, std::string const& name, T value)
  {
    struct Registrar {
      Registrar(std::string const& name, T value)
      {
        gDirectory->Append(new TParameter<T>(name.c_str(), value));
      }
    };
    delete destDir.make<Registrar>(name, value);
  }

} // local namespace

namespace phot {

  std::string const PhotonLibrary::OpChannelBranchName = "OpChannel";

  //------------------------------------------------------------
  PhotonLibrary::PhotonLibrary(art::TFileDirectory* pDir /* = nullptr */) : fDir(pDir) {}

  //------------------------------------------------------------

  void PhotonLibrary::StoreLibraryToFile(std::string,
                                         bool storeReflected,
                                         bool storeReflT0,
                                         size_t storeTiming) const
  {
    mf::LogInfo("PhotonLibrary") << "Writing photon library to file";

    if (!fDir) {
      throw cet::exception("PhotonLibrary")
        << "StoreLibraryToFile(): no ROOT file provided, can't store anything.\n";
    }

    TTree* tt = fDir->make<TTree>("PhotonLibraryData", "PhotonLibraryData");

    Int_t Voxel = 0;
    Int_t OpChannel = 0;
    Float_t Visibility = 0;
    Float_t ReflVisibility = 0;
    Float_t ReflTfirst = 0;
    Float_t* timing_par = nullptr;

    tt->Branch("Voxel", &Voxel, "Voxel/I");
    tt->Branch(OpChannelBranchName.c_str(), &OpChannel, (OpChannelBranchName + "/I").c_str());
    tt->Branch("Visibility", &Visibility, "Visibility/F");

    if (storeTiming != 0) {
      if (!hasTiming()) {
        // if this happens, you need to call CreateEmptyLibrary() with storeReflected set true
        throw cet::exception("PhotonLibrary")
          << "StoreLibraryToFile() requested to store the time propagation distribution "
             "parameters, which was not simulated.";
      }
      tt->Branch("timing_par", timing_par, Form("timing_par[%i]/F", size_t2int(storeTiming)));
    }

    if (storeReflected) {
      if (!hasReflected()) {
        // if this happens, you need to call CreateEmptyLibrary() with storeReflected set true
        throw cet::exception("PhotonLibrary")
          << "StoreLibraryToFile() requested to store reflected light, which was not simulated.";
      }
      tt->Branch("ReflVisibility", &ReflVisibility, "ReflVisibility/F");
    }
    if (storeReflT0) {
      if (!hasReflectedT0()) {
        // if this happens, you need to call CreateEmptyLibrary() with storeReflectedT0 set true
        throw cet::exception("PhotonLibrary") << "StoreLibraryToFile() requested to store "
                                                 "reflected light timing, which was not simulated.";
      }
      tt->Branch("ReflTfirst", &ReflTfirst, "ReflTfirst/F");
    }
    for (size_t ivox = 0; ivox != fNVoxels; ++ivox) {
      for (size_t ichan = 0; ichan != fNOpChannels; ++ichan) {
        Visibility = uncheckedAccess(ivox, ichan);
        if (storeReflected) ReflVisibility = uncheckedAccessRefl(ivox, ichan);
        if (storeReflT0) ReflTfirst = uncheckedAccessReflT(ivox, ichan);
        if (storeTiming != 0) {
          for (size_t i = 0; i < storeTiming; i++)
            timing_par[i] = uncheckedAccessTimingPar(ivox, ichan, i);
        }
        if (Visibility > 0 || ReflVisibility > 0) {
          Voxel = ivox;
          OpChannel = ichan;
          // visibility(ies) is(are) already set
          tt->Fill();
        }
      }
    }

    StoreMetadata();
  }

  //------------------------------------------------------------

  void PhotonLibrary::CreateEmptyLibrary(size_t NVoxels,
                                         size_t NOpChannels,
                                         bool storeReflected /* = false */,
                                         bool storeReflT0 /* = false */,
                                         size_t storeTiming /* = false */
  )
  {
    fNVoxels = NVoxels;
    fNOpChannels = NOpChannels;

    fHasReflected = storeReflected;
    fHasReflectedT0 = storeReflT0;
    fHasTiming = storeTiming;
    fTimingParNParameters = storeTiming;
    fTimingParFormula.clear(); // no formula: GetTimingTF1s() is not available

    resetStorage(storeReflected, storeReflT0, storeTiming);
  }

  //------------------------------------------------------------

  void PhotonLibrary::resetStorage(bool storeReflected, bool storeReflT0, size_t storeTiming)
  {
    if (fNVoxels > size_t(std::numeric_limits<std::int32_t>::max())) {
      throw cet::exception("PhotonLibrary")
        << "Photon library with " << fNVoxels << " voxels: at most "
        << std::numeric_limits<std::int32_t>::max() << " are supported.\n";
    }

    fVoxelTableIdx.assign(fNVoxels, NotFoundIdx);
    fVoxelTableIdx.shrink_to_fit();
    fNFilledVoxels = 0;

    fCounts.reset(fNOpChannels);
    fReflCounts.reset(storeReflected ? fNOpChannels : 0);
    fReflT0s.reset(storeReflT0 ? fNOpChannels : 0);
    fTimingPars.reset(fNOpChannels * storeTiming);

    fZeroRow.assign(
      std::max({fCounts.width(), fReflCounts.width(), fReflT0s.width(), fTimingPars.width()}),
      0.0f);

    fTimingTF1Row.clear();
    fTimingParRow.clear();
    invalidateTimingRows();
  }

  //------------------------------------------------------------

  size_t PhotonLibrary::tableIdxForWriting(size_t Voxel)
  {
    std::int32_t& tableIdx = fVoxelTableIdx[Voxel];
    if (tableIdx == NotFoundIdx) {
      if (fNFilledVoxels >= size_t(std::numeric_limits<std::int32_t>::max())) {
        throw cet::exception("PhotonLibrary") << "Too many filled voxels in the library.\n";
      }
      for (PropertyTable* table : {&fCounts, &fReflCounts, &fReflT0s, &fTimingPars})
        table->addRow(fNFilledVoxels);
      tableIdx = static_cast<std::int32_t>(fNFilledVoxels++);
    }
    return static_cast<size_t>(tableIdx);
  }

  //------------------------------------------------------------

  void PhotonLibrary::LoadLibraryFromFile(std::string LibraryFile,
                                          size_t NVoxels,
                                          bool getReflected,
                                          bool getReflT0,
                                          size_t getTiming,
                                          int maxRange)
  {
    mf::LogInfo("PhotonLibrary") << "Reading photon library from input file: "
                                 << LibraryFile.c_str() << std::endl;

    TFile* f = nullptr;
    TTree* tt = nullptr;
    TDirectory* pSrcDir = nullptr;

    try {
      f = TFile::Open(LibraryFile.c_str());
      tt = (TTree*)f->Get("PhotonLibraryData");
      if (tt) { pSrcDir = f; }
      else { // Library not in the top directory
        TKey* key = f->FindKeyAny("PhotonLibraryData");
        if (key) {
          tt = (TTree*)key->ReadObj();
          pSrcDir = key->GetMotherDir();
        }
        else {
          mf::LogError("PhotonLibrary") << "PhotonLibraryData not found in file" << LibraryFile;
        }
      }
    }
    catch (...) {
      throw cet::exception("PhotonLibrary")
        << "Error in ttree load, reading photon library: " << LibraryFile << "\n";
    }

    Int_t Voxel;
    Int_t OpChannel;
    Float_t Visibility;
    Float_t ReflVisibility;
    Float_t ReflTfirst;
    std::vector<Float_t> timing_par;

    tt->SetBranchAddress("Voxel", &Voxel);
    tt->SetBranchAddress("OpChannel", &OpChannel);
    tt->SetBranchAddress("Visibility", &Visibility);

    fHasTiming = getTiming;

    fHasReflected = getReflected;
    if (getReflected) tt->SetBranchAddress("ReflVisibility", &ReflVisibility);
    fHasReflectedT0 = getReflT0;
    if (getReflT0) tt->SetBranchAddress("ReflTfirst", &ReflTfirst);

    fNVoxels = NVoxels;
    fNOpChannels = PhotonLibrary::ExtractNOpChannels(tt); // EXPENSIVE!!!

    fTimingParNParameters = getTiming;
    fTimingMaxRange = maxRange;
    fTimingParFormula.clear();
    if (fHasTiming != 0) {
      timing_par.resize(getTiming);
      tt->SetBranchAddress("timing_par", timing_par.data());
      // should be pSrcDir->Get()? kept as is for backward compatibility
      TNamed* n = (TNamed*)f->Get("fTimingParFormula");
      if (!n) {
        throw cet::exception("PhotonLibrary")
          << "Error reading the photon propagation formula 'fTimingParFormula' from " << LibraryFile
          << ". Please check the photon library.\n";
      }
      fTimingParFormula = n->GetTitle();
      mf::LogInfo("PhotonLibrary")
        << "Time parametrization is activated. Using the formula: " << fTimingParFormula << " with "
        << fTimingParNParameters << " parameters." << std::endl;
    }

    // only the voxels with entries in the file are filled (see `fVoxelTableIdx`)
    resetStorage(getReflected, getReflT0, getTiming);

    size_t NEntries = tt->GetEntries();

    for (size_t i = 0; i != NEntries; ++i) {

      tt->GetEntry(i);

      if ((Voxel < 0) || (size_t(Voxel) >= fNVoxels) || (OpChannel < 0) ||
          (size_t(OpChannel) >= fNOpChannels)) {
        throw cet::exception("PhotonLibrary")
          << "Entry " << i << " of the photon library " << LibraryFile << " has voxel " << Voxel
          << " and channel " << OpChannel << ", outside the " << fNVoxels << " voxels x "
          << fNOpChannels << " channels expected.\n";
      }

      // Set the visibility at this optical channel
      uncheckedAccess(Voxel, OpChannel) = Visibility;

      if (fHasReflected) uncheckedAccessRefl(Voxel, OpChannel) = ReflVisibility;
      if (fHasReflectedT0) uncheckedAccessReflT(Voxel, OpChannel) = ReflTfirst;
      // the parameters are stored as they are; see `fillTimingTF1Row()` for their use
      for (size_t k = 0; k < fTimingParNParameters; ++k)
        uncheckedAccessTimingPar(Voxel, OpChannel, k) = timing_par[k];
    } // for entries

    if (fHasTiming != 0) makeTimingTF1Row();

    LoadMetadata(*pSrcDir);
    {
      size_t const allocated = fCounts.allocatedBytes() + fReflCounts.allocatedBytes() +
                               fReflT0s.allocatedBytes() + fTimingPars.allocatedBytes() +
                               fVoxelTableIdx.size() * sizeof(fVoxelTableIdx.front());
      mf::LogInfo log("PhotonLibrary");
      log << "Photon lookup table size : " << NVoxels << " voxels,  " << fNOpChannels
          << " channels";
      if (hasVoxelDef())
        log << "; " << GetVoxelDef();
      else
        log << " (no voxel geometry included)";
      log << "; " << fNFilledVoxels << " voxels filled, " << (allocated / 1048576)
          << " MiB allocated";
    }

    try {
      f->Close();
    }
    catch (...) {
      mf::LogError("PhotonLibrary") << "Error in closing file : " << LibraryFile;
    }
  }

  //----------------------------------------------------

  float PhotonLibrary::GetCount(size_t Voxel, size_t OpChannel) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccess(Voxel, OpChannel);
  }
  //----------------------------------------------------

  float PhotonLibrary::GetTimingPar(size_t Voxel, size_t OpChannel, size_t parnum) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccessTimingPar(Voxel, OpChannel, parnum);
  } //----------------------------------------------------

  float PhotonLibrary::GetReflCount(size_t Voxel, size_t OpChannel) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccessRefl(Voxel, OpChannel);
  }
  //----------------------------------------------------

  float PhotonLibrary::GetReflT0(size_t Voxel, size_t OpChannel) const
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      return 0;
    else
      return uncheckedAccessReflT(Voxel, OpChannel);
  }

  //----------------------------------------------------

  void PhotonLibrary::SetCount(size_t Voxel, size_t OpChannel, float Count)
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")
        << "Error - attempting to set count in voxel " << Voxel << " which is out of range";
    else
      uncheckedAccess(Voxel, OpChannel) = Count;
  }
  //----------------------------------------------------

  void PhotonLibrary::SetTimingPar(size_t Voxel, size_t OpChannel, float Count, size_t parnum)
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary") << "Error - attempting to set timing t0 count in voxel "
                                    << Voxel << " which is out of range";
    else if (parnum >= fTimingParNParameters)
      mf::LogError("PhotonLibrary") << "Error - attempting to set timing parameter " << parnum
                                    << " but the library has " << fTimingParNParameters;
    else
      uncheckedAccessTimingPar(Voxel, OpChannel, parnum) = Count;
  }
  //----------------------------------------------------

  void PhotonLibrary::SetTimingTF1(size_t Voxel, size_t OpChannel, TF1 func)
  {
    // Only the parameters are stored, in the layout of the library files:
    // the lower end of the range takes the place of parameter 0 (see `fillTimingTF1Row()`).
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary") << "Error - attempting to set a propagation function in voxel "
                                    << Voxel << " which is out of range";
    else if (fTimingParNParameters == 0)
      mf::LogError("PhotonLibrary")
        << "Error - attempting to set a propagation function in a library without timing";
    else {
      uncheckedAccessTimingPar(Voxel, OpChannel, 0) = func.GetXmin();
      for (size_t k = 1; k < fTimingParNParameters; ++k)
        uncheckedAccessTimingPar(Voxel, OpChannel, k) = func.GetParameter(k);
    }
  }
  //----------------------------------------------------

  void PhotonLibrary::SetReflCount(size_t Voxel, size_t OpChannel, float Count)
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")
        << "Error - attempting to set count in voxel " << Voxel << " which is out of range";
    else
      uncheckedAccessRefl(Voxel, OpChannel) = Count;
  }
  //----------------------------------------------------

  void PhotonLibrary::SetReflT0(size_t Voxel, size_t OpChannel, float Count)
  {
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")
        << "Error - attempting to set count in voxel " << Voxel << " which is out of range";
    else
      uncheckedAccessReflT(Voxel, OpChannel) = Count;
  }

  //----------------------------------------------------

  float const* PhotonLibrary::GetCounts(size_t Voxel) const
  {
    if (Voxel >= fNVoxels)
      return nullptr;
    else
      return valuesOf(fCounts, Voxel);
  }

  //----------------------------------------------------

  const std::vector<float>* PhotonLibrary::GetTimingPars(size_t Voxel) const
  {
    if ((Voxel >= fNVoxels) || (fHasTiming == 0)) return nullptr;
    fillTimingParRow(Voxel);
    return fTimingParRow.data();
  }

  //----------------------------------------------------

  TF1* PhotonLibrary::GetTimingTF1s(size_t Voxel) const
  {
    if ((Voxel >= fNVoxels) || (fHasTiming == 0)) return nullptr;
    if (fTimingTF1Row.empty()) {
      throw cet::exception("PhotonLibrary")
        << "GetTimingTF1s(): no time parametrization formula is available; the timing functions "
           "can be used only with a library loaded from a file with time parametrization.\n";
    }
    // the functions are not constant because `TF1::GetRandom()` is not (ROOT-9549)
    fillTimingTF1Row(Voxel);
    return fTimingTF1Row.data();
  }

  //----------------------------------------------------

  void PhotonLibrary::makeTimingTF1Row()
  {
    TF1 const prototype("PhotonLibraryTiming", fTimingParFormula.c_str(), 0.0, fTimingMaxRange);
    fTimingTF1Row.assign(fNOpChannels, prototype);
    fTimingTF1RowVoxel.reset();
  }

  //----------------------------------------------------

  void PhotonLibrary::fillTimingTF1Row(size_t Voxel) const
  {
    if (fTimingTF1RowVoxel == Voxel) return; // keeps the integrals cached by `GetRandom()`
    float const* pars = valuesOf(fTimingPars, Voxel);
    for (size_t channel = 0; channel < fNOpChannels; ++channel) {
      float const* channelPars = pars + channel * fTimingParNParameters;
      TF1& function = fTimingTF1Row[channel];
      // the first parameter is the lower end of the range, and the parameter itself is 0,
      // for compatibility with old libraries
      function.SetRange(channelPars[0], fTimingMaxRange);
      function.SetParameter(0, 0.0);
      for (size_t k = 1; k < fTimingParNParameters; ++k)
        function.SetParameter(k, channelPars[k]);
    }
    fTimingTF1RowVoxel = Voxel;
  }

  //----------------------------------------------------

  void PhotonLibrary::fillTimingParRow(size_t Voxel) const
  {
    if (fTimingParRowVoxel == Voxel) return;
    float const* pars = valuesOf(fTimingPars, Voxel);
    fTimingParRow.resize(fNOpChannels);
    for (size_t channel = 0; channel < fNOpChannels; ++channel) {
      float const* channelPars = pars + channel * fTimingParNParameters;
      fTimingParRow[channel].assign(channelPars, channelPars + fTimingParNParameters);
    }
    fTimingParRowVoxel = Voxel;
  }

  //----------------------------------------------------

  float const* PhotonLibrary::GetReflCounts(size_t Voxel) const
  {
    if ((Voxel >= fNVoxels) || !fHasReflected)
      return nullptr;
    else
      return valuesOf(fReflCounts, Voxel);
  }

  //----------------------------------------------------

  float const* PhotonLibrary::GetReflT0s(size_t Voxel) const
  {
    if ((Voxel >= fNVoxels) || !fHasReflectedT0)
      return nullptr;
    else
      return valuesOf(fReflT0s, Voxel);
  }

  //------------------------------------------------------------
  void PhotonLibrary::LoadMetadata(TDirectory& srcDir)
  {

    constexpr std::size_t NExpectedKeys = 9U;

    std::vector<std::string> missingKeys;

    MetadataReader<Int_t> readInt{srcDir, missingKeys};
    MetadataReader<Double_t> readDouble{srcDir, missingKeys};

    // FIXME: Not sure initializing to 0 is appropriate.
    double xMin{};
    double xMax{};
    int xN{};
    double yMin{};
    double yMax{};
    int yN{};
    double zMin{};
    double zMax{};
    int zN{};
    if (auto metaValue = readDouble("MinX")) xMin = *metaValue;
    if (auto metaValue = readDouble("MaxX")) xMax = *metaValue;
    if (auto metaValue = readInt("NDivX")) xN = *metaValue;
    if (auto metaValue = readDouble("MinY")) yMin = *metaValue;
    if (auto metaValue = readDouble("MaxY")) yMax = *metaValue;
    if (auto metaValue = readInt("NDivY")) yN = *metaValue;
    if (auto metaValue = readDouble("MinZ")) zMin = *metaValue;
    if (auto metaValue = readDouble("MaxZ")) zMax = *metaValue;
    if (auto metaValue = readInt("NDivZ")) zN = *metaValue;

    if (!missingKeys.empty()) {
      if (missingKeys.size() != NExpectedKeys) {
        mf::LogError log("PhotonLibrary");
        log << "Photon library at '" << srcDir.GetPath() << "' is missing " << missingKeys.size()
            << " metadata elements:";
        for (auto const& key : missingKeys)
          log << " '" << key << "'";
      }
      else {
        mf::LogTrace("PhotonLibrary") << "No voxel metadata found in '" << srcDir.GetPath() << "'";
      }
      return;
    } // if missing keys

    fVoxelDef.emplace(xMin, xMax, xN, yMin, yMax, yN, zMin, zMax, zN);

  } // PhotonLibrary::LoadMetadata()

  //------------------------------------------------------------
  void PhotonLibrary::StoreMetadata() const
  {

    assert(fDir);

    // total number of voxels and of optical detector channels in the library
    storeParameter<Int_t>(*fDir, "NVoxels", size_t2int(fNVoxels));
    storeParameter<Int_t>(*fDir, "NChannels", size_t2int(fNOpChannels));

    if (!hasVoxelDef()) return;
    sim::PhotonVoxelDef const& voxelDef = GetVoxelDef();

    // lower point covered by the library (world coordinates, cm)
    geo::Point_t const& lower = voxelDef.GetRegionLowerCorner();
    storeParameter<Double_t>(*fDir, "MinX", lower.X());
    storeParameter<Double_t>(*fDir, "MinY", lower.Y());
    storeParameter<Double_t>(*fDir, "MinZ", lower.Z());

    // upper point covered by the library (world coordinates, cm)
    geo::Point_t const& upper = voxelDef.GetRegionUpperCorner();
    storeParameter<Double_t>(*fDir, "MaxX", upper.X());
    storeParameter<Double_t>(*fDir, "MaxY", upper.Y());
    storeParameter<Double_t>(*fDir, "MaxZ", upper.Z());

    // size of a voxel on each direction (cm)
    geo::Vector_t const& stepSizes = voxelDef.GetVoxelSize();
    storeParameter<Double_t>(*fDir, "StepX", stepSizes.X());
    storeParameter<Double_t>(*fDir, "StepY", stepSizes.Y());
    storeParameter<Double_t>(*fDir, "StepZ", stepSizes.Z());

    // steps on each direction
    auto const& steps = voxelDef.GetSteps();
    storeParameter<Int_t>(*fDir, "NDivX", steps[0]);
    storeParameter<Int_t>(*fDir, "NDivY", steps[1]);
    storeParameter<Int_t>(*fDir, "NDivZ", steps[2]);

  } // PhotonLibrary::StoreMetadata()

  //----------------------------------------------------

  size_t PhotonLibrary::ExtractNOpChannels(TTree* tree)
  {
    TBranch* channelBranch = tree->GetBranch(OpChannelBranchName.c_str());
    if (!channelBranch) {
      throw art::Exception(art::errors::NotFound)
        << "Tree '" << tree->GetName() << "' has no branch 'OpChannel'";
    }

    // fix a new local address for the branch
    char* oldAddress = channelBranch->GetAddress();
    Int_t channel;
    channelBranch->SetAddress(&channel);
    Int_t maxChannel = -1;

    // read all the channel values and pick the largest one
    Long64_t iEntry = 0;
    while (channelBranch->GetEntry(iEntry++)) {
      if (channel > maxChannel) maxChannel = channel;
    } // while

    MF_LOG_DEBUG("PhotonLibrary") << "Detected highest channel to be " << maxChannel << " from "
                                  << iEntry << " tree entries";

    // restore the old branch address
    channelBranch->SetAddress(oldAddress);

    return size_t(maxChannel + 1);

  } // PhotonLibrary::ExtractNOpChannels()

  //----------------------------------------------------

}
