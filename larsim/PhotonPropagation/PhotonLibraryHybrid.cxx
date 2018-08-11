#include "larsim/PhotonPropagation/PhotonLibraryHybrid.h"

#include "larsim/Simulation/PhotonVoxels.h"

#include "larcore/Geometry/Geometry.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"

#include <iostream>

namespace
{
  bool fatal(const std::string& msg)
  {
    std::cerr << "FATAL: PhotonLibraryHybrid: " << msg << std::endl;
    abort();
  }
}

namespace phot
{
  //--------------------------------------------------------------------
  PhotonLibraryHybrid::PhotonLibraryHybrid(const std::string& fname,
                                           const sim::PhotonVoxelDef& voxdef)
    : fVoxDef(voxdef)
  {
    TFile f(fname.c_str());
    !f.IsZombie() || fatal("Could not open PhotonLibrary "+fname);

    for(int opdetIdx = 0; true; ++opdetIdx){
      const std::string dirname = TString::Format("opdet_%d", opdetIdx).Data();
      TDirectory* dir = (TDirectory*)f.Get(dirname.c_str());
      if(!dir) break; // Ran out of opdets

      OpDetRecord rec;

      TVectorD* fit = (TVectorD*)dir->Get("fit");
      fit || fatal("Didn't find "+dirname+"/fit in "+fname);
      rec.fit = FitFunc((*fit)[0], (*fit)[1]);

      TTree* tr = (TTree*)dir->Get("tr");
      tr || fatal("Didn't find "+dirname+"/tr in "+fname);
      int vox;
      float vis;
      tr->SetBranchAddress("vox", &vox);
      tr->SetBranchAddress("vis", &vis);

      rec.exceptions.reserve(tr->GetEntries());
      for(int i = 0; i < tr->GetEntries(); ++i){
        tr->GetEntry(i);
        vox < NVoxels() || fatal("Voxel out of range");
        rec.exceptions.emplace_back(vox, vis);
      }
      rec.exceptions.shrink_to_fit(); // In case we don't trust reserve()

      // In case they weren't sorted in the file
      std::sort(rec.exceptions.begin(), rec.exceptions.end());

      fRecords.push_back(rec);
    } // end for opdetIdx

    !fRecords.empty() || fatal("No opdet_*/ directories in "+fname);

    art::ServiceHandle<geo::Geometry> geom;
    geom->NOpDets() == fRecords.size() || fatal("Number of opdets mismatch");

    fRecords.shrink_to_fit(); // save memory
  }

  //--------------------------------------------------------------------
  PhotonLibraryHybrid::~PhotonLibraryHybrid()
  {
  }

  //--------------------------------------------------------------------
  int PhotonLibraryHybrid::NVoxels() const
  {
    return fVoxDef.GetNVoxels();
  }

  //--------------------------------------------------------------------
  const float* PhotonLibraryHybrid::GetCounts(size_t vox) const
  {
    // The concept of GetCounts() conflicts fairly badly with how this library
    // works. This is probably the best we can do.

    static float* counts = 0;
    if(!counts) counts = new float[NOpChannels()];
    for(int od = 0; od < NOpChannels(); ++od) counts[od] = GetCount(vox, od);
    return counts;
  }

  //--------------------------------------------------------------------
  float PhotonLibraryHybrid::GetCount(size_t vox, size_t opchan) const
  {
    int(vox) < NVoxels() || fatal("GetCount(): Voxel out of range");
    int(opchan) < NOpChannels() || fatal("GetCount(): OpChan out of range");

    const OpDetRecord& rec = fRecords[opchan];

    auto it2 = std::lower_bound(rec.exceptions.begin(),
                                rec.exceptions.end(),
                                vox);
    if(it2->vox == vox) return it2->vis;

    // Otherwise, we use the interpolation, which requires a distance

    static art::ServiceHandle<geo::Geometry> geom;
    const geo::OpDetGeo& opdet = geom->OpDetGeoFromOpDet(opchan);

    const TVector3 voxvec = fVoxDef.GetPhotonVoxel(vox).GetCenter();
    const double xyzvox[] = {voxvec.X(), voxvec.Y(), voxvec.Z()};

    const double dist = opdet.DistanceToPoint(xyzvox);

    return rec.fit.Eval(dist);
  }
}
