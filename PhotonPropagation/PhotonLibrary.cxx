
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "Geometry/Geometry.h"

#include "PhotonPropagation/PhotonLibrary.h"
#include "Simulation/PhotonVoxels.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"


namespace phot{
  
  std::string const PhotonLibrary::OpChannelBranchName = "OpChannel";
  
  //------------------------------------------------------------

  PhotonLibrary::PhotonLibrary()
  {
    fLookupTable.clear();
  }

  
  //------------------------------------------------------------  

  PhotonLibrary::~PhotonLibrary()
  {
    fLookupTable.clear();
  }
  
  //------------------------------------------------------------
  
  void PhotonLibrary::StoreLibraryToFile(std::string LibraryFile)
  {
    mf::LogInfo("PhotonLibrary") << "Writing photon library to input file: " << LibraryFile.c_str()<<std::endl;

    art::ServiceHandle<art::TFileService> tfs;

    TTree *tt = tfs->make<TTree>("PhotonLibraryData","PhotonLibraryData");
 
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;


    tt->Branch("Voxel",      &Voxel,      "Voxel/I");
    tt->Branch(OpChannelBranchName.c_str(),  &OpChannel,  (OpChannelBranchName + "/I").c_str());
    tt->Branch("Visibility", &Visibility, "Visibility/F");
 

    for(size_t ivox=0; ivox!= fNVoxels; ++ivox)
      {
	for(size_t ichan=0; ichan!= fNOpChannels; ++ichan)
	  {
	    Visibility = uncheckedAccess(ivox, ichan);
	    if (Visibility > 0)
	      {
		Voxel      = ivox;
		OpChannel  = ichan;
		// visibility is already set
		tt->Fill();
	      }
	  }	
      }
  }


  //------------------------------------------------------------

  void PhotonLibrary::CreateEmptyLibrary( size_t NVoxels, size_t NOpChannels)
  {
    fLookupTable.clear();

    fNVoxels     = NVoxels;
    fNOpChannels = NOpChannels;

    fLookupTable.resize(LibrarySize(), 0.);

  }


  //------------------------------------------------------------

  void PhotonLibrary::LoadLibraryFromFile(std::string LibraryFile, size_t NVoxels)
  {
    fLookupTable.clear();
    
    mf::LogInfo("PhotonLibrary") << "Reading photon library from input file: " << LibraryFile.c_str()<<std::endl;

    TFile *f = nullptr;
    TTree *tt = nullptr;
      
    try
      {
	f  =  TFile::Open(LibraryFile.c_str());
	tt =  (TTree*)f->Get("PhotonLibraryData");
        if (!tt) { // Library not in the top directory
            TKey *key = f->FindKeyAny("PhotonLibraryData");
            if (key) 
                tt = (TTree*)key->ReadObj();
            else {
                mf::LogError("PhotonLibrary") << "PhotonLibraryData not found in file" <<LibraryFile;
            }
        }
      }
    catch(...)
      {
	mf::LogError("PhotonLibrary") << "Error in ttree load, reading photon library: " << LibraryFile.c_str()<<std::endl;
      }
    
    Int_t     Voxel;
    Int_t     OpChannel;
    Float_t   Visibility;

    tt->SetBranchAddress("Voxel",      &Voxel);
    tt->SetBranchAddress("OpChannel",  &OpChannel);
    tt->SetBranchAddress("Visibility", &Visibility);



    
    fNVoxels     = NVoxels;
    fNOpChannels = PhotonLibrary::ExtractNOpChannels(tt); // EXPENSIVE!!!

    
    fLookupTable.resize(LibrarySize(), 0.);


    size_t NEntries = tt->GetEntries();

    for(size_t i=0; i!=NEntries; ++i) {
      tt->GetEntry(i);

      // Set the visibility at this optical channel
      uncheckedAccess(Voxel, OpChannel) = Visibility;
      
    } // for entries

    mf::LogInfo("PhotonLibrary") <<"Photon lookup table size : "<<  NVoxels << " voxels,  " << fNOpChannels<<" channels";


    try
      {
	f->Close();
      }
    catch(...)
      {
	mf::LogError("PhotonLibrary") << "Error in closing file : " << LibraryFile.c_str()<<std::endl;
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

  void PhotonLibrary::SetCount(size_t Voxel, size_t OpChannel, float Count) 
  { 
    if ((Voxel >= fNVoxels) || (OpChannel >= fNOpChannels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range"; 
    else
      uncheckedAccess(Voxel, OpChannel) = Count; 
  }

  //----------------------------------------------------

  float const* PhotonLibrary::GetCounts(size_t Voxel) const
  { 
    if (Voxel >= fNVoxels) return nullptr;
    else return fLookupTable.data() + uncheckedIndex(Voxel, 0);
  }

  //----------------------------------------------------
  
  size_t PhotonLibrary::ExtractNOpChannels(TTree* tree) {
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
    
    LOG_DEBUG("PhotonLibrary")
      << "Detected highest channel to be " << maxChannel << " from " << iEntry
      << " tree entries";
    
    // restore the old branch address
    channelBranch->SetAddress(oldAddress);
    
    return size_t(maxChannel + 1);
    
  } // PhotonLibrary::ExtractNOpChannels()
  
}
