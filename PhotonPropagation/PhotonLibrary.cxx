
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
  
  std::vector<float> PhotonLibrary::EmptyChannelsList; // used for invalid return value
  
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
    tt->Branch("OpChannel",  &OpChannel,  "OpChannel/I");
    tt->Branch("Visibility", &Visibility, "Visibility/F");
 

    for(size_t ivox=0; ivox!=fLookupTable.size(); ++ivox)
      {
	for(size_t ichan=0; ichan!=fLookupTable.at(ivox).size(); ++ichan)
	  {
	    if(fLookupTable[ivox].at(ichan) > 0)
	      {
		Voxel      = ivox;
		OpChannel  = ichan;
		Visibility = fLookupTable[ivox][ichan];
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

    fLookupTable.resize(NVoxels);    

    for(size_t ivox=0; ivox!=NVoxels; ivox++)
      {
        fLookupTable[ivox].resize(NOpChannels,0);
      }
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

    
    fLookupTable.resize(NVoxels);    


    size_t NEntries = tt->GetEntries();

    for(size_t i=0; i!=NEntries; ++i) {
      tt->GetEntry(i);

      // Set # of optical channels to 1 more than largest one seen
      if (OpChannel >= fNOpChannels)
        fNOpChannels = OpChannel+1;

      // Expand this voxel's vector if needed
      if (fLookupTable[Voxel].size() < fNOpChannels)
        fLookupTable[Voxel].resize(fNOpChannels, 0);

      // Set the visibility at this optical channel
      fLookupTable[Voxel].at(OpChannel) = Visibility;
    }

    // Go through the table and fill in any missing 0's
    for(size_t ivox=0; ivox!=NVoxels; ivox++)
    {
      if (fLookupTable[ivox].size() < fNOpChannels)
	fLookupTable[ivox].resize(fNOpChannels,0);
    }
    
    
    mf::LogInfo("PhotonLibrary") <<"Photon lookup table size : "<<  NVoxels << " voxels,  " << NOpChannels<<" channels";


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

  float PhotonLibrary::GetCount(size_t Voxel, size_t OpChannel) 
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels)||/*(OpChannel<0)||*/(OpChannel>=fNOpChannels))
      return 0;   
    else
      return fLookupTable[Voxel].at(OpChannel); 
  }

  //----------------------------------------------------

  void PhotonLibrary::SetCount(size_t Voxel, size_t OpChannel, float Count) 
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels))
      mf::LogError("PhotonLibrary")<<"Error - attempting to set count in voxel " << Voxel<<" which is out of range"; 
    else
      fLookupTable[Voxel].at(OpChannel) = Count; 
  }

  //----------------------------------------------------

  const std::vector<float>* PhotonLibrary::GetCounts(size_t Voxel) const
  { 
    if(/*(Voxel<0)||*/(Voxel>=fNVoxels))
      return EmptyList(); // FIXME!!! better to throw an exception!
    else 
      return &(fLookupTable[Voxel]);
  }

  
  
}
