TChain * CreateChainFromList_opt(std::string ListFileName, std::string ChainName, std::string BaseDirectory, bool DoCheck=false);
void MakeCombinedFile(std::string FileList, std::string BaseDirectory, std::string OutputName);


void AssembleSingleFile(std::string FileList, std::string BaseDirectory, std::string OutputName)
{
  
  TFile *f = TFile::Open(OutputName.c_str(),"RECREATE");
  TChain * ch = CreateChainFromList_opt(FileList.c_str(),BaseDirectory.c_str(),"pmtresponse/PhotonLibraryData",false);
  
  Int_t Voxel, OpChannel;
  Float_t Visibility;
  ch->SetBranchAddress("Voxel",      &Voxel);
  ch->SetBranchAddress("OpChannel",  &OpChannel);
  ch->SetBranchAddress("Visibility", &Visibility);

  TTree * tt = new TTree("PhotonLibraryData","PhotonLibraryData");
  tt->Branch("Voxel",      &Voxel,      "Voxel/I");
  tt->Branch("OpChannel",  &OpChannel,  "OpChannel/I");
  tt->Branch("Visibility", &Visibility, "Visibility/F");
  
  for(int i=0; i!=ch->GetEntries(); ++i)
    {
      ch->GetEntry(i);
      tt->Fill();
    }
  
  f->Write();
  f->Close();
}




TChain * CreateChainFromList_opt(std::string ListFileName, std::string BaseDirectory, std::string ChainName, bool DoCheck)
{
  ifstream InputFile(ListFileName.c_str());
  std::string FileName;

  TChain * TheChain = new TChain(ChainName.c_str());
  if(!DoCheck)
    {
      while(getline(InputFile, FileName))
        {
	  FileName=BaseDirectory+FileName;
	  std::cout<<FileName.c_str()<<std::endl;
          TheChain->Add(FileName.c_str());
        }
    }
  else
    {
      while(getline(InputFile, FileName))
        {
	  FileName=BaseDirectory+FileName;
          TFile*f=TFile::Open(FileName.c_str());
          if(f->Get(ChainName.c_str()))
            {
	      std::cout<<FileName.c_str()<<std::endl;
              TheChain->Add(FileName.c_str());
            }
          else std::cout<<"Chain " <<ChainName << " not found in file " << FileName<<std::endl;
        }


    }
  return TheChain;
}



