////////////////////////////////////////////////////////////////////////
//
// SimWire class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/Utilities/LArFFT.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH2.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


///Detector simulation of raw signals on wires
namespace detsim {

  /// Base class for creation of raw signals on wires. 
  class SimWireAna : public art::EDAnalyzer {
    
  public:
        
    explicit SimWireAna(fhicl::ParameterSet const& pset); 
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();

  private:

    std::string            fDetSimModuleLabel;///< name of module that produced the digits
    TH1F*                  fDiffs;            ///< histogram of Raw tdc to tdc differences

    TH1F*                  fCompressErr;      ///< histogram of difference between original 
                                              ///<tdc value and compressesed value
    TH1F*                  fCompressFactor;   ///< compression factor 

    TH2F*                  fRawVsCompress;    ///< histogram of original tdc value vs compressesed value
    TH2F*                  fCompressErr2D;    ///< histogram of original tdc value vs compressesed value
    

  }; // class SimWire

}

namespace detsim{

  //-------------------------------------------------
  SimWireAna::SimWireAna(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset)
    , fDetSimModuleLabel{pset.get< std::string >("DetSimModuleLabel")}
  {}

  //-------------------------------------------------
  void SimWireAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fDiffs          = tfs->make<TH1F>("One timestamp diffs", ";#Delta ADC;",        40,   -19.5,  20.5);
    fCompressErr    = tfs->make<TH1F>("compressErr",         ";Raw-Compressed;",    1000, -495.5, 500.5);
    fCompressFactor = tfs->make<TH1F>("compressFactor",      ";Compression;",       500,     0.,    1.);

    fCompressErr2D  = tfs->make<TH2F>("compressErr2D",       ";Raw;Raw-Compressed", 100, -50., 50., 1000,  -495.5, 500.5);
    fRawVsCompress  = tfs->make<TH2F>("rawVsCompress",       ";Raw;Compressed",     100, -50., 50.,  100,   -50.,  50.);
  
    return;

  }

  //-------------------------------------------------
  void SimWireAna::analyze(const art::Event& evt)
  {

    // loop over the raw digits and get the adc vector for each, then compress it and uncompress it

    art::Handle< std::vector<raw::RawDigit> > rdHandle;
    evt.getByLabel(fDetSimModuleLabel,rdHandle);

    art::PtrVector<raw::RawDigit> rdvec;
    for(unsigned int i = 0; i < rdHandle->size(); ++i){
      art::Ptr<raw::RawDigit> r(rdHandle,i);
      rdvec.push_back(r);
    }

    /// loop over all the raw digits
    for(unsigned int rd = 0; rd < rdvec.size(); ++rd){

      std::vector<short> adc;
      std::vector<short> uncompressed(rdvec[rd]->Samples());
      for(unsigned int t = 1; t < rdvec[rd]->Samples(); ++t){
	fDiffs->Fill(rdvec[rd]->ADC(t) - rdvec[rd]->ADC(t-1));
	adc.push_back(rdvec[rd]->ADC(t-1));
      }
    
      //get the last one for the adc vector
      adc.push_back(rdvec[rd]->ADC(rdvec[rd]->Samples()-1));
    
      raw::Compress(adc, raw::kHuffman);
    
      fCompressFactor->Fill((1.*adc.size())/(1.*rdvec[rd]->Samples()));
    
      raw::Uncompress(adc, uncompressed, raw::kHuffman);
    
      if(uncompressed.size() != rdvec[rd]->Samples()){ 
	cet::exception("WrongSizeUncompress") 
	  << "uncompression does not produce same size vector as original: " 
	  << "original = " << rdvec[rd]->Samples() << " uncompress = " 
	  << uncompressed.size() << "\n";
      }
    
      for(unsigned int t = 0; t <  uncompressed.size(); ++t){
	//std::cout << t << " " << rdFE->ADC(t) << " " << uncompressed[t] << std::endl;
	if(uncompressed[t]-rdvec[rd]->ADC(t) > 1) 
	  mf::LogWarning("SimWireAna") << "problem with event "
				       << " time " << t << " ADC " << rdvec[rd]->ADC(t) 
				       << " uncompress " << uncompressed[t] 
				       << " channel " << rdvec[rd]->Channel();
	
	fCompressErr->Fill(uncompressed[t]-rdvec[rd]->ADC(t));
	fCompressErr2D->Fill(rdvec[rd]->ADC(t), uncompressed[t]-rdvec[rd]->ADC(t));
	fRawVsCompress->Fill(rdvec[rd]->ADC(t), uncompressed[t]);
      }
    }//end loop over digits

    return;
  }//end analyze method


}//end namespace

namespace detsim{

  DEFINE_ART_MODULE(SimWireAna)

}
