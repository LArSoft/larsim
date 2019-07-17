////////////////////////////////////////////////////////////////////////
//
// WienerFilterAna class designed to calculate the optimum filter for an event
//  (based strongly on CalWireAna)
//  andrzej.szelc@yale.edu
//
//
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <string>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include <TH1.h>

namespace detsim {

  /// Base class for creation of raw signals on wires.
  class WienerFilterAna : public art::EDAnalyzer {

  public:

    explicit WienerFilterAna(fhicl::ParameterSet const& pset);

    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void endJob();

  private:

    std::string            fDetSimModuleLabel; //< name of module that produced the digits

    TH1F*		   fCnoise[10][10][5];
    TH1F*		   fCsignal[10][10][5];

    TH1F*		   fCnoise_av[10][10][5];
    TH1F*		   fCsignal_av[10][10][5];

    TH1F*		   fFilter_av[10][10][5];

    TH1*		   ff;
    TH1F*		   hh;
    int 		  fNBins;
  }; // class WienerFilterAna

} // End caldata namespace.

namespace detsim{

  //-------------------------------------------------
  WienerFilterAna::WienerFilterAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fDetSimModuleLabel(pset.get< std::string >("DetSimModuleLabel"))
  {

  }

  //-------------------------------------------------
  void WienerFilterAna::beginJob()
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;

    art::ServiceHandle<util::LArFFT const> fFFT;
    int fNTicks = fFFT->FFTSize();
    fNBins=fNTicks/2+1;
    const detinfo::DetectorProperties* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double samprate=detp->SamplingRate();
    double sampfreq=1./samprate *1e6; // in kHz
    art::ServiceHandle<geo::Geometry const> geo;
    unsigned int fNPlanes = geo->Nplanes();
    unsigned int fNCryostats=geo->Ncryostats();
    unsigned int fNTPC=geo->NTPC();


   for(unsigned int icstat=0;icstat<fNCryostats;icstat++){
    for(unsigned int itpc=0;itpc<fNTPC;itpc++){
     for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
	{
	fCnoise[icstat][itpc][iplane]=tfs->make<TH1F>(Form("fft_noise_%d_%d_%d",icstat,itpc,iplane), 					Form("fft_of_unfiltered_noise_%d_%d_%d",icstat,itpc,iplane),fNBins,0,sampfreq/2);

	fCsignal[icstat][itpc][iplane]=tfs->make<TH1F>(Form("fft_signal_%d_%d_%d",icstat,itpc,iplane),Form("fft_of_unfiltered_noise_and_signal_%d_%d_%d",icstat,itpc,iplane),fNBins,0,sampfreq/2);

	fCnoise_av[icstat][itpc][iplane]=tfs->make<TH1F>(Form("fft_noise_%d_%d_%d_av",icstat,itpc,iplane),Form("fft_of_unfiltered_noise_%d_%d_%d_av",icstat,itpc,iplane),fNBins,0,sampfreq/2);
	fCsignal_av[icstat][itpc][iplane]=tfs->make<TH1F>(Form("fft_signal_%d_%d_%d_av",icstat,itpc,iplane),Form("fft_of_unfiltered_noise_and_signal_%d_%d_%d_av",icstat,itpc,iplane),fNBins,0,sampfreq/2);

	fFilter_av[icstat][itpc][iplane]=tfs->make<TH1F>(Form("fft_filter_%d_%d_%d_av",icstat,itpc,iplane),Form("fft_filter_%d_%d_%d_av",icstat,itpc,iplane),fNBins,0,sampfreq/2);
	}
      }
    }

    //ff=tfs->make<TH1>(Form("fftwaveform"),Form("fftwaveform"),4096,0,4096);
    hh=tfs->make<TH1F>(Form("waveform"),Form("waveform"),fNTicks,0,fNTicks);


    return;

  }

  //-------------------------------------------------
  void WienerFilterAna::endJob()
  {

     art::ServiceHandle<geo::Geometry const> geom;
    unsigned int nplanes=geom->Nplanes();
    unsigned int fNCryostats=geom->Ncryostats();
    unsigned int fNTPC=geom->NTPC();

     // calculate filters

   for(unsigned int icstat=0;icstat<fNCryostats;icstat++){
    for(unsigned int itpc=0;itpc<fNTPC;itpc++){
    for (unsigned int pp=0;pp<nplanes;pp++)
    {

      for(int ii=1;ii<fCsignal_av[icstat][itpc][pp]->GetNbinsX();ii++)
      {


      double diff=( (fCsignal_av[icstat][itpc][pp]->GetBinContent(ii)*fCsignal_av[icstat][itpc][pp]->GetBinContent(ii)-fCnoise_av[icstat][itpc][pp]->GetBinContent(ii)*fCnoise_av[icstat][itpc][pp]->GetBinContent(ii)) >=0 ) ? fCsignal_av[icstat][itpc][pp]->GetBinContent(ii)*fCsignal_av[icstat][itpc][pp]->GetBinContent(ii)-fCnoise_av[icstat][itpc][pp]->GetBinContent(ii)*fCnoise_av[icstat][itpc][pp]->GetBinContent(ii) : 0;




      if(fCsignal_av[icstat][itpc][pp]->GetBinContent(ii)>0)
	  fFilter_av[icstat][itpc][pp]->SetBinContent(ii,(double)((diff)/(fCsignal_av[icstat][itpc][pp]->GetBinContent(ii)*fCsignal_av[icstat][itpc][pp]->GetBinContent(ii))) );
      else
   	  fFilter_av[icstat][itpc][pp]->SetBinContent(ii,0);



	} // end loop on Csignal
      } // end loop on planes
    } // end loop on TPCs
   } // end loop on cryostats



  }

  //-------------------------------------------------
  void WienerFilterAna::analyze(const art::Event& evt)
  {

    // loop over the raw digits and get the adc vector for each, then compress it and uncompress it

    art::Handle< std::vector<raw::RawDigit> > rdHandle;
    evt.getByLabel(fDetSimModuleLabel,rdHandle);
    mf::LogInfo("WienerFilterMicroBooNE") << " readout  Wiener " << rdHandle->size() <<  std::endl;
   // return;
     if (!rdHandle->size())  return;
    mf::LogInfo("WienerFilterMicroBooNE") << "WienerFilterMicroBooNE:: rdHandle size is " << rdHandle->size();


    // Read in the digit List object(s).



    // Use the handle to get a particular (0th) element of collection.
    art::PtrVector<raw::RawDigit> rdvec;
    for(unsigned int i = 0; i < rdHandle->size(); ++i){
      art::Ptr<raw::RawDigit> r(rdHandle,i);
      rdvec.push_back(r);
     // std::cout << " i, rdvec: "<<i <<" " << r->ADC(0) << " "<< rdvec[i]->ADC(0)<< std::endl;
    }


    art::ServiceHandle<geo::Geometry const> geom;
    art::ServiceHandle<util::LArFFT const> fft;

    for(unsigned int rd = 0; rd < rdvec.size(); ++rd){

    std::vector<double> adc(fft->FFTSize());

      for(unsigned int t = 1; t < rdvec[rd]->Samples(); t++){
	adc[t-1]=rdvec[rd]->ADC(t-1);
	hh->SetBinContent(t,rdvec[rd]->ADC(t));
      }

      geo::WireID wireid= geom->ChannelToWire(rdvec[rd]->Channel())[0];


      // this is hardcoded for the time being. Should be automatized.
      unsigned int plane = wireid.Plane;   /// \todo  Need to change hardcoded values to an automatic determination of noise vs. signal
      unsigned int wire = wireid.Wire;
      unsigned int cstat = wireid.Cryostat;
      unsigned int tpc = wireid.TPC;
	  ff=hh->FFT(NULL,"MAG M");
	  if(wire>=50 && wire <250) {
	      for(int ii=0;ii<fNBins;ii++)  {
		fCnoise_av[cstat][tpc][plane]->AddBinContent(ii,ff->GetBinContent(ii));
		  if(wire==150)
		   fCnoise[cstat][tpc][plane]->SetBinContent(ii,ff->GetBinContent(ii));
	      }
	  }
	  else if(wire>=700 && wire <900) {
		for(int ii=0;ii<fNBins;ii++)  {
		  fCsignal_av[cstat][tpc][plane]->AddBinContent(ii,ff->GetBinContent(ii));
		  if(wire==800)
		    fCsignal[cstat][tpc][plane]->SetBinContent(ii,ff->GetBinContent(ii));
		}

	  }


//
     }//end loop over rawDigits

    return;
  }//end analyze method

}//end namespace

namespace detsim{

  DEFINE_ART_MODULE(WienerFilterAna)

}
