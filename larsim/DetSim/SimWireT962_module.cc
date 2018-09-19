////////////////////////////////////////////////////////////////////////
//
// SimWireT962 class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////
#ifndef SIMWIRET962_H
#define SIMWIRET962_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// ROOT includes
#include <TMath.h>
#include <TH1D.h>
#include <TFile.h>
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"

// C++ includes
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// ArgoNeuT includes
#include "T962Utilities/SignalShapingServiceT962.h"


namespace art {
  class Event;
  class ParameterSet;
}

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires. 
  class SimWireT962 : public art::EDProducer {
    
  public:
        
    explicit SimWireT962(fhicl::ParameterSet const& pset); 
    virtual ~SimWireT962();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:

    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    raw::Compress_t        fCompression;      ///< compression type to use
    double                 fNoiseFact;        ///< noise scale factor 
    double                 fNoiseWidth;       ///< exponential noise width (kHz) 
    double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
    int                    fNTicks;           ///< number of ticks of the clock

    double  fSampleRate;       ///< sampling rate in ns
    std::vector<std::vector<float>> fNoise;///< noise on each channel for each time
    void  GenNoise(std::vector<float>& array);
    TH1D* fNoiseDist;        ///< distribution of noise counts

  }; // class SimWireT962

}

namespace detsim{

  //-------------------------------------------------
  SimWireT962::SimWireT962(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();
 
    fCompression = raw::kNone;
    std::string compression(pset.get< std::string >("CompressionType"));
    if(compression.compare("Huffman") == 0) fCompression = raw::kHuffman;    

    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, pset, "Seed");
  }

  //-------------------------------------------------
  SimWireT962::~SimWireT962()
  {
  }

  //-------------------------------------------------
  void SimWireT962::reconfigure(fhicl::ParameterSet const& p) 
  {

    fDriftEModuleLabel= p.get<std::string>("DriftEModuleLabel");
    fNoiseFact        = p.get<double>("NoiseFact");
    fNoiseWidth       = p.get<double>("NoiseWidth");
    fLowCutoff        = p.get<double>("LowCutoff");
//    double ADCperFC   = p.get<double>("ADCperFC");
//    fADCperElectron = ADCperFC / 6250;
    return;
  }

  //-------------------------------------------------
  void SimWireT962::beginJob() 
  { 
    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    fNoiseDist = tfs->make<TH1D>("Noise", ";Noise (ADC);", 1000, -10., 10.);
    art::ServiceHandle<geo::Geometry> geo;
    fNoise.resize(geo->Nchannels());
    auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSampleRate = detprop->SamplingRate();
    // GenNoise() will further resize each channel's 
    // noise vector to fNTicks long.
    for(unsigned int p = 0; p < geo->Nchannels(); ++p){
      GenNoise(fNoise[p]);
      for(int i = 0; i < fNTicks; ++i) fNoiseDist->Fill(fNoise[p][i]);
    }// end loop over wires
  } // beginJob
  
  //-------------------------------------------------
  void SimWireT962::endJob() 
  {
  }

  //-------------------------------------------------
  void SimWireT962::produce(art::Event& evt)
  {

    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    const detinfo::DetectorClocks* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    unsigned int signalSize = fNTicks/2;
    // vectors for working
    std::vector<short> adcvec(signalSize, 0);	
    std::vector<double> chargeWork(fNTicks);

    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);
    
    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(geo->Nchannels());
    for(size_t c = 0; c < chanHandle.size(); ++c) channels[chanHandle[c]->Channel()] = chanHandle[c];

    // make an unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr<std::vector<raw::RawDigit>>  digcol(new std::vector<raw::RawDigit>);
	  
    art::ServiceHandle<util::SignalShapingServiceT962> sss;

    // Add all channels  
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    std::map<int,double>::iterator mapIter;      
    for(unsigned int chan = 0; chan < geo->Nchannels(); chan++) {
      
      std::fill(chargeWork.begin(), chargeWork.end(), 0.);

      if( channels[chan] ){
        // get the sim::SimChannel for this channel
        const sim::SimChannel* sc = channels[chan];

        // loop over the tdcs and grab the number of electrons for each
        for(size_t t = 0; t < chargeWork.size(); ++t){
          int tdc = ts->TPCTick2TDC(t);
          if (tdc < 0) continue;
          if (tdc >= int(chargeWork.size())) continue;
          chargeWork[t] = sc->Charge(tdc);
        } // t
        sss->Convolute(chan, chargeWork);
      } // channels[chan]
      
      if(chan == 99) {
        unsigned short lobin = 172;
        for(unsigned short bin = lobin; bin < lobin+100; ++bin) mf::LogVerbatim("SignalShapingServiceT962")<<chargeWork[bin];
      }

      // noise was already generated for each wire in the event
      // raw digit vec is already in channel order
      // pick a new "noise channel" for every channel  - this makes sure    
      // the noise has the right coherent characteristics to be on one channel
      int noisechan = TMath::Nint(flat.fire()*(1.*(fNoise.size()-1)+0.1));
      for(unsigned int i = 0; i < signalSize; ++i) {
        adcvec[i] = (short)TMath::Nint(fNoise[noisechan][i] + chargeWork[i]);
      } // i

      // compress the adc vector using the desired compression scheme,
      // if raw::kNone is selected nothing happens to adcvec
      // This shrinks adcvec, if fCompression is not kNone.
      raw::Compress(adcvec, fCompression);
      
      raw::RawDigit rd(chan, signalSize, adcvec, fCompression);
      // Then, resize adcvec back to full length!
      adcvec.clear();
      adcvec.resize(signalSize,0.0);
      // add this digit to the collection
      digcol->push_back(rd);
    }// end loop over channels      

    evt.put(std::move(digcol));
  } // produce
  
  //-------------------------------------------------
  void SimWireT962::GenNoise(std::vector<float>& noise)
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    noise.resize(fNTicks);
    std::fill(noise.begin(), noise.end(), 0.);
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);///<noise in frequency space

    double pval = 0.; 
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);
    for(int i=0; i< fNTicks/2+1; ++i){
      // exponential noise spectrum 
      pval = fNoiseFact*exp(-(double)i*binWidth/fNoiseWidth);
      // low frequency cutoff     
      lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
      // randomize 10%
      flat.fireArray(2,rnd,0,1);
      pval *= lofilter*(0.9+0.2*rnd[0]);
      // random pahse angle
      phase = rnd[1]*2.*TMath::Pi();

      TComplex tc(pval*cos(phase),pval*sin(phase));
      noiseFrequency[i] += tc;
    }
  
    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);

    noiseFrequency.clear();

    // multiply each noise value by fNTicks as the InvFFT 
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    for(unsigned int i = 0; i < noise.size(); ++i) noise[i] *= 1.*fNTicks;

  } // GenNoise

}

namespace detsim{

  DEFINE_ART_MODULE(SimWireT962)

}

#endif // SIMWIRET962_H
