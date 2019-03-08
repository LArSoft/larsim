////////////////////////////////////////////////////////////////////////
//
// SimWire class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////

// ROOT includes
#include <TMath.h>
#include <TH1D.h>
#include <TFile.h>
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>
#include <vector>
#include <string>


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"

namespace art {
  class Event;
  class ParameterSet;
}

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires. 
  class SimWire : public art::EDProducer {
  public:
        
    explicit SimWire(fhicl::ParameterSet const& pset); 
    virtual ~SimWire();
    
  private:

    // read/write access to event
    void produce(art::Event& evt) override;
    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const& p);

    void         ConvoluteResponseFunctions(); ///< convolute electronics and field response
    
    void         SetFieldResponse();           ///< response of wires to field
    void         SetElectResponse();           ///< response of electronics
    
    void         GenNoise(std::vector<float>& array, CLHEP::HepRandomEngine& engine);

    bool                   fResponseSet;      ///< flag of whether to set the response functions or not
    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    std::string            fResponseFile;     ///< response file for induction planes
    raw::Compress_t        fCompression;      ///< compression type to use

    double                 fNoiseFact;        ///< noise scale factor 
    double                 fNoiseWidth;       ///< exponential noise width (kHz) 
    double                 fLowCutoff;        ///< low frequency filter cutoff (kHz)
    int                    fNTicks;           ///< number of ticks of the clock
    int                    fNFieldBins;       ///< number of bins for field response
    double                 fSampleRate;       ///< sampling rate in ns
    unsigned int           fNSamplesReadout;  ///< number of ADC readout samples in 1 readout frame
    double                 fCol3DCorrection;  ///< correction factor to account for 3D path of 
                                              ///< electrons thru wires
    double                 fInd3DCorrection;  ///< correction factor to account for 3D path of 
                                              ///< electrons thru wires
    double                 fColFieldRespAmp;  ///< amplitude of response to field 
    double                 fIndFieldRespAmp;  ///< amplitude of response to field 
    std::vector<double>    fShapeTimeConst;   ///< time constants for exponential shaping
    int                    fTriggerOffset;    ///< (units of ticks) time of expected neutrino event
    unsigned int           fNElectResp;       ///< number of entries from response to use
    
    std::vector<double>    fColFieldResponse; ///< response function for the field @ collection plane
    std::vector<double>    fIndFieldResponse; ///< response function for the field @ induction plane
    std::vector<TComplex>  fColShape;         ///< response function for the field @ collection plane
    std::vector<TComplex>  fIndShape;         ///< response function for the field @ induction plane
    std::vector<double>    fChargeWork;
    std::vector<double>    fElectResponse;    ///< response function for the electronics
    std::vector< std::vector<float> > fNoise;///< noise on each channel for each time
    
    TH1D*                fIndFieldResp;     ///< response function for the field @ induction plane
    TH1D*                fColFieldResp;     ///< response function for the field @ collection plane
    TH1D*                fElectResp;        ///< response function for the electronics
    TH1D*                fColTimeShape;     ///< convoluted shape for field x electronics @ col plane
    TH1D*                fIndTimeShape;     ///< convoluted shape for field x electronics @ ind plane
    TH1D*                fNoiseDist;        ///< distribution of noise counts

    CLHEP::HepRandomEngine& fEngine;        ///< Random-number engine owned by art
  }; // class SimWire

}

namespace detsim{

  //-------------------------------------------------
  SimWire::SimWire(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    , fEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, pset, "Seed"))
  {
    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();
 
    fCompression = raw::kNone;
    std::string compression(pset.get< std::string >("CompressionType"));
    if(compression.compare("Huffman") == 0) fCompression = raw::kHuffman;
  }

  //-------------------------------------------------
  SimWire::~SimWire()
  {
    fColFieldResponse.clear();
    fIndFieldResponse.clear();
    fColShape.clear();
    fIndShape.clear();
    fChargeWork.clear();
    fElectResponse.clear();
 
    for(unsigned int i = 0; i < fNoise.size(); ++i) fNoise[i].clear();
    fNoise.clear();

    MF_LOG_WARNING("SimWire") << "SimWire is an example module that works for the "
			   << "MicroBooNE detector.  Each experiment should implement "
			   << "its own version of this module to simulate electronics "
			   << "response.";

  }

  //-------------------------------------------------
  void SimWire::reconfigure(fhicl::ParameterSet const& p) 
  {
    std::string fResponseFile = p.get<std::string>("ResponseFile", "");
    
    fResponseSet      = !fResponseFile.empty();
    if (fResponseSet) {
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(p.get<std::string>("ResponseFile"), fResponseFile);
    }
    else fResponseFile.clear();
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");
    fNoiseFact        = p.get< double              >("NoiseFact");
    fNoiseWidth       = p.get< double              >("NoiseWidth");
    fLowCutoff        = p.get< double              >("LowCutoff");
    fNFieldBins       = p.get< int                 >("FieldBins");
    fCol3DCorrection  = p.get< double              >("Col3DCorrection");
    fInd3DCorrection  = p.get< double              >("Ind3DCorrection");
    fColFieldRespAmp  = p.get< double              >("ColFieldRespAmp");
    fIndFieldRespAmp  = p.get< double              >("IndFieldRespAmp");
    fShapeTimeConst   = p.get< std::vector<double> >("ShapeTimeConst");

    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fSampleRate       = detprop->SamplingRate();
    fTriggerOffset    = detprop->TriggerOffset();
    fNSamplesReadout  = detprop->NumberTimeSamples();
  }

  //-------------------------------------------------
  void SimWire::beginJob() 
  { 
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist      = tfs->make<TH1D>("Noise", ";Noise (ADC);", 1000, -10., 10.);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();
    fChargeWork.resize(fNTicks, 0.);

    // Note the magic 100 here. Argo and uBooNe use NChannels.
    fNoise.resize(100);
    // GenNoise() will further resize each channel's 
    // fNoise vector to fNTicks long.

    for(int p = 0; p < 100; ++p){
      GenNoise(fNoise[p], fEngine);
      for(int i = 0; i < fNTicks; ++i){
	fNoiseDist->Fill(fNoise[p][i]);
      }
    }//end loop over wires

    ///set field response and electronics response, then convolute them
    SetFieldResponse();
    SetElectResponse();
    ConvoluteResponseFunctions();
  }

  //-------------------------------------------------
  void SimWire::produce(art::Event& evt)
  {
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    unsigned int signalSize = fNTicks;
    // vectors for working
    std::vector<short>    adcvec(signalSize, 0);	

    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(geo->Nchannels());
    for(size_t c = 0; c < chanHandle.size(); ++c){
      channels[chanHandle[c]->Channel()] = chanHandle[c];
    }
    
    // make an unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcol(new std::vector<raw::RawDigit>);
	  
    unsigned int chan = 0; 
    fChargeWork.clear();
    fChargeWork.resize(fNTicks, 0.);
	  
    art::ServiceHandle<util::LArFFT> fFFT;

    // Add all channels  
    CLHEP::RandFlat flat(fEngine);

    std::map<int,double>::iterator mapIter;      
    for(chan = 0; chan < geo->Nchannels(); chan++) {    
      //       std::cout << "on channel " << chan << std::endl;
      
      fChargeWork.clear();    
      fChargeWork.resize(fNTicks, 0.);    

      if( channels[chan] ){      

	// get the sim::SimChannel for this channel
	const sim::SimChannel* sc = channels[chan];

	// loop over the tdcs and grab the number of electrons for each
	for(size_t t = 0; t < fChargeWork.size(); ++t) 
	  fChargeWork[t] = sc->Charge(t);      

        //Convolve charge with appropriate response function 
	if(geo->SignalType(chan) == geo::kInduction)        
	  fFFT->Convolute(fChargeWork,fIndShape);      
	else fFFT->Convolute(fChargeWork,fColShape);      
      }

      // noise was already generated for each wire in the event
      // raw digit vec is already in channel order
      // pick a new "noise channel" for every channel  - this makes sure    
      // the noise has the right coherent characteristics to be on one channel   
      int noisechan = TMath::Nint(flat.fire()*(1.*(fNoise.size()-1)+0.1));
      for(unsigned int i = 0; i < signalSize; ++i){
	adcvec[i] = (short)TMath::Nint(fNoise[noisechan][i] + fChargeWork[i]);
      }
      adcvec.resize(fNSamplesReadout);

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
    }//end loop over channels      

    evt.put(std::move(digcol));

    return;
  }

  //-------------------------------------------------
  void SimWire::ConvoluteResponseFunctions()
  {
    std::vector<double> col(fNTicks, 0.);
    std::vector<double> ind(fNTicks, 0.);

    unsigned int mxbin = TMath::Min(fNTicks, (int)fNElectResp + fNFieldBins);

    double sumCol = 0.;
    double sumInd = 0.;
  
    for(unsigned int i = 1; i < mxbin; ++i){
      sumCol = 0.;
      sumInd = 0.;
      for(unsigned int j = 0; j < (unsigned int)fNFieldBins; ++j){
	unsigned int k = i - j;
	if(k == 0) break;
	sumCol += fElectResponse[k]*fColFieldResponse[j];
	sumInd += fElectResponse[k]*fIndFieldResponse[j];
      }
      col[i] = sumCol;
      ind[i] = sumInd;
    
    }//end loop over bins;

    ///pad out the rest of the vector with 0.  
    ind.resize(fNTicks, 0.);  
    col.resize(fNTicks, 0.);  

    // write the shapes out to a file
    art::ServiceHandle<art::TFileService> tfs;
    fColTimeShape = tfs->make<TH1D>("ConvolutedCollection",";ticks; Electronics#timesCollection",fNTicks,0,fNTicks);
    fIndTimeShape = tfs->make<TH1D>("ConvolutedInduction",";ticks; Electronics#timesInduction",fNTicks,0,fNTicks);    

    fIndShape.resize(fNTicks/2+1);  
    fColShape.resize(fNTicks/2+1);  

    ///do the FFT of the shapes  
    std::vector<double> delta(fNTicks);  
    delta[0] = 1.0;  
    delta[fNTicks-1]=1.0;  

    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->AlignedSum(ind,delta,false);  
    fFFT->AlignedSum(col,delta,false);  
    fFFT->DoFFT(ind, fIndShape);  
    fFFT->DoFFT(col, fColShape);

    ///check that you did the right thing  
    for(unsigned int i = 0; i < ind.size(); ++i){    
      fColTimeShape->Fill(i, col[i]);    
      fIndTimeShape->Fill(i, ind[i]);  
    }  

    fColTimeShape->Write();
    fIndTimeShape->Write();


    return;
  }

  //-------------------------------------------------
  void SimWire::GenNoise(std::vector<float>& noise, CLHEP::HepRandomEngine& engine)
  {
    CLHEP::RandFlat flat(engine);

    noise.clear();
    noise.resize(fNTicks, 0.);
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);///<noise in frequency space

    double pval = 0.; 
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    //width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);
    for(int i=0; i< fNTicks/2+1; ++i){
      //exponential noise spectrum 
      pval = fNoiseFact*exp(-(double)i*binWidth/fNoiseWidth);
      //low frequency cutoff     
      lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
      //randomize 10%
      flat.fireArray(2,rnd,0,1);
      pval *= lofilter*(0.9+0.2*rnd[0]);
      //random pahse angle
      phase = rnd[1]*2.*TMath::Pi();

      TComplex tc(pval*cos(phase),pval*sin(phase));
      noiseFrequency[i] += tc;
    }
  
    //std::cout << "filled noise freq" << std::endl;

    //inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);

    noiseFrequency.clear();

    ///multiply each noise value by fNTicks as the InvFFT 
    ///divides each bin by fNTicks assuming that a forward FFT
    ///has already been done.
    for(unsigned int i = 0; i < noise.size(); ++i) noise[i] *= 1.*fNTicks;

    return;
  }

  //-------------------------------------------------
  void SimWire::SetFieldResponse()
  {

    //     std::cerr << "SetFieldResponse" << std::endl;

    art::ServiceHandle<geo::Geometry> geo;
 
    double xyz1[3] = {0.};
    double xyz2[3] = {0.};
    double xyzl[3] = {0.};
    ///< should always have at least 2 planes
    geo->Plane(0).LocalToWorld(xyzl, xyz1);
    geo->Plane(1).LocalToWorld(xyzl, xyz2);

    ///this assumes all planes are equidistant from each other,
    ///probably not a bad assumption
    double pitch = xyz2[0] - xyz1[0]; ///in cm

    fColFieldResponse.resize(fNFieldBins, 0.);
    fIndFieldResponse.resize(fNFieldBins, 0.);

    ///set the response for the collection plane first
    ///the first entry is 0

    // write out the response functions to the file
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    fIndFieldResp = tfs->make<TH1D>("InductionFieldResponse",";t (ns);Induction Response",fNTicks,0,fNTicks);
    fColFieldResp = tfs->make<TH1D>("CollectionFieldResponse",";t (ns);Collection Response",fNTicks,0,fNTicks);
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    double driftvelocity=detprop->DriftVelocity(detprop->Efield(),detprop->Temperature())/1000.;  
    int nbinc = TMath::Nint(fCol3DCorrection*(std::abs(pitch))/(driftvelocity*fSampleRate)); ///number of bins //KP
  
    double integral = 0.;
    for(int i = 1; i < nbinc; ++i){
      fColFieldResponse[i] = fColFieldResponse[i-1] + 1.0;
      integral += fColFieldResponse[i];
    }

    for(int i = 0; i < nbinc; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
      fColFieldResp->Fill(i, fColFieldResponse[i]);
    }

    ///now the induction plane
    
    int nbini = TMath::Nint(fInd3DCorrection*(std::abs(pitch))/(driftvelocity*fSampleRate));//KP
    for(int i = 0; i < nbini; ++i){
      fIndFieldResponse[i] = fIndFieldRespAmp/(1.*nbini);
      fIndFieldResponse[nbini+i] = -fIndFieldRespAmp/(1.*nbini);
    
      fIndFieldResp->Fill(i, fIndFieldResponse[i]);
      fIndFieldResp->Fill(nbini+i, fIndFieldResponse[nbini+i]);
    
    }

    fColFieldResp->Write();
    fIndFieldResp->Write();

    return;
  }

  //-------------------------------------------------
  void SimWire::SetElectResponse()
  {

    //     std::cerr << "SetElectResponse" << std::endl;

    art::ServiceHandle<geo::Geometry> geo;
  
    fElectResponse.resize(fNTicks, 0.);
    std::vector<double> time(fNTicks,0.);

    double norm = fShapeTimeConst[1]*TMath::Pi();
    norm /= sin(fShapeTimeConst[1]*TMath::Pi()/fShapeTimeConst[0])/fSampleRate;
    
    double peak = 0.;
    
    for(int i = 0; i < fNTicks; ++i){
      time[i] = (1.*i - 0.33333*fNTicks)*fSampleRate;
      
      // The 120000 is an arbitrary scaling to get displays for microboone
      fElectResponse[i] = 120000.0*exp(-time[i]/fShapeTimeConst[0])/(1. + exp(-time[i]/fShapeTimeConst[1]))/norm;  
      
      if(fElectResponse[i] > peak){
	peak = fElectResponse[i];
      }
    }///end loop over time buckets
    
    ///remove all values of fElectResponse and time where fElectResponse < 0.01*peak 
    peak *= 0.01;
    std::vector<double>::iterator eitr = fElectResponse.begin();
    std::vector<double>::iterator titr = time.begin();
    while(eitr != fElectResponse.end()){
      if(*eitr < peak){
	fElectResponse.erase(eitr);
	time.erase(titr);
	
      }
      else{
	++eitr;
	++titr;
      }
    }//end loop to remove low response values
    
    fNElectResp = fElectResponse.size();

    // write the response out to a file
    art::ServiceHandle<art::TFileService> tfs;
    fElectResp = tfs->make<TH1D>("ElectronicsResponse",";t (ns);Electronics Response",fNElectResp,0,fNElectResp);
    for(unsigned int i = 0; i < fNElectResp; ++i){
      //mf::LogInfo("SimWire") <<"checking ElectResponse: i=  "<< i << "  time[i]=  " << time[i] << "  fElectResponse[i]=  " << fElectResponse[i];
      fElectResp->Fill(i, fElectResponse[i]);
    }

    fElectResp->Write();

    return;
  }

}

namespace detsim{

  DEFINE_ART_MODULE(SimWire)

}
