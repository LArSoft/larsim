//////////////////////////////////////////////////////////////////////////
//
// SimWireLBNE10kt class designed to simulate signal on a wire in the TPC
//
// katori@fnal.gov
//
// jti3@fnal.gov
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
/////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Utilities/LArFFT.h"
#include "RawData/raw.h"
#include "Utilities/LArProperties.h"
#include "Utilities/SignalShapingServiceLBNE10kt.h"
#include "Geometry/Geometry.h"

#include "Simulation/sim.h"
#include "Simulation/SimChannel.h"
#include "RawData/RawDigit.h"
#include "Utilities/DetectorProperties.h"

#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

///Detector simulation of raw signals on wires
namespace detsim {

  // Base class for creation of raw signals on wires. 
  class SimWireLBNE10kt : public art::EDProducer {
    
  public:
        
    explicit SimWireLBNE10kt(fhicl::ParameterSet const& pset); 
    virtual ~SimWireLBNE10kt();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void         GenNoise(std::vector<float>& array);

    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    raw::Compress_t        fCompression;      ///< compression type to use
    unsigned int           fNoiseOn;          ///< noise turned on or off for debugging; default is on
    float                  fNoiseFact;        ///< noise scale factor
    float                  fNoiseWidth;       ///< exponential noise width (kHz)
    float                  fLowCutoff;        ///< low frequency filter cutoff (kHz)
    float                  fNoiseFactZ;        ///< noise scale factor for Z (collection) plane
    float                  fNoiseWidthZ;       ///< exponential noise width (kHz)  for Z (collection) plane
    float                  fLowCutoffZ;        ///< low frequency filter cutoff (kHz) for Z (collection) plane
    float                  fNoiseFactU;        ///< noise scale factor  for U plane
    float                  fNoiseWidthU;       ///< exponential noise width (kHz)   for U plane
    float                  fLowCutoffU;        ///< low frequency filter cutoff (kHz)  for U plane
    float                  fNoiseFactV;        ///< noise scale factor   for V plane
    float                  fNoiseWidthV;       ///< exponential noise width (kHz)   for V plane
    float                  fLowCutoffV;        ///< low frequency filter cutoff (kHz)  for V plane
    unsigned int           fZeroThreshold;    ///< Zero suppression threshold
    int                    fNearestNeighbor;         ///< Maximum distance between hits above threshold before they are separated into different blocks
    int                    fNTicks;           ///< number of ticks of the clock
    double                 fSampleRate;       ///< sampling rate in ns
    unsigned int           fNSamplesReadout;  ///< number of ADC readout samples in 1 readout frame
    unsigned int           fNTimeSamples;     ///< number of ADC readout samples in all readout frames (per event)
    unsigned int           fNoiseArrayPoints; ///< number of  points in randomly generated noise array
  
    std::vector<double>    fChargeWork;
    //std::vector< std::vector<float> > fNoise;///< noise on each channel for each time
    std::vector< std::vector<float> > fNoiseZ;///< noise on each channel for each time for Z (collection) plane
    std::vector< std::vector<float> > fNoiseU;///< noise on each channel for each time for U plane
    std::vector< std::vector<float> > fNoiseV;///< noise on each channel for each time for V plane
    
    TH1D*                fNoiseDist;          ///< distribution of noise counts

  }; // class SimWireLBNE10kt

  DEFINE_ART_MODULE(SimWireLBNE10kt)

  //-------------------------------------------------
  SimWireLBNE10kt::SimWireLBNE10kt(fhicl::ParameterSet const& pset)
  {

    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();
    if(fNSamplesReadout != fNTimeSamples ) {
      produces< std::vector<raw::RawDigit>   >("preSpill");
      produces< std::vector<raw::RawDigit>   >("postSpill");
    } 

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;    
    if(compression.Contains("ZeroSuppression",TString::kIgnoreCase)) fCompression = raw::kZeroSuppression;      

    // get the random number seed, use a random default if not specified    
    // in the configuration file.  
    unsigned int seed = pset.get< unsigned int >("Seed", sim::GetRandomNumberSeed());

    createEngine(seed);
  }

  //-------------------------------------------------
  SimWireLBNE10kt::~SimWireLBNE10kt()
  {

    fChargeWork.clear();
 
    for(unsigned int i = 0; i < fNoiseZ.size(); ++i) fNoiseZ[i].clear();
    fNoiseZ.clear();
   
    for(unsigned int i = 0; i < fNoiseU.size(); ++i) fNoiseU[i].clear();
    fNoiseU.clear();
   
    for(unsigned int i = 0; i < fNoiseV.size(); ++i) fNoiseV[i].clear();
    fNoiseV.clear();

  }

  //-------------------------------------------------
  void SimWireLBNE10kt::reconfigure(fhicl::ParameterSet const& p) 
  {
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");


    fNoiseFactZ        = p.get< double              >("NoiseFactZ");
    fNoiseWidthZ       = p.get< double              >("NoiseWidthZ");
    fLowCutoffZ        = p.get< double              >("LowCutoffZ");
    fNoiseFactU        = p.get< double              >("NoiseFactU");
    fNoiseWidthU       = p.get< double              >("NoiseWidthU");
    fLowCutoffU        = p.get< double              >("LowCutoffU");
    fNoiseFactV        = p.get< double              >("NoiseFactV");
    fNoiseWidthV       = p.get< double              >("NoiseWidthV");
    fLowCutoffV        = p.get< double              >("LowCutoffV");
    fZeroThreshold    = p.get< unsigned int        >("ZeroThreshold");
    fNearestNeighbor         = p.get< int                 >("NearestNeighbor");
    fNoiseArrayPoints = p.get< unsigned int        >("NoiseArrayPoints");
    fNoiseOn           = p.get< unsigned int       >("NoiseOn");
    art::ServiceHandle<util::DetectorProperties> detprop;
    fSampleRate       = detprop->SamplingRate();
    fNSamplesReadout  = detprop->ReadOutWindowSize();
    fNTimeSamples  = detprop->NumberTimeSamples();
    
    return;
  }

  //-------------------------------------------------
  void SimWireLBNE10kt::beginJob() 
  { 

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist  = tfs->make<TH1D>("Noise", ";Noise  (ADC);", 1000,   -10., 10.);

    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();

    if ( fNTicks%2 != 0 ) 
      LOG_DEBUG("SimWireLBNE10kt") << "Warning: FFTSize not a power of 2. "
				   << "May cause issues in (de)convolution.\n";

    if ( (int)fNSamplesReadout > fNTicks ) 
      mf::LogError("SimWireLBNE10kt") << "Cannot have number of readout samples "
				      << "greater than FFTSize!";

    fChargeWork.resize(fNTicks, 0.);
    art::ServiceHandle<geo::Geometry> geo;
    
    //Generate noise if selected to be on
    if(fNoiseOn){

      //fNoise.resize(geo->Nchannels());
      fNoiseZ.resize(fNoiseArrayPoints);
      fNoiseU.resize(fNoiseArrayPoints);
      fNoiseV.resize(fNoiseArrayPoints);
      
      // GenNoise() will further resize each channel's 
      // fNoise vector to fNoiseArrayPoints long.
      
      for(unsigned int p = 0; p < fNoiseArrayPoints; ++p){
	fNoiseFact = fNoiseFactZ;
	fNoiseWidth = fNoiseWidthZ;
	fLowCutoff = fLowCutoffZ;

	GenNoise(fNoiseZ[p]);
	for(int i = 0; i < fNTicks; ++i)
	  fNoiseDist->Fill(fNoiseZ[p][i]);
	
	fNoiseFact = fNoiseFactU;
	fNoiseWidth = fNoiseWidthU;
	fLowCutoff = fLowCutoffU;

	GenNoise(fNoiseU[p]);
	for(int i = 0; i < fNTicks; ++i)	 
	  fNoiseDist->Fill(fNoiseU[p][i]);


	fNoiseFact = fNoiseFactV;
	fNoiseWidth = fNoiseWidthV;
	fLowCutoff = fLowCutoffV;
 
    
	GenNoise(fNoiseV[p]);
	for(int i = 0; i < fNTicks; ++i)
	  fNoiseDist->Fill(fNoiseV[p][i]);	
	
      }// end loop over wires
    } 
    return;

  }

  //-------------------------------------------------
  void SimWireLBNE10kt::endJob() 
  {
  }

  //-------------------------------------------------
  void SimWireLBNE10kt::produce(art::Event& evt)
  {
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    unsigned int signalSize = fNTicks;

    // vectors for working
    std::vector<short>    adcvec(signalSize, 0);	
    std::vector<short>    adcvecPreSpill(signalSize, 0);	
    std::vector<short>    adcvecPostSpill(signalSize, 0);	
    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);

    //Get fIndShape and fColShape from SignalShapingService, on the fly
    art::ServiceHandle<util::SignalShapingServiceLBNE10kt> sss;

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
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcolPreSpill(new std::vector<raw::RawDigit>);
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcolPostSpill(new std::vector<raw::RawDigit>);

	  
    unsigned int chan = 0; 
    fChargeWork.clear();
    fChargeWork.resize(fNTicks, 0.);
	  
    std::vector<double> fChargeWorkPreSpill, fChargeWorkPostSpill;

    art::ServiceHandle<util::LArFFT> fFFT;

    // Add all channels  
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    std::map<int,double>::iterator mapIter;      

    bool prepost = false;  // whether or not to do the prespill and postspill digitization
    if(fNSamplesReadout != fNTimeSamples ) { prepost = true; }  

    for(chan = 0; chan < geo->Nchannels(); chan++) {    
      
      fChargeWork.clear();    
      //      fChargeWork.resize(fNTicks, 0.);    
      fChargeWork.resize(fNTimeSamples, 0.);    


      if (prepost) {
        fChargeWorkPreSpill.clear();
        fChargeWorkPostSpill.clear();
        fChargeWorkPreSpill.resize(fNTicks,0);
        fChargeWorkPostSpill.resize(fNTicks,0);
      }

      // get the sim::SimChannel for this channel
      const sim::SimChannel* sc = channels[chan];

      const geo::View_t view = geo->View(chan);

      if( sc ){      

	// loop over the tdcs and grab the number of electrons for each
	for(size_t t = 0; t < fChargeWork.size(); ++t) 
	  fChargeWork[t] = sc->Charge(t);      


        // Convolve charge with appropriate response function 
	if(prepost) {
	  fChargeWorkPreSpill.assign(fChargeWork.begin(), fChargeWork.begin()+fNSamplesReadout);
	  fChargeWorkPostSpill.assign(fChargeWork.begin()+2*fNSamplesReadout, fChargeWork.end());

	  fChargeWork.erase( fChargeWork.begin()+2*fNSamplesReadout, fChargeWork.end() );
	  fChargeWork.erase( fChargeWork.begin(), fChargeWork.begin()+fNSamplesReadout );

	  fChargeWorkPreSpill.resize(fNTicks,0);
	  fChargeWorkPostSpill.resize(fNTicks,0);
	  sss->Convolute(chan, fChargeWorkPreSpill);
	  sss->Convolute(chan, fChargeWorkPostSpill);
	}
	fChargeWork.resize(fNTicks,0);
	sss->Convolute(chan,fChargeWork);


      }

      // noise was already generated for each wire in the event
      // raw digit vec is already in channel order
      // pick a new "noise channel" for every channel  - this makes sure    
      // the noise has the right coherent characteristics to be on one channel

    
      int noisechan = nearbyint(flat.fire()*(1.*(fNoiseArrayPoints-1)+0.1));
      int noisechanpre = nearbyint(flat.fire()*(1.*(fNoiseArrayPoints-1)+0.1));
      int noisechanpost = nearbyint(flat.fire()*(1.*(fNoiseArrayPoints-1)+0.1));
    
      // optimize for speed -- access vectors as arrays 

      double *fChargeWork_a=0;
      double *fChargeWorkPreSpill_a=0;
      double *fChargeWorkPostSpill_a=0;
      short *adcvec_a=0;
      short *adcvecPreSpill_a=0;
      short *adcvecPostSpill_a=0;
      float *noise_a_U=0;
      float *noise_a_V=0;
      float *noise_a_Z=0;
      float *noise_a_Upre=0;
      float *noise_a_Vpre=0;
      float *noise_a_Zpre=0;
      float *noise_a_Upost=0;
      float *noise_a_Vpost=0;
      float *noise_a_Zpost=0;

      if (signalSize>0)	{
	fChargeWork_a = fChargeWork.data();
	adcvec_a = adcvec.data();
	if (prepost) {
	  fChargeWorkPreSpill_a = fChargeWorkPreSpill.data();
	  fChargeWorkPostSpill_a = fChargeWorkPostSpill.data();
	  adcvecPreSpill_a = adcvecPreSpill.data();
	  adcvecPostSpill_a = adcvecPostSpill.data();
	}
	if (fNoiseOn) {
          noise_a_U=(fNoiseU[noisechan]).data();
	  noise_a_V=(fNoiseV[noisechan]).data();
	  noise_a_Z=(fNoiseZ[noisechan]).data();
	  if (prepost) {
	    noise_a_Upre=(fNoiseU[noisechanpre]).data();
	    noise_a_Vpre=(fNoiseV[noisechanpre]).data();
	    noise_a_Zpre=(fNoiseZ[noisechanpre]).data();
	    noise_a_Upost=(fNoiseU[noisechanpost]).data();
	    noise_a_Vpost=(fNoiseV[noisechanpost]).data();
	    noise_a_Zpost=(fNoiseZ[noisechanpost]).data();
	  }
	}
      }

      float tmpfv=0;  // this is here so we do our own rounding from floats to short ints (saves CPU time)
      float tnoise=0;
      float tnoisepre=0;
      float tnoisepost=0;

      if (view != geo::kU && view != geo::kV && view != geo::kZ) {
	mf::LogError("SimWireLBNE10kt") << "ERROR: CHANNEL NUMBER " << chan << " OUTSIDE OF PLANE";
      }

      if(fNoiseOn) {	      
	for(unsigned int i = 0; i < signalSize; ++i){
	  if(view==geo::kU)       { tnoise = noise_a_U[i]; }
	  else if (view==geo::kV) { tnoise = noise_a_V[i]; }
	  else                    { tnoise = noise_a_Z[i]; }
	  tmpfv = tnoise + fChargeWork_a[i];
	  adcvec_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	}
	if (prepost) {
	  for(unsigned int i = 0; i < signalSize; ++i){
	    if(view==geo::kU)      { tnoisepre = noise_a_Upre[i]; tnoisepost = noise_a_Upost[i];  }
	    else if(view==geo::kV) { tnoisepre = noise_a_Vpre[i]; tnoisepost = noise_a_Vpost[i]; }
	    else                   { tnoisepre = noise_a_Zpre[i]; tnoisepost = noise_a_Zpost[i]; }

	    tmpfv = tnoisepre + fChargeWorkPreSpill_a[i];
	    adcvecPreSpill_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	    tmpfv = tnoisepost + fChargeWorkPostSpill_a[i];
	    adcvecPostSpill_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	  }
	}
      }
      else {   // no noise, so just round the values to nearest short ints and store them
	for(unsigned int i = 0; i < signalSize; ++i){
	  tmpfv = fChargeWork_a[i];
	  adcvec_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	}
	if (prepost) {
	  for(unsigned int i = 0; i < signalSize; ++i){
	    tmpfv = fChargeWorkPreSpill_a[i];
	    adcvecPreSpill_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	    tmpfv = fChargeWorkPostSpill_a[i];
	    adcvecPostSpill_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	  }
	}
      }

      // resize the adcvec to be the correct number of time samples, 
      // just drop the extra samples

      // compress the adc vector using the desired compression scheme,
      // if raw::kNone is selected nothing happens to adcvec
      // This shrinks adcvec, if fCompression is not kNone.

      //std::cout << "Channel view: " << view << std::endl;

      adcvec.resize(fNSamplesReadout);
      raw::Compress(adcvec, fCompression, fZeroThreshold, fNearestNeighbor); 
      raw::RawDigit rd(chan, fNSamplesReadout, adcvec, fCompression);
      adcvec.resize(signalSize);        // Then, resize adcvec back to full length.  Do not initialize to zero (slow)
      digcol->push_back(rd);            // add this digit to the collection

      // do all this for the prespill and postspill samples if need be
      if (prepost) {
        adcvecPreSpill.resize(fNSamplesReadout);
        adcvecPostSpill.resize(fNSamplesReadout);
        raw::Compress(adcvecPreSpill, fCompression, fZeroThreshold, fNearestNeighbor); 
        raw::Compress(adcvecPostSpill, fCompression, fZeroThreshold, fNearestNeighbor); 
        raw::RawDigit rdPreSpill(chan, fNSamplesReadout, adcvecPreSpill, fCompression);
        raw::RawDigit rdPostSpill(chan, fNSamplesReadout, adcvecPostSpill, fCompression);
        adcvecPreSpill.resize(signalSize);
        adcvecPostSpill.resize(signalSize);
        digcolPreSpill->push_back(rdPreSpill);
        digcolPostSpill->push_back(rdPostSpill);
      }

    }// end loop over channels      

    evt.put(std::move(digcol));
    if(prepost) {
      evt.put(std::move(digcolPreSpill), "preSpill");
      evt.put(std::move(digcolPostSpill), "postSpill");
    }

    return;
  }

  //-------------------------------------------------
  void SimWireLBNE10kt::GenNoise(std::vector<float>& noise)
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    noise.clear();
    noise.resize(fNTicks,0.0);
    // noise in frequency space
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);

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
  
    // std::cout << "filled noise freq" << std::endl;

    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);

    noiseFrequency.clear();

    // multiply each noise value by fNTicks as the InvFFT 
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    for(unsigned int i = 0; i < noise.size(); ++i) noise[i] *= 1.*fNTicks;

    return;
  }
}
