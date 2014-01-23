
#ifndef TrigAnaMicroBoone_CC
#define TrigAnaMicroBoone_CC

// LArSoft includes
#include "Geometry/Geometry.h"
#include "TriggerAlgoMicroBoone.h"
#include "OpticalDetectorData/FIFOChannel.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT include
#include <TH1D.h>

namespace trigana {
 
  class TrigAnaMicroBoone : public art::EDAnalyzer{
  public:
 
    TrigAnaMicroBoone(const fhicl::ParameterSet&);
    virtual ~TrigAnaMicroBoone();

    void beginJob();

    void analyze (const art::Event&); 

  private:
    
    void add_fifo_waveform(optdata::TimeSlice_t time, 
			   std::vector<optdata::ADC_Count_t> &output,
			   const optdata::FIFOChannel input);

    art::ServiceHandle<trigger::TriggerAlgoMicroBoone> _trig_mod;

    std::string _name_opt_module;

  };

} 

namespace trigana {
  DEFINE_ART_MODULE(TrigAnaMicroBoone)
}


namespace trigana {

  //-----------------------------------------------------------------------
  // Constructor
  TrigAnaMicroBoone::TrigAnaMicroBoone(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    _name_opt_module   = pset.get< std::string >("ModuleName_PMT");
    
  }

  //-----------------------------------------------------------------------
  // Destructor
  TrigAnaMicroBoone::~TrigAnaMicroBoone(){}
   
  //-----------------------------------------------------------------------
  void TrigAnaMicroBoone::beginJob(){}
   

  //-----------------------------------------------------------------------
  void TrigAnaMicroBoone::analyze(const art::Event& event) 
  {

    //
    // PMT FEM modules chops a big waveform into multiple small readout windows.
    // Here we patch chopped-waveforms that are triggered for readout, and recover a (part of) big waveform.
    //
    // Strategy is the following:
    //
    // (1) run trigger simulation, obtain a list of waveforms triggered for readout
    // (2) prepare waveform container per channel per gain (i.e. high/low)
    // (3) loop over the readout waveforms and store in the container.

    art::ServiceHandle<art::TFileService> fileService;

    //
    // (1) run trigger sim., get a list of waveforms triggered for readout
    //
    // Use trigger algorithm service module to apply trigger logic on data
    _trig_mod->RunTriggerSim(event);

    // Retrieve a set of indexes that tells us which input data is triggered for readout
    const std::set<size_t>* index_fifo = _trig_mod->GetIndexes_FIFO();

    //
    // (2) prepare waveform container per channel per gain.
    //     I think having all channels with the same length makes sense though that takes 
    //     up memory depending on the simulation config (readout window size in particular).
    //     But this is an example code, so why not prioritize the clarity? Let's use same 
    //     length for all channels.

    //     
    //     For simplicity let me first look up how long these vector containers should be.
    //
    size_t max_waveform_length = 0;
    
    // Make an attempt to obtain PMT data ... FIFOChannel
    std::vector<const optdata::FIFOChannel*> pmtFIFOArray;
    try {
      event.getView(_name_opt_module, pmtFIFOArray);
    }
    catch ( art::Exception const& e ) {
      if (e.categoryCode() != art::errors::ProductNotFound) {
	throw;
      }
    }

    // Loop over readout FIFO channel array
    for(auto index : *index_fifo){
      
      const optdata::FIFOChannel* fifo(pmtFIFOArray.at(index));

      if( max_waveform_length < (fifo->size() + (size_t)(fifo->TimeSlice())) )

	max_waveform_length = (fifo->size() + (size_t)(fifo->TimeSlice()) );

    }
    mf::LogInfo("TrigAnaMicroBoone")<<"Waveform length: " <<max_waveform_length;
    
    // Define waveform containers ... high and low gain
    typedef std::map<optdata::Channel_t,std::vector<optdata::ADC_Count_t> > Patched_Waveform_t;
    Patched_Waveform_t highWaveform;
    Patched_Waveform_t lowWaveform;

    //
    // (3) Store readout waveforms into the container
    //
    for(auto index : *index_fifo){

      const optdata::FIFOChannel* fifo(pmtFIFOArray.at(index));

      switch(fifo->Category()){
	
      case optdata::kUndefined:
      case optdata::kBeamPMTTrigger:
      case optdata::kCosmicPMTTrigger:
	break;
      case optdata::kHighGain:
      case optdata::kFEMCosmicHighGain:
      case optdata::kFEMBeamHighGain:

	if(highWaveform.find(fifo->ChannelNumber()) == highWaveform.end())

	  highWaveform[fifo->ChannelNumber()] = std::vector<optdata::ADC_Count_t>(max_waveform_length,0);
	    
	add_fifo_waveform(fifo->TimeSlice(),
			  ((*(highWaveform.find(fifo->ChannelNumber()))).second),
			  *fifo);
	break;
      case optdata::kLowGain:
      case optdata::kFEMCosmicLowGain:
      case optdata::kFEMBeamLowGain:

	if(lowWaveform.find(fifo->ChannelNumber()) == lowWaveform.end())

	  lowWaveform[fifo->ChannelNumber()] = std::vector<optdata::ADC_Count_t>(max_waveform_length,0);

	add_fifo_waveform(fifo->TimeSlice(),
			  ((*(lowWaveform.find(fifo->ChannelNumber()))).second),
			  *fifo);
	break;
      }

    }
    
    
    for(Patched_Waveform_t::const_iterator iter(highWaveform.begin());
	iter != highWaveform.end();
	++iter) {
      
      TH1D* h=fileService->make<TH1D>(Form("hWF_High_Ch%02d",(*iter).first),
				      Form("Waveform for Ch. %-2d", (*iter).first),
				      (*iter).second.size(),-0.5,(*iter).second.size() - 0.5);
      
      for(size_t i=0; i < (*iter).second.size(); i++)
	
	h->SetBinContent(i+1,(*iter).second[i]);

    }

    for(Patched_Waveform_t::const_iterator iter(lowWaveform.begin());
	iter != lowWaveform.end();
	++iter) {
      
      TH1D* h=fileService->make<TH1D>(Form("hWF_Low_Ch%02d",(*iter).first),
				      Form("Waveform for Ch. %-2d", (*iter).first),
				      (*iter).second.size(),-0.5,(*iter).second.size() - 0.5);
      
      for(size_t i=0; i < (*iter).second.size(); i++)
	
	h->SetBinContent(i+1,(*iter).second[i]);
      
    }



  }

  void TrigAnaMicroBoone::add_fifo_waveform(optdata::TimeSlice_t time, 
					    std::vector<optdata::ADC_Count_t> &output,
					    const optdata::FIFOChannel input)
  {

    if( (time + input.size()) > output.size() )

      output.resize( (time + input.size()) );
    
    for(size_t index=0; index<input.size(); ++index)

      output[index + time] += input.at(index);

  }
					    

} // namespace opdet


#endif
