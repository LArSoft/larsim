////////////////////////////////////////////////////////////////////////
//
// \file LArSimChannelAna_module.cc
//
// \author dmckee@phys.ksu.edu
//
// \brief Build some histograms based on the som::SimChannels created
// \brief by LArVoxelReadout
////////////////////////////////////////////////////////////////////////

// C++ std library includes
#include <vector>
#include <string>
#include <map>

// Root Includes
#include "TH1.h"

// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Simulation/SimChannel.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace larg {

  /// Base class for creation of raw signals on wires. 
  class LArSimChannelAna : public art::EDAnalyzer {
    
  public:
        
    explicit LArSimChannelAna(fhicl::ParameterSet const& pset); 
    virtual ~LArSimChannelAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob(){};
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

    // intilize the histograms
    //
    // Can't be done in Begin job because I want to use LArProperties
    // which used the database, so I test and run on each
    // event. Wasteful and silly, but at least it *works*.
    void ensureHists();

  private:

    std::string            fLArG4ModuleLabel;
    bool                   fDumpChannels; ///< whether to print the content of all channels

    // Flag for initialization done, because we set up histograms the
    // first time through beginRun() so that we can use the
    // database...
    bool initDone;

    TH1D * fChargeXpos;   ///< position of the MC Truth charge deposition 
    TH1D * fChargeYpos;   ///< position of the MC Truth charge deposition 
    TH1D * fChargeZpos;   ///< position of the MC Truth charge deposition 

    TH1D * fTDC;          ///< Which TDCs have activity

    TH1D * fTDCsPerChannel; ///< Number of TDCs with activity
    TH1D * fIDEsPerChannel;

    TH1D * fElectrons;    ///< Electrons in the whole channel entry
    TH1D * fEnergy;       ///< Energy in the whole channel entry

    TH1D * fElectronsPerTDC;
    TH1D * fEnergyPerTDC;

    TH1D * fElectronsPerIDE; 
    TH1D * fEnergyPerIDE;



  }; // class LArSimChannelAna


  //-------------------------------------------------
  LArSimChannelAna::LArSimChannelAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , initDone(false)
    , fChargeXpos()
    , fChargeYpos()
    , fChargeZpos()
    , fTDC()
    , fTDCsPerChannel()
    , fIDEsPerChannel()
    , fElectrons()
    , fEnergy()
    , fElectronsPerTDC()
    , fEnergyPerTDC()
    , fElectronsPerIDE()
    , fEnergyPerIDE()
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  LArSimChannelAna::~LArSimChannelAna()
  {
  }

  void LArSimChannelAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fLArG4ModuleLabel         = p.get< std::string >("LArGeantModuleLabel");
    fDumpChannels             = p.get< bool >("DumpChannels", false);
  }
  //-------------------------------------------------
  void LArSimChannelAna::ensureHists() {
    using CLHEP::cm;
    
    if (initDone) return; // Bail if we've already done this.
    initDone = true; // Insure that we bail later on

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    // geometry data.
    art::ServiceHandle<geo::Geometry> geom; 
    // detector specific properties
    art::ServiceHandle<util::DetectorProperties> detprop;
    // assumes all TPCs are the same
  //  double width  = 2 * geom->TPC(0).HalfWidth(); 
  //  double halfHeight = geom->TPC(0).HalfHeight(); 
  //  double length = geom->TPC(0).Length();
    double boundaries[6];
    geom->AllTPCactiveBoundaries(boundaries);
    double width  = boundaries[1] - boundaries[0];
    double height = boundaries[3] - boundaries[2];
    double length = boundaries[5] - boundaries[4];

    LOG_DEBUG("LArSimChannelAna") << "Plotting a volume of "
      << (width/cm) << " x " << (height/cm) << " x " << (length/cm) << " cm^3";
    
    // Assumes microboone dimensions. Ideally we'd fix this later...
    // charge deposition bins are about, but not exactly, 10 cm each
    fChargeXpos  = tfs->make<TH1D>("hChargeXpos",
                                   "X charge depositions;X (/10cm);Events",
                                   width / (10.*cm), boundaries[0], boundaries[1]);
    fChargeYpos  = tfs->make<TH1D>("hChargeYpos",
                                   "Y charge depositions;Y (/10cm);Events",
                                   height / (10.*cm), boundaries[2], boundaries[3]);
    fChargeZpos  = tfs->make<TH1D>("hChargeZpos",
                                   "Z charge depositions;Z (/10cm);Events",
                                   length / (10.*cm), boundaries[4], boundaries[5]);
    fTDC         = tfs->make<TH1D>("hTDC",
                                   "Active TDC;TDCs;Events;",
                                   detprop->NumberTimeSamples(), 0,
                                   detprop->NumberTimeSamples());
    fTDCsPerChannel =tfs->make<TH1D>("hTDCsPerChannel",
                                     "TDCs per channel entry;# TDCs;Events",
                                     128, 0, detprop->NumberTimeSamples());
    fIDEsPerChannel =tfs->make<TH1D>("hIDEsPerChannel",
                                     "IDE per channel entry;# IDEs;Events",
                                     100,0,20000);
    fElectrons      =tfs->make<TH1D>("hElectrons",
                                     "Electrons per channel;Electrons;Events",
                                     100,0,2e7);
    fEnergy         =tfs->make<TH1D>("hEnergy",
                                     "Energy per channel;energy;Events",
                                     100,0,2500);
    fElectronsPerIDE=tfs->make<TH1D>("hElectronsPerIDE",
                                     "Electrons per IDE;Electrons;Events",
                                     100,0,10000);
    fEnergyPerIDE   =tfs->make<TH1D>("hEnergyPerIDE",
                                     "Energy per IDE;energy;Events",
                                     100,0,50);
    fElectronsPerTDC=tfs->make<TH1D>("hElectronsPerTDC",
                                     "Electrons per TDC;Electrons;Events",
                                     100,0,10000);
    fEnergyPerTDC   =tfs->make<TH1D>("hEnergyPerTDC",
                                     "Energy per YDC;energy;Events",
                                     100,0,50);
    return;

  }

  //-------------------------------------------------
  void LArSimChannelAna::endJob() {}

  //-------------------------------------------------
  void LArSimChannelAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()) {
      throw cet::exception("LArSimChannelAna") << "Not for use on Data yet...\n";
    }
    
    ensureHists();

    art::ServiceHandle<geo::Geometry> geom; 

    art::Handle< std::vector<sim::SimChannel> > chanHandle;
    evt.getByLabel(fLArG4ModuleLabel,chanHandle);
    const std::vector<sim::SimChannel>& scVec(*chanHandle);

    if (fDumpChannels)
      mf::LogInfo("DumpChannels") << "We have " << scVec.size() << " channels";
    
    //++++++++++
    // Loop over the Channels, dump information and fill histograms
    //++++++++++
    unsigned int totalIDEs = 0;
    double totalElectrons = 0;
    double totalEnergy   = 0;
    for (const auto& sc : scVec ) {
      std::map<unsigned short, std::vector<sim::IDE>> tdcidemap=sc.TDCIDEMap();
      fTDCsPerChannel->Fill(tdcidemap.size());
      
      if (fDumpChannels) {
        mf::LogVerbatim("DumpChannels") << "  channel " << sc.Channel()
          << " recorded on " << tdcidemap.size() << " TDC ticks:";
      }

      double channelElectrons = 0.;
      double channelEnergy = 0.;
      
      for (const auto& tdcide : tdcidemap) {
        unsigned int tdc = tdcide.first;
        const std::vector<sim::IDE>& ideVec = tdcide.second;
        totalIDEs += ideVec.size();
        double tdcElectrons=0;
        double tdcEnergy=0;

        fTDC->Fill(tdc);
        
        if (fDumpChannels) {
          mf::LogVerbatim("DumpChannels") << "    TDC #" << tdc
            << " with " << ideVec.size() << " ionization drifting electrons:";
        } // dump TDC info
        
        for (const auto& ide : ideVec) {
          totalElectrons += ide.numElectrons;
          totalEnergy    += ide.energy;
          tdcElectrons   += ide.numElectrons;
          tdcEnergy      += ide.energy;
          
          fChargeXpos->Fill(ide.x);
          fChargeYpos->Fill(ide.y);
          fChargeZpos->Fill(ide.z);
          fElectronsPerIDE->Fill(ide.numElectrons);
          fEnergyPerIDE->Fill(ide.energy);
          
          if (fDumpChannels) {
            mf::LogVerbatim("DumpChannels") << "      " << ide.numElectrons
              << " for " << ide.energy << " MeV at ("
              << ide.x << ", " << ide.y << ", " << ide.z
              << "), Geant track ID:" << ide.trackID;
          } // dump IDE info
          
        } // for IDEs
        fElectronsPerTDC->Fill(tdcElectrons);
        fEnergyPerTDC->Fill(tdcEnergy);
        
        channelElectrons += tdcElectrons;
        channelEnergy += tdcEnergy;
        
        if (fDumpChannels && (ideVec.size() > 1)) {
          mf::LogVerbatim("DumpChannels") << "    => TDC #" << tdc << " has "
            << tdcElectrons << " electrons for " << tdcEnergy << " MeV";
        } // dump TDC summary

      } // for TDCs in channel
      
      if (fDumpChannels) {
        mf::LogVerbatim("DumpChannels") << "  => channel " << sc.Channel()
          << " collected " << channelElectrons << " electrons for "
          << channelEnergy << " MeV";
      } // dump channel summary

    } // for all sim channels
    
    fIDEsPerChannel->Fill(totalIDEs);
    fElectrons->Fill(totalElectrons);
    fEnergy->Fill(totalEnergy);
    
    if (fDumpChannels) {
      mf::LogVerbatim("DumpChannels") << "=> collected " << totalElectrons
        << " electrons for " << totalEnergy << " MeV";
    }
    
  } // end analyze method

  DEFINE_ART_MODULE(LArSimChannelAna)

} // end of hit namespace

