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
#include <string>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Root Includes
#include "TH1.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

///Detector simulation of raw signals on wires
namespace larg {

  /// Base class for creation of raw signals on wires.
  class LArSimChannelAna : public art::EDAnalyzer {

  public:

    explicit LArSimChannelAna(fhicl::ParameterSet const& pset);

    /// read/write access to event
    void analyze (const art::Event& evt);

    // intilize the histograms
    //
    // Can't be done in Begin job because I want to use LArProperties
    // which used the database, so I test and run on each
    // event. Wasteful and silly, but at least it *works*.
    void ensureHists();

  private:

    std::string            fLArG4ModuleLabel;

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
    , fLArG4ModuleLabel{pset.get< std::string >("LArGeantModuleLabel")}
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
  {}

  //-------------------------------------------------
  void LArSimChannelAna::ensureHists() {
    if (initDone) return; // Bail if we've already done this.
    initDone = true; // Insure that we bail later on

    // get access to the TFile service
    art::ServiceHandle<art::TFileService const> tfs;
    // geometry data.
    art::ServiceHandle<geo::Geometry const> geom;
    // detector specific properties
    const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // assumes all TPCs are the same
    double width  = 2 * geom->TPC(0).HalfWidth();
    double halfHeight = geom->TPC(0).HalfHeight();
    double length = geom->TPC(0).Length();

    // Assumes microboone dimensions. Ideally we'd fix this later...
    fChargeXpos  = tfs->make<TH1D>("hChargeXpos",
				   "X charge depositions;X (cm);Events",
				   101, 0.0, width);
    fChargeYpos  = tfs->make<TH1D>("hChargeYpos",
				   "Y charge depositions;Y (cm);Events",
				   101, -halfHeight, halfHeight);
    fChargeZpos  = tfs->make<TH1D>("hChargeZpos",
				   "Z charge depositions;Z (cm);Events",
				   101, 0.0, length);
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
  void LArSimChannelAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()) {
      throw cet::exception("LArSimChannelAna") << "Not for use on Data yet...\n";
    }

    ensureHists();

    art::ServiceHandle<geo::Geometry const> geom;

    art::Handle< std::vector<sim::SimChannel> > chanHandle;
    evt.getByLabel(fLArG4ModuleLabel,chanHandle);
    const std::vector<sim::SimChannel>& scVec(*chanHandle);

    //++++++++++
    // Loop over the Chnnels and fill histograms
    //++++++++++
    unsigned int totalIDEs = 0;
    double totalElectrons = 0;
    double totalEnergy   = 0;
    for (const auto& sc : scVec ) {
      const auto & tdcidemap=sc.TDCIDEMap();
      fTDCsPerChannel->Fill(tdcidemap.size());

      for (const auto& tdcide : tdcidemap) {
	unsigned int tdc = tdcide.first;
	const std::vector<sim::IDE>& ideVec = tdcide.second;
	totalIDEs += ideVec.size();
	double tdcElectrons=0;
	double tdcEnergy=0;

	fTDC->Fill(tdc);

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
	}
	fElectronsPerTDC->Fill(tdcElectrons);
	fEnergyPerTDC->Fill(tdcEnergy);
      }
    }
    fIDEsPerChannel->Fill(totalIDEs);
    fElectrons->Fill(totalElectrons);
    fEnergy->Fill(totalEnergy);
    return;
  }//end analyze method

  DEFINE_ART_MODULE(LArSimChannelAna)

} // end of hit namespace
