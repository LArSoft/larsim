/**
 * @file SimDriftElectrons_module.cxx
 *
 * @brief Transports energy depositions in the LAr TPC to the TPC
 * channels.
 *
 * author: 
 * This module was prepared by William Seligman (me), based on code
 * that had been in
 * `LArG4::LArVoxelReadout::DriftIonizationElectrons`. However, though
 * I wrote the original LArVoxelReadout code, I have no idea who added
 * DriftIonizationElectrons. I probably will not be able to answer any
 * questions about how this code works.
 * 
 * This module acts on sim::SimEnergyDeposit, the single energy
 * depositions from the detector simulation (LArG4), and simulates the
 * transport of the ensuing ionization electrons to the readout
 * channels:
 * 
 * 1. the number of ionisation electrons is read from the current
 *   `larg4::IonizationAndScintillation` instance
 * 2. space charge displacement is optionally applied
 * 3. lifetime correction is applied
 * 4. charge is split in small electron clusters
 * 5. each cluster is subject to longitudinal and transverse diffusion
 * 6. each cluster is assigned to one TPC channel for each wire plane
 * 7. optionally, charge is forced to stay on the planes; otherwise charge
 *    drifting outside the plane is lost
 * 
 * For each energy deposition, entries on the appropriate
 * `sim::SimChannel` are added, with the information of the position
 * where the energy deposit happened (in global coordinates,
 * centimeters), the ID of the track in the detector simulation which
 * produced the deposition, and the quantized time of arrival to the
 * channel (in global TDC tick units).  At most one entry is added for
 * each electron cluster, but entries from the same energy deposit can
 * be compacted if falling on the same TDC tick.
 * 
 * Options
 * --------
 * 
 * A few optional behaviours are supported:
 * 
 * * lead off-plane charge to the planes: regulated by
 *   `RecoverOffPlaneDeposit()`, if charge which reaches a wire plane
 *   is actually off it by less than the chosen margin, it's accounted for by
 *   that plane; by default the margin is 0 and all the charge off the plane
 *   is lost (with a warning)
 * 
 */

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/RandomUtils/NuRandomService.h"

// External libraries
#include "CLHEP/Random/RandGauss.h"
#include "TMath.h"

// C++ includes
#include <vector>
#include <map>
#include <algorithm> // std::find
#include <cmath>

//stuff from wes
#include "larsim/IonizationScintillation/ISCalculationSeparate.h"


namespace detsim {

  // Base class for creation of raw signals on wires. 
  class SimDriftElectrons : public art::EDProducer {
    
  public:
        
    explicit SimDriftElectrons(fhicl::ParameterSet const& pset);
    virtual ~SimDriftElectrons() {};
    
    // Methods that that are available for a module derived from
    // art::EDProducer.
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    // The label of the module that created the sim::SimEnergyDeposit
    // objects (as of Oct-2017, this is probably "largeant"). 
    art::InputTag fSimModuleLabel; 

    const detinfo::DetectorClocks* fTimeService;
    std::unique_ptr<CLHEP::RandGauss> fRandGauss;

    double fElectronLifetime;
    double fElectronClusterSize;
    int    fMinNumberOfElCluster;
    double fSampleRate;
    int    fTriggerOffset;
    double fLongitudinalDiffusion;
    double fTransverseDiffusion;

    double fLifetimeCorr_const;
    double fLDiff_const;
    double fTDiff_const;
    double fRecipDriftVel[3];
    //double fOffPlaneMargin;
    
    // In order to create the associations, for each channel we create
    // we have to keep track of its index in the output vector, and the 
    // indexes of all the steps that contributed to it. 
    typedef struct {
      size_t              channelIndex;
      std::vector<size_t> stepList;
    } ChannelBookKeeping_t;

    // Define type: channel -> sim::SimChannel's bookkeeping. 
    typedef std::map<raw::ChannelID_t, ChannelBookKeeping_t> ChannelMap_t;

    // Array of maps of channel data indexed by [cryostat,tpc]
    std::vector< std::vector<ChannelMap_t> > fChannelMaps; 
    // The above ensemble may be thought of as a 3D array of
    // ChannelBookKeepings: e.g., SimChannel[cryostat,tpc,channel ID].

    // Save the number of cryostats, and the number of TPCs within
    // each cryostat.
    size_t fNCryostats;
    std::vector<size_t> fNTPCs;

    // Per-cluster information. 
    std::vector< double > fXDiff;
    std::vector< double > fYDiff;
    std::vector< double > fZDiff;
    std::vector< double > fnElDiff;
    std::vector< double > fnEnDiff;

    double fDriftClusterPos[3];

    // Utility routine.
    //geo::vect::Vector_t RecoverOffPlaneDeposit (geo::vect::Vector_t const&, geo::PlaneGeo const&) const;

    art::ServiceHandle<geo::Geometry> fGeometry;  ///< Handle to the Geometry service
    ::detinfo::ElecClock              fClock;     ///< TPC electronics clock


    //IS calculationg
    larg4::ISCalculationSeparate fISAlg;

  }; // class SimDriftElectrons

  //-------------------------------------------------
  SimDriftElectrons::SimDriftElectrons(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<sim::SimChannel> >();
    //produces< art::Assns<sim::SimChannel, sim::SimEnergyDeposit> >();
    
    // create a default random engine; obtain the random seed from
    // NuRandomService, unless overridden in configuration with key
    // "Seed"
    art::ServiceHandle<rndm::NuRandomService>()
      ->createEngine(*this, pset, "Seed");
    
    /**
     * @brief Sets the margin for recovery of charge drifted off-plane.
     * @param margin the extent of the margin on each frame coordinate [cm]
     *
     * This method sets the margin for the recovery of off-plane ionization
     * charge. See `RecoverOffPlaneDeposit()` for a description of that feature.
     *
     */
    //fOffPlaneMargin = pset.get< double >("ChargeRecoveryMargin",0.0);
    // Protection against a silly value.
    //fOffPlaneMargin = std::max(fOffPlaneMargin,0.0);
  }

  //-------------------------------------------------
  void SimDriftElectrons::reconfigure(fhicl::ParameterSet const& p) 
  {
    std::string label= p.get< std::string >("SimulationLabel");
    fSimModuleLabel = art::InputTag(label);

    return;
  }

  //-------------------------------------------------
  void SimDriftElectrons::beginJob() 
  { 
    fTimeService = lar::providerFrom<detinfo::DetectorClocksService>();
    fClock = fTimeService->TPCClock();
    
    // Set up the gaussian generator. 
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine();
    fRandGauss = std::unique_ptr<CLHEP::RandGauss>(new CLHEP::RandGauss(engine));

    // Define the physical constants we'll use.

    auto const * detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    fElectronLifetime      = detprop->ElectronLifetime();
    for (int i = 0; i<3; ++i) {
      double driftVelocity = detprop->DriftVelocity(detprop->Efield(i),
						    detprop->Temperature())/1000.;
      fRecipDriftVel[i] = 1./driftVelocity;
    }      

    // To-do: Move the parameters we fetch from "LArG4" to detector
    // properties.
    art::ServiceHandle<sim::LArG4Parameters> paramHandle;
    fElectronClusterSize   = paramHandle->ElectronClusterSize();
    fMinNumberOfElCluster  = paramHandle->MinNumberOfElCluster();
    fLongitudinalDiffusion = paramHandle->LongitudinalDiffusion();
    fTransverseDiffusion   = paramHandle->TransverseDiffusion();

    LOG_DEBUG("SimDriftElectrons")  << " e lifetime: "        << fElectronLifetime
				    << "\n Temperature: "     << detprop->Temperature()
				    << "\n Drift velocity: "  << 1./fRecipDriftVel[0]
				    <<" "<<1./fRecipDriftVel[1]<<" "<<1./fRecipDriftVel[2];

    fLifetimeCorr_const = -1000. * fElectronLifetime;
    fLDiff_const        = std::sqrt(2.*fLongitudinalDiffusion);
    fTDiff_const        = std::sqrt(2.*fTransverseDiffusion);

    // For this detector's geometry, save the number of cryostats and
    // the number of TPCs within each cryostat.
    fNCryostats = fGeometry->Ncryostats();
    fNTPCs.resize(fNCryostats);
    for ( size_t n = 0; n < fNCryostats; ++n )
      fNTPCs[n] = fGeometry->NTPC(n);


    fISAlg.Initialize(lar::providerFrom<detinfo::LArPropertiesService>(),
		      detprop,
		      &(*paramHandle),
		      lar::providerFrom<spacecharge::SpaceChargeService>());


    return;
  }

  //-------------------------------------------------
  void SimDriftElectrons::endJob() 
  {
  }

  //-------------------------------------------------
  void SimDriftElectrons::produce(art::Event& event)
  {
    // Fetch the SimEnergyDeposit objects for this event.
    typedef art::Handle< std::vector<sim::SimEnergyDeposit> > energyDepositHandle_t;
    energyDepositHandle_t energyDepositHandle;
    // If there aren't any energy deposits for this event, don't
    // panic. It's possible someone is doing a study with events
    // outside the TPC, or where there are only non-ionizing
    // particles, or something like that.
    if (!event.getByLabel(fSimModuleLabel, energyDepositHandle)) 
      return;

    // Define the container for the SimChannel objects that will be
    // transferred to the art::Event after the put statement below.
    std::unique_ptr< std::vector<sim::SimChannel> > channels(new std::vector<sim::SimChannel>);

    // Clear the channel maps from the last event. Remember,
    // fChannelMaps is an array[cryo][tpc] of maps.
    size_t cryo = 0;
    fChannelMaps.resize(fNCryostats);
    for (auto& cryoData: fChannelMaps) { // each, a vector of maps
      cryoData.resize(fNTPCs[cryo++]);
      for (auto& channelsMap: cryoData) channelsMap.clear(); // each, a map
    }

    // We're going through the input vector by index, rather than by
    // iterator, because we need the index number to compute the
    // associations near the end of this method.
    auto const& energyDeposits = *energyDepositHandle;
    auto energyDepositsSize = energyDeposits.size();

    // For each energy deposit in this event
    for ( size_t edIndex = 0; edIndex < energyDepositsSize; ++edIndex )
      {
	auto const& energyDeposit = energyDeposits[edIndex];

	// "xyz" is the position of the energy deposit in world
	// coordinates. Note that the units of distance in
	// sim::SimEnergyDeposit are supposed to be cm.
	auto const mp = energyDeposit.MidPoint();
	double const xyz[3] = { mp.X(), mp.Y(), mp.Z() };

	// From the position in world coordinates, determine the
	// cryostat and tpc. If somehow the step is outside a tpc
	// (e.g., cosmic rays in rock) just move on to the next one.
	unsigned int cryostat = 0;
	try {
	  fGeometry->PositionToCryostat(xyz, cryostat);
	}
	catch(cet::exception &e){
	  mf::LogWarning("SimDriftElectrons") << "step "// << energyDeposit << "\n"
	  				      << "cannot be found in a cryostat\n"
	  				      << e;
	  continue;
	}

	unsigned int tpc = 0;
	try {
	  fGeometry->PositionToTPC(xyz, tpc, cryostat);
	}
	catch(cet::exception &e){
	  mf::LogWarning("SimDriftElectrons") << "step "// << energyDeposit << "\n"
					      << "cannot be found in a TPC\n"
					      << e;
	  continue;
	}

	const geo::TPCGeo& tpcGeo = fGeometry->TPC(tpc, cryostat);

	// X drift distance - the drift direction can be either in
	// the positive or negative direction, so use std::abs

	/// \todo think about effects of drift between planes 
	double XDrift = std::abs(xyz[0] - tpcGeo.PlaneLocation(0)[0]);
	//std::cout<<tpcGeo.DriftDirection()<<std::endl;
	if (tpcGeo.DriftDirection() == geo::kNegX)
	  XDrift = xyz[0] - tpcGeo.PlaneLocation(0)[0];
	else if (tpcGeo.DriftDirection() == geo::kPosX)
	  XDrift = tpcGeo.PlaneLocation(0)[0] - xyz[0];
      
	if(XDrift < 0.) continue;

	// Space-charge effect (SCE): Get SCE {x,y,z} offsets for
	// particular location in TPC
	std::vector<double> posOffsets{0.0,0.0,0.0};
	auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
	if (SCE->EnableSimSpatialSCE() == true)
	  {
	    posOffsets = SCE->GetPosOffsets(mp.X(),mp.Y(),mp.Z());
	  }
	posOffsets[0] = (-1*posOffsets[0]);
	XDrift += posOffsets[0];
	  
	// Space charge distortion could push the energy deposit beyond the wire
	// plane (see issue #15131). Given that we don't have any subtlety in the
	// simulation of this region, bringing the deposit exactly on the plane
	// should be enough for the time being.
	if (XDrift < 0.) XDrift = 0.;
      
	// Drift time (nano-sec)
	double TDrift = XDrift * fRecipDriftVel[0];
	if (tpcGeo.Nplanes() == 2){// special case for ArgoNeuT (plane 0 is the second wire plane)
	  TDrift = ((XDrift - tpcGeo.PlanePitch(0,1)) * fRecipDriftVel[0] 
		    + tpcGeo.PlanePitch(0,1) * fRecipDriftVel[1]);
	}
          
	auto efieldoffsets = SCE->GetEfieldOffsets(mp.X(),mp.Y(),mp.Z());
	fISAlg.Reset();
	fISAlg.CalculateIonizationAndScintillation(energyDeposit,
						   efieldoffsets);

	//std::cout << "Got " << fISAlg.NumberIonizationElectrons() << "." << std::endl;

	const double lifetimecorrection = TMath::Exp(TDrift / fLifetimeCorr_const);
	const int    nIonizedElectrons  = fISAlg.NumberIonizationElectrons();
	const double energy             = energyDeposit.Energy();
      
	// if we have no electrons (too small energy or too large recombination)
	// we are done already here
	if (nIonizedElectrons <= 0) {
	  LOG_DEBUG("SimDriftElectrons")
	    << "step "// << energyDeposit << "\n"
	    << "No electrons drifted to readout, " << energy << " MeV lost.";
	  continue;
	}

	// includes the effect of lifetime
	const double nElectrons = nIonizedElectrons * lifetimecorrection;
	//std::cout << "After lifetime, " << nElectrons << " electrons." << std::endl;

	// Longitudinal & transverse diffusion sigma (cm)
	double SqrtT    = std::sqrt(TDrift);
	double LDiffSig = SqrtT * fLDiff_const;
	double TDiffSig = SqrtT * fTDiff_const;
	double electronclsize = fElectronClusterSize;
      
	// Number of electron clusters.
	int nClus = (int) std::ceil(nElectrons / electronclsize);
	if (nClus < fMinNumberOfElCluster)
	  {
	    electronclsize = nElectrons / fMinNumberOfElCluster; 
	    if (electronclsize < 1.0)
	      {
		electronclsize = 1.0;
	      }
	    nClus = (int) std::ceil(nElectrons / electronclsize);
	  }
      
	// Empty and resize the electron-cluster vectors.
	fXDiff.clear();
	fYDiff.clear();
	fZDiff.clear();
	fnElDiff.clear();
	fnEnDiff.clear();
	fXDiff.resize(nClus);
	fYDiff.resize(nClus);
	fZDiff.resize(nClus);
	fnElDiff.resize(nClus, electronclsize);
	fnEnDiff.resize(nClus);

	// fix the number of electrons in the last cluster, that has a smaller size
	fnElDiff.back() = nElectrons - (nClus-1)*electronclsize;
      
	for(size_t xx = 0; xx < fnElDiff.size(); ++xx){
	  if(nElectrons > 0) fnEnDiff[xx] = energy/nElectrons*fnElDiff[xx];
	  else               fnEnDiff[xx] = 0.;
	}

	//std::cout << "Split into, " << nClus << " clusters." << std::endl;
	  
	double const avegageYtransversePos
	  = xyz[1] + posOffsets[1];
	double const avegageZtransversePos
	  = xyz[2] + posOffsets[2];

	// Smear drift times by x position and drift time
	if (LDiffSig > 0.0)
	  fRandGauss->fireArray( nClus, &fXDiff[0], 0., LDiffSig);
	else
	  fXDiff.assign(nClus, 0.0);
      
	if (TDiffSig > 0.0) {
	  // Smear the Y,Z position by the transverse diffusion
	  fRandGauss->fireArray( nClus, &fYDiff[0], avegageYtransversePos, TDiffSig);
	  fRandGauss->fireArray( nClus, &fZDiff[0], avegageZtransversePos, TDiffSig);
	}
	else {
	  fYDiff.assign(nClus, avegageYtransversePos);
	  fZDiff.assign(nClus, avegageZtransversePos);
	}

	//std::cout << "Smeared the " << nClus << " clusters." << std::endl;
	
	// make a collection of electrons for each plane
	for(size_t p = 0; p < tpcGeo.Nplanes(); ++p){
	  
	  //std::cout << "Doing plane " << p << std::endl;

	  geo::PlaneGeo const& plane = tpcGeo.Plane(p);

	  double Plane0Pitch = tpcGeo.Plane0Pitch(p);
        
	  // "-" sign is because Plane0Pitch output is positive. Andrzej
	  fDriftClusterPos[0] = tpcGeo.PlaneLocation(0)[0] - Plane0Pitch;

	  // Drift nClus electron clusters to the induction plane
	  for(int k = 0; k < nClus; ++k){

	    //std::cout << "\tCluser " << k << " diffs are " 
	    //	      << fXDiff[k] << " " << fYDiff[k] << " " << fZDiff[k]
	    //	      << std::endl;

	    // Correct drift time for longitudinal diffusion and plane
	    double TDiff = TDrift + fXDiff[k] * fRecipDriftVel[0];
	    // Take into account different Efields between planes
	    // Also take into account special case for ArgoNeuT where Nplanes = 2.
	    for (size_t ip = 0; ip<p; ++ip){
	      TDiff += tpcGeo.PlanePitch(ip,ip+1) * fRecipDriftVel[tpcGeo.Nplanes()==3?ip+1:ip+2];
	    }
	    fDriftClusterPos[1] = fYDiff[k];
	    fDriftClusterPos[2] = fZDiff[k];
          
	    /// \todo think about effects of drift between planes
          
	    // grab the nearest channel to the fDriftClusterPos position
	    try{
	      /*
	      if (fOffPlaneMargin != 0) {
		// get the effective position where to consider the charge landed;
		// 
		// Some optimisations are possible; in particular, this method
		// could be extended to inform us if the point was too far.
		// Currently, if that is the case the code will proceed, find the
		// point is off plane, emit a warning and skip the deposition.
		//
		auto const landingPos
		  = RecoverOffPlaneDeposit({ fDriftClusterPos[0], fDriftClusterPos[1], fDriftClusterPos[2] }, plane);
		fDriftClusterPos[0] = landingPos.X();
		fDriftClusterPos[1] = landingPos.Y();
		fDriftClusterPos[2] = landingPos.Z();
              
	      } // if charge lands off plane
	      */
	      raw::ChannelID_t channel = fGeometry->NearestChannel(fDriftClusterPos, p, tpc, cryostat);

	      //std::cout << "\tgot channel " << channel << " for cluster " << k << std::endl;
            
	      /// \todo check on what happens if we allow the tdc value to be
	      /// \todo beyond the end of the expected number of ticks
	      // Add potential decay/capture/etc delay effect, simTime.
	      auto const simTime = energyDeposit.Time();
	      unsigned int tdc = fClock.Ticks(fTimeService->G4ToElecTime(TDiff + simTime));
        
	      // Find whether we already have this channel in our map. 
	      ChannelMap_t& channelDataMap = fChannelMaps[cryostat][tpc];
	      auto search = channelDataMap.find(channel);
	      
	      // We will find (or create) the pointer to a
	      // sim::SimChannel.
	      //sim::SimChannel* channelPtr = NULL;
	      size_t channelIndex=0;

	      // Have we created the sim::SimChannel corresponding to
	      // channel ID?
	      if (search == channelDataMap.end())
		{
		  //std::cout << "\tHaven't done this channel before." << std::endl;

		  // We haven't. Initialize the bookkeeping information
		  // for this channel.
		  ChannelBookKeeping_t bookKeeping;
		  
		  // Add a new channel to the end of the list we'll
		  // write out after we've processed this event.
		  bookKeeping.channelIndex = channels->size();
		  channels->emplace_back( channel );
		  channelIndex = bookKeeping.channelIndex;

		  // Save the pointer to the newly-created
		  // sim::SimChannel.
		  //channelPtr = &(channels->back());
		  //bookKeeping.channelPtr = channelPtr;
		  
		  // Initialize a vector with the index of the step that
		  // created this channel.
		  bookKeeping.stepList.push_back( edIndex );
		  
		  // Save the bookkeeping information for this channel.
		  channelDataMap[channel] = bookKeeping;
		}
	      else {
		// We've created this SimChannel for a previous energy
		// deposit. Get its address.

		//std::cout << "\tHave seen this channel before." << std::endl;

		auto& bookKeeping = search->second;
		channelIndex = bookKeeping.channelIndex;
		//channelPtr = bookKeeping.channelPtr;
		
		// Has this step contributed to this channel before?
		auto& stepList = bookKeeping.stepList;
		auto stepSearch = std::find(stepList.begin(), stepList.end(), edIndex );
		if ( stepSearch == stepList.end() ) {
		  // No, so add this step's index to the list.
		  stepList.push_back( edIndex );
		}
	      }
	      sim::SimChannel* channelPtr = &(channels->at(channelIndex));
	      //std::cout << "\tAdding electrons to SimChannel" << std::endl;
	      //std::cout << "\t\t" 
	      //	<< energyDeposit.TrackID() << " " << tdc
	      //	<< " " << xyz[0] << " " << xyz[1] << " " << xyz[2]
	      //	<< " " << fnEnDiff[k] << " " << fnElDiff[k]
	      //	<< std::endl;

	      //if(!channelPtr) std::cout << "\tUmm...ptr is NULL?" << std::endl;
	      //else std::cout << "\tChannel is " << channelPtr->Channel() << std::endl;
	      // Add the electron clusters and energy to the
	      // sim::SimChannel
	      channelPtr->AddIonizationElectrons(energyDeposit.TrackID(),
						 tdc,
						 fnElDiff[k],
						 xyz,
						 fnEnDiff[k]);

	      //std::cout << "\tAdded the electrons." << std::endl;

 	    }
	    catch(cet::exception &e) {
	      mf::LogWarning("SimDriftElectrons") << "unable to drift electrons from point ("
						  << xyz[0] << "," << xyz[1] << "," << xyz[2]
						  << ") with exception " << e;
	    } // end try to determine channel
	  } // end loop over clusters
	} // end loop over planes
      } // for each sim::SimEnergyDeposit
    /*
    // Now that we've processed the information for all the
    // sim::SimEnergyDeposit objects into sim::SimChannel objects,
    // create the associations between them. 

    // Define the container for the associations between the channels
    // and the energy deposits (steps). Note it's possible for an
    // energy deposit to be associated with more than one channel (if
    // its electrons drift to multiple wires), and a channel will
    // almost certainly have multiple energy deposits.
    std::unique_ptr< art::Assns<sim::SimEnergyDeposit, sim::SimChannel> > 
      step2channel (new art::Assns<sim::SimEnergyDeposit, sim::SimChannel>);

    // For every element in the 3-D fChannelMaps array...
    for ( auto i = fChannelMaps.begin(); i != fChannelMaps.end(); ++i ) {
      for ( auto j = i->begin(); j != i->end(); ++j ) {
	for ( auto k = j->begin(); k != j->end(); ++k ) {
	  const ChannelBookKeeping_t& bookKeeping = (*k).second;
	  const size_t channelIndex = bookKeeping.channelIndex;

	  // Create a one-to-one association between each channel and
	  // each step that created it.
	  for ( size_t m = 0; m < bookKeeping.stepList.size(); ++m)
	    {
	      const size_t edIndex = bookKeeping.stepList[m];
	      // Props to me for figuring out the following two
	      // statements. You have to look deeply in the
	      // documentation for art::Ptr and util::Associations to
	      // put this together. 
	      art::Ptr<sim::SimEnergyDeposit> energyDepositPtr( energyDepositHandle, edIndex );
	      util::CreateAssn(*this, event, *channels, energyDepositPtr, *step2channel, channelIndex);
	    }
	}
      }
    }
    */
    // Write the sim::SimChannel collection. 
    event.put(std::move(channels));
    // ... and its associations.
    //event.put(std::move(step2channel));

    return;
  }

  //----------------------------------------------------------------------------
  /*
  geo::vect::Vector_t SimDriftElectrons::RecoverOffPlaneDeposit
  (geo::vect::Vector_t const& pos, geo::PlaneGeo const& plane) const
  {
    //
    // translate the landing position on the two frame coordinates
    // ("width" and "depth")
    //
    auto const landingPos = plane.PointWidthDepthProjection(pos);
    
    //
    // compute the distance of the landing position on the two frame
    // coordinates ("width" and "depth");
    // keep the point within 10 micrometers (0.001 cm) from the border
    //
    auto const offPlane = plane.DeltaFromActivePlane(landingPos, 0.001);
    
    //
    // if both the distances are below the margin, move the point to
    // the border
    //
    
    // nothing to recover: landing is inside
    if ((offPlane.X() == 0.0) && (offPlane.Y() == 0.0)) return pos;
    
    // won't recover: too far
    if ((std::abs(offPlane.X()) > fOffPlaneMargin)
	|| (std::abs(offPlane.Y()) > fOffPlaneMargin))
      return pos;
    
    // we didn't fully decompose because it might be unnecessary;
    // now we need the full thing
    auto const distance = plane.DistanceFromPlane(pos);
    
    return plane.ComposePoint(distance, landingPos + offPlane);
    
  } // SimDriftElectrons::RecoverOffPlaneDeposit()
  */
} // namespace detsim

DEFINE_ART_MODULE(detsim::SimDriftElectrons)
