////////////////////////////////////////////////////////////////////////
// Class:       PDFastSimPVS
// Plugin Type: producer
// File:        PDFastSimPVS_module.cc
// Description:
// - acts on sim::SimEnergyDeposit from LArG4Main, 
// - simulate (fast, photon visibility service) the OpDet response to optical photons
// Input: 'sim::SimEnergyDeposit'
// Output: 'sim::OpDetBacktrackerRecord'
//Fast simulation of propagating the photons created from SimEnergyDeposits.

//This module does a fast simulation of propagating the photons created from SimEnergyDeposits,
//This simulation is done using the PhotonLibrary, which stores the visibilities of each optical channel
//with respect to each optical voxel in the TPC volume, to avoid propagating single photons using Geant4.
//At the end of this module a collection of the propagated photons either as
//'sim::OpDetBacktrackerRecord' are placed into the art event.

//The steps this module takes are:
//  - to take number of photon and the vertex information from 'sim::SimEnergyDeposits',
//  - use the PhotonLibrary (visibilities) to determine the amount of visible photons at each optical channel,
//  - visible photons: the number of photons times the visibility at the middle of the Geant4 step for a given optical channel.
//  - other photon information is got from 'sim::SimEnergyDeposits'
//  - add 'sim::OpDetBacktrackerRecord' to event
// Aug. 19 by Mu Wei
////////////////////////////////////////////////////////////////////////

// Art libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"

// Random number engine
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
//#include "CLHEP/Random/RandGauss.h"

namespace phot
{
    class PDFastSimPVS : public art::EDProducer
    {
    public:
        explicit PDFastSimPVS(fhicl::ParameterSet const&);
        void produce(art::Event&) override;
        void AddOpDetBTR(std::vector< sim::OpDetBacktrackerRecord > & opbtr,
                             std::map<int, int> & ChannelMap,
                             sim::OpDetBacktrackerRecord btr);
    private:
        bool                          fDoSlowComponent;
        art::InputTag                 simTag;
        std::unique_ptr<ScintTime>    fScintTime;        // Tool to retrive timinig of scintillation        
        CLHEP::HepRandomEngine&       fPhotonEngine;
        CLHEP::HepRandomEngine&       fScintTimeEngine;
        std::map<int, int>            PDChannelToSOCMap; //Where each OpChan is.
    };
    
    //......................................................................    
    PDFastSimPVS::PDFastSimPVS(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fDoSlowComponent{pset.get<bool>("DoSlowComponent")}
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fScintTime{art::make_tool<ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", pset, "SeedScintTime"))
    {
        std::cout << "PDFastSimPVS Module Construct" << std::endl;
        
        produces< std::vector<sim::SimPhotonsLite> >("pvs");
        produces< std::vector<sim::OpDetBacktrackerRecord> >("pvs");     
    }
    
    //......................................................................    
    void PDFastSimPVS::produce(art::Event& event)
    {
        std::cout << "PDFastSimPVS Module Producer" << std::endl;
        
        art::ServiceHandle<PhotonVisibilityService const> pvs;
        //unused auto const* larp = lar::providerFrom<detinfo::LArPropertiesService>();
        auto const nOpChannels = pvs->NOpChannels();
        
        CLHEP::RandPoissonQ randpoisphot{fPhotonEngine};
//        CLHEP::RandFlat randflatscinttime{fScintTimeEngine};
        
        std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr  (new std::vector<sim::OpDetBacktrackerRecord>);
        std::unique_ptr< std::vector< sim::SimPhotonsLite> >          phlit  (new std::vector<sim::SimPhotonsLite>);
        
        auto& photonLiteCollection (*phlit);
        photonLiteCollection.resize(nOpChannels);
        for (unsigned int i = 0; i < nOpChannels; i ++)
        {
            photonLiteCollection[i].OpChannel = i;
        }
        
        art::Handle< std::vector<sim::SimEnergyDeposit> > edepHandle;
        if (!event.getByLabel(simTag, edepHandle))
        {
            std::cout << "PDFastSimPVS Module Cannot getByLabel: " << simTag << std::endl;
            return;
        }
                
        art::ServiceHandle<geo::Geometry> geom;
        auto const& edeps = edepHandle;
        
        int num_points    = 0;
        int num_fastph    = 0;
        int num_slowph    = 0;
        int num_tot       = 0;
        int num_fastdp    = 0;
        int num_slowdp    = 0;
        
        float vis_scale   = 1.0; // to scale the visibility fraction, for test only;
        
        for (auto const& edepi: *edeps)
        {
            num_points ++;
            
            auto const& prt          = edepi.MidPoint();
            auto const& Visibilities = pvs->GetAllVisibilities(prt);
            if (!Visibilities)
            {
                //throw cet::exception("PDFastSimPVS")
                std::cout << "There is no entry in the PhotonLibrary for this position in space. "
                "Position: " << edepi.MidPoint();
                std::cout << "\n Move to next point" << std::endl;
                continue;
            }
            
            int trackID       = edepi.TrackID();
            int nphot         = edepi.NumPhotons();
            double edeposit   = edepi.Energy()/nphot;
            double pos[3]     = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
            
            int nphot_fast    = edepi.NumFPhotons();
            int nphot_slow    = edepi.NumSPhotons();

            num_tot    += nphot;
            num_fastph += nphot_fast;
            num_slowph += nphot_slow;
            
            for (unsigned int channel = 0; channel < nOpChannels; ++ channel)
            {
                auto visibleFraction = Visibilities[channel] * vis_scale;  // to scale the visibility fraction;
                
                if (visibleFraction == 0.0)
                {
                    continue; //voxel is not visible at this optical channel.
                }
                
                sim::OpDetBacktrackerRecord tmpbtr(channel);
                
                if (nphot_fast > 0)
                {
                    //random number, poisson distribution, mean: the amount of photons visible at this channel
                    auto n = static_cast<int>(randpoisphot.fire(nphot_fast * visibleFraction));
                    num_fastdp += n;                    
                    for (long i = 0; i < n; ++i) 
                    {
                        //calculates the time at which the photon was produced
                        fScintTime->GenScintTime(true, fScintTimeEngine);
                        auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                        ++ photonLiteCollection[channel].DetectedPhotons[time];
                        tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                        
                    }
                }
                
                if ((nphot_slow > 0) && fDoSlowComponent) 
                {
                    auto n = static_cast<int>(randpoisphot.fire(nphot_slow * visibleFraction));
                    num_slowdp += n;                    
                    for (long i = 0; i < n; ++i) 
                    {
                        fScintTime->GenScintTime(false, fScintTimeEngine);
                        auto time = static_cast<int>(edepi.StartT()+ fScintTime->GetScintTime());
                        ++ photonLiteCollection[channel].DetectedPhotons[time];
                        tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                    }
                }
                
                AddOpDetBTR(*opbtr, PDChannelToSOCMap, tmpbtr);
            }
        }
        
        std::cout << "Total points: " << num_points << ", total fast photons: " << num_fastph << ", total slow photons: " << num_slowph << ", total photon: " << num_tot << std::endl;
        std::cout << "Detected fast photons: " << num_fastdp << ", detected slow photons: " << num_slowdp << std::endl;
        
        PDChannelToSOCMap.clear();
        event.put(move(phlit), "pvs");
        event.put(move(opbtr), "pvs");
        
        return;
    }
    
    //......................................................................    
    void PDFastSimPVS::AddOpDetBTR(std::vector< sim::OpDetBacktrackerRecord > & opbtr,
                                         std::map<int, int> & ChannelMap,
                                         sim::OpDetBacktrackerRecord btr) 
    {
        int iChan = btr.OpDetNum();
        std::map<int, int>::iterator channelPosition = ChannelMap.find(iChan);
        
        if (channelPosition == ChannelMap.end() )
        {
            ChannelMap[iChan] = opbtr.size();
            opbtr.emplace_back(std::move(btr));
        }
        else
        {
            unsigned int idtest = channelPosition->second;
            auto const& timePDclockSDPsMap = btr.timePDclockSDPsMap();
            
            for(auto const& timePDclockSDP : timePDclockSDPsMap)
            {
                for(auto const& sdp : timePDclockSDP.second)
                {
                    double xyz[3] = {sdp.x, sdp.y, sdp.z};
                    opbtr.at(idtest).AddScintillationPhotons(sdp.trackID,
                                                             timePDclockSDP.first,
                                                             sdp.numPhotons,
                                                             xyz,
                                                             sdp.energy);
                }
            }
        }    
    }  

/*
template <typename Point> MappedCounts_t GetAllVisibilities(Point const& p, bool wantReflected=false ) const
{
    return doGetAllVisibilities(geo::vect::toPoint(p), wantReflected); 
}
auto PhotonVisibilityService::doGetAllVisibilities(geo::Point_t const& p, bool wantReflected) const -> MappedCounts_t
{
    phot::IPhotonLibrary::Counts_t data{};
    
    // first we fill a container of visibilities in the library index space
    // (it is directly the values of the library unless interpolation is
    //  requested)
    if(fInterpolate)
    {
        // this is a punch into multithreading face:
        static std::vector<float> ret;
        ret.resize(fMapping->libraryMappingSize(p));
        for(std::size_t libIndex = 0; libIndex < ret.size(); ++libIndex) 
        {
            ret[libIndex] = doGetVisibilityOfOpLib(p, LibraryIndex_t(libIndex), wantReflected);
        }
        data = &ret.front();
    }
    else
    {
        auto const VoxID = VoxelAt(p);
        data = GetLibraryEntries(VoxID, wantReflected);
    }
    return fMapping->applyOpDetMapping(p, data);
}

float PhotonVisibilityService::doGetVisibilityOfOpLib(geo::Point_t const& p, LibraryIndex_t libIndex, bool wantReflected ) const
{
    if(!fInterpolate) 
    {
        return GetLibraryEntry(VoxelAt(p), libIndex, wantReflected);
    }
    
    // In case we're outside the bounding box we'll get a empty optional list.
    auto const neis = fVoxelDef.GetNeighboringVoxelIDs(LibLocation(p));
    if (!neis) return 0.0;
    
    // Sum up all the weighted neighbours to get interpolation behaviour
    float vis = 0.0;
    for(const sim::PhotonVoxelDef::NeiInfo& n: neis.value()) 
    {
        if (n.id < 0) continue;
        vis += n.weight * GetLibraryEntry(n.id, libIndex, wantReflected);
    }
    
    return vis;
}

float PhotonVisibilityService::GetLibraryEntry(int VoxID, OpDetID_t libOpChannel, bool wantReflected) const
{
    if(fTheLibrary == 0)
    LoadLibrary();
    
    if(!wantReflected)
    return fTheLibrary->GetCount(VoxID, libOpChannel);
    else
    return fTheLibrary->GetReflCount(VoxID, libOpChannel);
}

void PhotonVisibilityService::LoadLibrary() const
{
// Don't do anything if the library has already been loaded.

if(fTheLibrary == 0) {

if((!fLibraryBuildJob)&&(!fDoNotLoadLibrary)) {
std::string LibraryFileWithPath;
cet::search_path sp("FW_SEARCH_PATH");

if( !sp.find_file(fLibraryFile, LibraryFileWithPath) )
throw cet::exception("PhotonVisibilityService") << "Unable to find photon library in "  << sp.to_string() << "\n";

if(!fParameterization) {
art::ServiceHandle<geo::Geometry const> geom;

mf::LogInfo("PhotonVisibilityService") << "PhotonVisibilityService Loading photon library from file "
                       << LibraryFileWithPath
                       << " for "
                       << GetVoxelDef().GetNVoxels()
                       << " voxels and "
                       << geom->NOpDets()
                       << " optical detectors."
                       << std::endl;

if(fHybrid){
fTheLibrary = new PhotonLibraryHybrid(LibraryFileWithPath,
                        GetVoxelDef());
}
else{
PhotonLibrary* lib = new PhotonLibrary;
fTheLibrary = lib;

size_t NVoxels = GetVoxelDef().GetNVoxels();
lib->LoadLibraryFromFile(LibraryFileWithPath, NVoxels, fStoreReflected, fStoreReflT0, fParPropTime_npar, fParPropTime_MaxRange);
}
}
}
else {
art::ServiceHandle<geo::Geometry const> geom;

size_t NOpDets = geom->NOpDets();
size_t NVoxels = GetVoxelDef().GetNVoxels();
mf::LogInfo("PhotonVisibilityService") << " Vis service running library build job.  Please ensure "
                     << " job contains LightSource, LArG4, SimPhotonCounter"<<std::endl;
PhotonLibrary* lib = new PhotonLibrary;
fTheLibrary = lib;

lib->CreateEmptyLibrary(NVoxels, NOpDets, fStoreReflected, fStoreReflT0, fParPropTime_npar);
}

}
}
*/
} // namespace

DEFINE_ART_MODULE(phot::PDFastSimPVS)
