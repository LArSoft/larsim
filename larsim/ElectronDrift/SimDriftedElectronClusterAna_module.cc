////////////////////////////////////////////////////////////////////////
// Class:        SimDriftedElectronClusterAna
// Plugin Type: analyzer (art v2_05_00)
// File:         DriftedElectronsClustersSAna_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art_root_io/TFileService.h"
#include "TNtuple.h"


#include "lardataobj/Simulation/SimDriftedElectronCluster.h"



namespace detsim {

class SimDriftedElectronClusterAna : public art::EDAnalyzer {
public:
  explicit  SimDriftedElectronClusterAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
   SimDriftedElectronClusterAna( SimDriftedElectronClusterAna const &) = delete;
   SimDriftedElectronClusterAna( SimDriftedElectronClusterAna &&) = delete;
   SimDriftedElectronClusterAna & operator = ( SimDriftedElectronClusterAna const &) = delete;
   SimDriftedElectronClusterAna & operator = ( SimDriftedElectronClusterAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:


  art::InputTag         fEDepTag;

  TNtuple* fNtuple;
};


SimDriftedElectronClusterAna:: SimDriftedElectronClusterAna(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
  , fEDepTag{p.get<art::InputTag>("EDepModuleLabel")}
{}

void SimDriftedElectronClusterAna::analyze(art::Event const & e)
{
  auto const& edep_handle = e.getValidHandle< std::vector<sim::SimDriftedElectronCluster> >(fEDepTag);
  auto const& edep_vec(*edep_handle);
  std::cout<< "=====================edep"<<edep_vec.size()<<std::endl;
  for(auto const& edep : edep_vec){
    fNtuple->Fill(e.run(),e.event(),
		  edep.NumberOfElectrons(),
		  edep.Time());
  }
}

void SimDriftedElectronClusterAna::beginJob()
{
  art::ServiceHandle<art::TFileService const> tfs;
  fNtuple = tfs->make<TNtuple>("nt_is","EDep IS Calc Ntuple","run:event:ne:t");
}
} // namespace detsim

DEFINE_ART_MODULE(detsim::SimDriftedElectronClusterAna)
