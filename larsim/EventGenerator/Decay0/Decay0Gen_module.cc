/// \file  Decay0Gen_module.cc
/// \brief Generator for radiological decays
/// Module designed to produce a set list of particles for MC to model radiological decays
/// \author  plasorak@FNAL.GOV
///          April 2020 PLasorak

// C++ includes.
#include <string>
#include <regex>
#include <cmath>
#include <memory>
#include <iterator>
#include <sys/stat.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoNode.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "cetlib/exempt_ptr.h"

// nurandom includes
#include "nurandom/RandomUtils/NuRandomService.h"

// nusimdata, nugen includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nugen/EventGeneratorBase/evgenbase.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"

// root includes

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TVector3.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"


#include <bxdecay0/std_random.h>       // Wrapper for the standard random PRNG
#include <bxdecay0/event.h>            // Decay event data model
#include <bxdecay0/decay0_generator.h> // Decay0 generator with OOP interface
#include <bxdecay0/particle.h>

namespace simb { class MCTruth; }

namespace evgen {
  /// Module to generate particles created by radiological decay, patterend off of SingleGen
  /// Currently it generates only in rectangular prisms oriented along the x,y,z axes

  class Decay0Gen : public art::EDProducer {
  public:
    explicit Decay0Gen(fhicl::ParameterSet const& pset);

  private:
    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);


    
    int GetNDecays();
    bool GetGoodPosition(TVector3& position);
    bool GetGoodPositionInGeoVolume(TVector3& position);
    
    std::string m_isotope;     ///< isotope to simulate.  Example:  "Ar39"
    std::string m_decay_chain; ///< decay chain to simulate.  Example:  "Rn222", this can also be the complete decay chain
    std::string m_material;    ///< regex of materials in which to generate the decays.  Example: "LAr"
    std::string m_volume;      ///< The volume in which to generate the decays
    double      m_Bq;          ///< Radioactivity in Becquerels (decay per sec) per cubic cm.
    double      m_rate;        ///< Radioactivity in Becquerels (decay per sec) use either of this of Bq
    double      m_T0;          ///< Beginning of time window to simulate in ns
    double      m_T1;          ///< End of time window to simulate in ns
    double      m_X0;          ///< Bottom corner x position (cm) in world coordinates
    double      m_Y0;          ///< Bottom corner y position (cm) in world coordinates
    double      m_Z0;          ///< Bottom corner z position (cm) in world coordinates
    double      m_X1;          ///< Top corner x position (cm) in world coordinates
    double      m_Y1;          ///< Top corner y position (cm) in world coordinates
    double      m_Z1;          ///< Top corner z position (cm) in world coordinates

    art::ServiceHandle<geo::Geometry const> m_geo_service;
    TGeoManager* m_geo_manager;
    TGeoVolume* m_geo_volume;
    std::unique_ptr<CLHEP::RandFlat       > m_random_flat;
    std::unique_ptr<CLHEP::RandExponential> m_random_exponential;
    std::unique_ptr<CLHEP::RandPoisson    > m_random_poisson;

    // Declare a Decay0 generator:
    std::unique_ptr<bxdecay0::decay0_generator> m_decay0_generator;
    std::unique_ptr<bxdecay0::std_random> m_random_decay0;

    bool m_single_isotope_mode;
    bool m_geo_volume_mode;
    bool m_rate_mode;
    
    std::regex  m_regex_material;
    int m_nevent;
    TH2D* m_pos_xy_TH2D;
    TH2D* m_pos_xz_TH2D;
    TH1D* m_dir_x_TH1D ;
    TH1D* m_dir_y_TH1D ;
    TH1D* m_dir_z_TH1D ;
    TH1D* m_pdg_TH1D   ;
    TH1D* m_mom_TH1D   ;
    TH1D* m_time_TH1D  ;
    TH1D* m_timediff_TH1D;
  };
}

namespace evgen{

  Decay0Gen::Decay0Gen(fhicl::ParameterSet const& pset):
    EDProducer(pset) {

    m_single_isotope_mode = pset.get_if_present<std::string>("isotope", m_isotope);
    
    if (not m_single_isotope_mode) {
      m_decay_chain = pset.get<std::string>("decay_chain");
    }

    m_material = pset.get<std::string>("material", "*");
    m_regex_material = (std::regex)m_material;

    m_rate_mode = pset.get_if_present<double>("rate", m_rate);
    if (not m_rate_mode)
      m_Bq = pset.get<double>("BqPercc");
    
    bool timed_mode = pset.get_if_present<double>("T0", m_T0);
    if (timed_mode) {
      m_T1 = pset.get<double>("T1");
    } else {
      throw cet::exception("Decay0Gen") << "for now, you have to specify T0, later, it should generate for all time in the events using detector property or so.";
    }

    m_geo_volume_mode = pset.get_if_present<std::string>("volume", m_volume);

    m_geo_manager = m_geo_service->ROOTGeoManager();

    if (not m_geo_volume_mode) {
      m_X0 = pset.get<double>("X0");
      m_Y0 = pset.get<double>("Y0");
      m_Z0 = pset.get<double>("Z0");
      m_X1 = pset.get<double>("X1");
      m_Y1 = pset.get<double>("Y1");
      m_Z1 = pset.get<double>("Z1");
    } else {
      m_geo_volume = nullptr;
      m_geo_volume = m_geo_manager->FindVolumeFast(m_volume.c_str());
      if (not m_geo_volume) {
        throw cet::exception("Decay0Gen") << "Cannot find volume " << m_volume << " in the geometry";
      }
    }
    auto& Seeds = *(art::ServiceHandle<rndm::NuRandomService>());

    // declare an engine; NuRandomService associates an (unknown) engine, in
    // the current module and an instance name, with a seed (returned)
    auto const seed = Seeds.declareEngine();

    // now create the engine (for example, use art); seed will be set
    auto& engine = createEngine(seed, "HepJamesRandom", "Decay0GenModule");

    // finally, complete the registration; seed will be set again
    Seeds.defineEngine(engine);
    m_random_flat        = std::make_unique<CLHEP::RandFlat       >(engine);
    m_random_exponential = std::make_unique<CLHEP::RandExponential>(engine);
    m_random_poisson     = std::make_unique<CLHEP::RandPoisson    >(engine);


    int seed_std = m_random_flat->fire(UINT_MAX);
    std::default_random_engine generator(seed_std); // Standard PRNG
    m_random_decay0 = std::make_unique<bxdecay0::std_random>(generator);       // PRNG wrapper
    m_decay0_generator = std::make_unique<bxdecay0::decay0_generator>();
    
    // Configure the Decay0 generator:
    m_decay0_generator->set_decay_category(bxdecay0::decay0_generator::DECAY_CATEGORY_BACKGROUND);
    m_decay0_generator->set_decay_isotope(m_isotope.c_str());
    m_decay0_generator->initialize(*m_random_decay0);

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
  }

  //____________________________________________________________________________
  void Decay0Gen::beginRun(art::Run& run)
  {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(m_geo_service->DetectorName()));
  }

  //____________________________________________________________________________
  void Decay0Gen::produce(art::Event& evt)
  {
    m_nevent++;

    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    int n_decay = GetNDecays();
    std::cout << "\033[32mn_decay : " << n_decay << "\033[0m\n";
    for (int iDecay=0; iDecay<n_decay; ++iDecay) {
      if (iDecay%100==0) std::cout << "\033[32miDecay : " << iDecay << "\033[0m\n";
      bxdecay0::event gendecay;     // Declare an empty decay event
      m_decay0_generator->shoot(*m_random_decay0, gendecay); // Randomize the decay event
      std::vector<bxdecay0::particle> part = gendecay.get_particles();
      for (auto const& p: part) {
        
      }
    }
    
    MF_LOG_DEBUG("Decay0Gen") << truth;
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
  }

  bool Decay0Gen::GetGoodPositionInGeoVolume(TVector3& position) {

    const TGeoShape *shape = m_geo_volume->GetShape();
    m_geo_volume->SetAsTopVolume();
    TGeoBBox *box = (TGeoBBox *)shape;
    double dx = box->GetDX();
    double dy = box->GetDY();
    double dz = box->GetDZ();
    double ox = (box->GetOrigin())[0];
    double oy = (box->GetOrigin())[1];
    double oz = (box->GetOrigin())[2];
    double xyz[3];
    TGeoNode *node = 0;
    int i=0;

    int npoints=10000;
    
    while (i<npoints) {
      xyz[0] = ox-dx+2*dx*m_random_flat->fire();
      xyz[1] = oy-dy+2*dy*m_random_flat->fire();
      xyz[2] = oz-dz+2*dz*m_random_flat->fire();

      m_geo_manager->SetCurrentPoint(xyz);
      node = m_geo_manager->FindNode();

      if (!node) continue;
      if (node->IsOverlapping()) continue;
      i++;

      position.SetXYZ(xyz[0],xyz[1],xyz[2]);

      std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
      bool flag = std::regex_match(volmaterial, m_regex_material);

      return flag;
    }
    throw art::Exception(art::errors::LogicError) << "Couldn't find a point in the volume: "
                                                  << m_volume << " which has material: "
                                                  << m_material << " after 10000 attempts.\n";
    return false;
  }

  bool Decay0Gen::GetGoodPosition(TVector3& position) {
    
    if (m_geo_volume_mode) {
      return GetGoodPositionInGeoVolume(position);
    }

    position.SetXYZ(m_X0 + m_random_flat->fire()*(m_X1 - m_X0),
                    m_Y0 + m_random_flat->fire()*(m_Y1 - m_Y0),
                    m_Z0 + m_random_flat->fire()*(m_Z1 - m_Z0));

    std::string volmaterial = m_geo_manager->FindNode(position.X(),position.Y(),position.Z())->GetMedium()->GetMaterial()->GetName();
    bool flag = std::regex_match(volmaterial, m_regex_material);

    return flag;
  }


  //____________________________________________________________________________
  // Generate radioactive decays per isotope per volume according to the FCL parameters
  int Decay0Gen::GetNDecays() {

    if (m_rate_mode) {
      return m_random_poisson->fire(m_rate);
    } else if (m_geo_volume_mode) {
      const TGeoShape *shape = m_geo_volume->GetShape();
      TGeoBBox *box = (TGeoBBox *)shape;
      double dx = box->GetDX();
      double dy = box->GetDY();
      double dz = box->GetDZ();
      double rate = abs(m_Bq * (m_T1-m_T0) * dx * dy * dz / 1.0E9);
      std::cout << "\033[32mrate : " << rate << "\033[0m\n";
      int n = m_random_poisson->fire(rate);
      std::cout << "\033[32mn : " << n << "\033[0m\n";
      return n;
    } else {
      double rate = abs(m_Bq * (m_T1-m_T0) * (m_X1-m_X0) * (m_Y1-m_Y0) * (m_Z1-m_Z0) / 1.0E9);
      return m_random_poisson->fire(rate);
    }

  }

}//end namespace evgen

DEFINE_ART_MODULE(evgen::Decay0Gen)
