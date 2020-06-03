/// \file  Decay0Gen_module.cc
/// \brief Generator for radiological decays
/// Module designed to produce a set list of particles for MC to model radiological decays
/// \author  plasorak@FNAL.GOV
///          April 2020 PLasorak

// C++ includes.
#include <string>
#include <regex>
#include <cmath>
#include <math.h>
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
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


#include <bxdecay0/i_random.h>
#include <bxdecay0/event.h>            // Decay event data model
#include <bxdecay0/decay0_generator.h> // Decay0 generator with OOP interface
#include <bxdecay0/particle.h>

namespace simb { class MCTruth; }
namespace evgen { class clhep_random; }
namespace evgen {
  /// Module to generate particles created by radiological decay, patterend off of SingleGen
  /// Currently it generates only in rectangular prisms oriented along the x,y,z axes

  class Decay0Gen : public art::EDProducer {
  public:
    explicit Decay0Gen(fhicl::ParameterSet const& pset);
    ~Decay0Gen() {
      for (auto& i: m_decay0_generator) {
        i->reset();
      }
    }
  private:
    // This is called for each event.
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    bool findNode(const TGeoNode* curnode, std::string& tgtnname, 
                  const TGeoNode* & targetnode);
    bool findMotherNode(const TGeoNode* cur_node, std::string& daughter_name,
                        const TGeoNode* & mother_node);

    int GetNDecays();
    bool GetGoodPositionTime(TLorentzVector& position);
    
    std::vector<std::string> m_isotope; ///< isotope to simulate.  Example:  "Ar39"
    std::string m_decay_chain; ///< decay chain to simulate.  Example:  "Rn222", this can also be the complete decay chain
    std::string m_material;    ///< regex of materials in which to generate the decays.  Example: "LAr"
    std::string m_volume_rand; ///< The volume in which to generate the decays
    std::string m_volume_gen;  ///< The volume in which to generate the decays
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
    double m_volume_cc;
    std::unique_ptr<CLHEP::RandFlat       > m_random_flat;
    std::unique_ptr<CLHEP::RandExponential> m_random_exponential;
    std::unique_ptr<CLHEP::RandPoisson    > m_random_poisson;
    std::shared_ptr<clhep_random          > m_random_decay0;

    // Declare a Decay0 generator:
    std::vector<std::unique_ptr<bxdecay0::decay0_generator>> m_decay0_generator;
    
    bool m_single_isotope_mode;
    bool m_geo_volume_mode;
    bool m_rate_mode;
    
    std::regex  m_regex_material;
    std::regex  m_regex_volume;
    int m_nevent;
    
    std::map<int,TH2D*> m_pos_xy_TH2D;
    std::map<int,TH2D*> m_pos_xz_TH2D;
    std::map<int,TH1D*> m_dir_x_TH1D;
    std::map<int,TH1D*> m_dir_y_TH1D;
    std::map<int,TH1D*> m_dir_z_TH1D;
    std::map<int,TH1D*> m_mom_TH1D;
    std::map<int,TH1D*> m_ke_TH1D;
    std::map<int,TH1D*> m_time_TH1D;
    TH1D* m_timediff_TH1D;
    TH1D* m_pdg_TH1D;
    
  };
    /// \brief Wrapper functor for a standard random number generator
  struct clhep_random : public bxdecay0::i_random{
    /// Constructor
    clhep_random(CLHEP::HepRandomEngine& gen):
      m_generator(gen),
      m_rand_flat(gen) { }
    
    /// Main operator
    virtual double operator()(){
      double v = m_rand_flat.fire(0.,1.);
      return v;
    }
    CLHEP::HepRandomEngine& m_generator; 
    CLHEP::RandFlat m_rand_flat;

  };

}

namespace evgen{

  Decay0Gen::Decay0Gen(fhicl::ParameterSet const& pset):
    EDProducer(pset)
  {
    std::string isotope="";
    m_single_isotope_mode = pset.get_if_present<std::string>("isotope", isotope);
    
    
    if (not m_single_isotope_mode) {
      fhicl::ParameterSet decay_chain = pset.get<fhicl::ParameterSet>("decay_chain");
      int index=0;
     
      while (decay_chain.get_if_present<std::string>("isotope_"+std::to_string(index++), isotope)) {
        m_isotope.push_back(isotope);
      }
      
    } else {
      m_isotope.push_back(isotope);
    }

    m_material = pset.get<std::string>("material", ".*");
    m_regex_material = (std::regex)m_material;

    m_volume_gen = pset.get<std::string>("volume_gen", ".*");
    m_regex_volume = (std::regex)m_volume_gen;

    m_rate_mode = pset.get_if_present<double>("rate", m_rate);
    if (not m_rate_mode)
      m_Bq = pset.get<double>("BqPercc");
    
    bool timed_mode = pset.get_if_present<double>("T0", m_T0);
    if (timed_mode) {
      m_T1 = pset.get<double>("T1");
    } else {
      const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
      int nsample = detprop->NumberTimeSamples();
      double rate = detprop->SamplingRate();
      m_T0 = -nsample * rate;
      m_T1 = -m_T0;
    }

    m_geo_volume_mode = pset.get_if_present<std::string>("volume_rand", m_volume_rand);

    m_geo_manager = m_geo_service->ROOTGeoManager();
    
    if (not m_geo_volume_mode) {
      m_X0 = pset.get<double>("X0");
      m_Y0 = pset.get<double>("Y0");
      m_Z0 = pset.get<double>("Z0");
      m_X1 = pset.get<double>("X1");
      m_Y1 = pset.get<double>("Y1");
      m_Z1 = pset.get<double>("Z1");
    } else {
      const TGeoNode* world = gGeoManager->GetNode(0);
      world->GetVolume()->SetAsTopVolume();
      const TGeoNode* node_to_throw = nullptr; // 
      bool found = findNode(world, m_volume_rand, node_to_throw);

      if (not found) {
        throw cet::exception("Decay0Gen") << "Didn't find the node " << m_volume_rand << " exiting because I cannot generate events in this volume.";
      }
  
      std::vector<const TGeoNode*> mother_nodes;
      const TGeoNode* current_node=node_to_throw;
      std::string daughter_name = node_to_throw->GetName();
      int nmax = 20;
      int iter=0;
      while (current_node != world and iter++<nmax) {
        const TGeoNode* mother_node = nullptr;
        daughter_name =current_node->GetName();
        bool found_mum = findMotherNode(world, daughter_name, mother_node);
        if(not found_mum) {
          throw cet::exception("Decay0Gen") << "Didn't find the mum of the following node: " << daughter_name;
        }
        mother_nodes.push_back(mother_node);
        current_node = mother_node;
      }
      
  
      TGeoVolume* vol   = node_to_throw->GetVolume();
      TGeoShape*  shape = vol->GetShape();
      TGeoBBox*   bbox  = (TGeoBBox*)shape;
  
      double dx = bbox->GetDX();
      double dy = bbox->GetDY();
      double dz = bbox->GetDZ();
  
      double halfs[3] = { dx, dy, dz };
      double posmin[3] = {  1.0e30,  1.0e30,  1.0e30 };
      double posmax[3] = { -1.0e30, -1.0e30, -1.0e30 };

      const double* origin = bbox->GetOrigin();
      for ( int ix = -1; ix <= 1; ix += 2) {
        for ( int iy = -1; iy <= 1; iy += 2) {
          for ( int iz = -1; iz <= 1; iz += 2) {
            double local[3];
            local[0] = origin[0] + (double)ix*halfs[0];
            local[1] = origin[1] + (double)iy*halfs[1];
            local[2] = origin[2] + (double)iz*halfs[2];
            double master[3];
            node_to_throw->LocalToMaster(local,master);//FIXME
            for (auto const& mum: mother_nodes) {
              local[0] = master[0];
              local[1] = master[1];
              local[2] = master[2];
              mum->LocalToMaster(local, master);
            }
            for ( int j = 0; j < 3; ++j ) {
              posmin[j] = TMath::Min(posmin[j],master[j]);
              posmax[j] = TMath::Max(posmin[j],master[j]);
            }
          }
        }
      }
      
      m_X0 = posmin[0];
      m_Y0 = posmin[1];
      m_Z0 = posmin[2];
      m_X1 = posmax[0];
      m_Y1 = posmax[1];
      m_Z1 = posmax[2];
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
    m_random_decay0      = std::make_shared<clhep_random          >(engine);

    for (auto const& isotope: m_isotope) {
      auto generator = std::make_unique<bxdecay0::decay0_generator>();
      generator->reset();
    
      // Configure the Decay0 generator:
      generator->set_decay_category(bxdecay0::decay0_generator::DECAY_CATEGORY_BACKGROUND);
      generator->set_decay_isotope(isotope.c_str());
      try{
        generator->initialize(*m_random_decay0);
      } catch (...) {
        throw cet::exception("Decay0Gen") << "The inialisation of Decay0 failed. Maybe the isotope " << isotope << " doesn't exists?\n";
      }
      m_decay0_generator.push_back(std::move(generator));
    }
    
    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

    m_volume_cc = (m_X1-m_X0) * (m_Y1-m_Y0) * (m_Z1-m_Z0);
    
    if (m_material != ".*" || m_volume_gen != ".*") {
      std::cout << "Calculating the proportion of " << m_material << " and the volume " << m_volume_gen << " in the specified volume " << m_volume_rand << ".\n";
      int nfound=0;
      int ntries=0;
      int npoint=10000; // 1% statistical uncertainty
      int nmax_tries=npoint*1000;
      double xyz[3];
      TGeoNode* node = nullptr;

      while (nfound<npoint and ntries<nmax_tries) {
        ntries++;
        node = nullptr;
        xyz[0] = m_X0 + (m_X1 - m_X0) * m_random_flat->fire(0, 1.);
        xyz[1] = m_Y0 + (m_Y1 - m_Y0) * m_random_flat->fire(0, 1.);
        xyz[2] = m_Z0 + (m_Z1 - m_Z0) * m_random_flat->fire(0, 1.);
        m_geo_manager->SetCurrentPoint(xyz);
        node = m_geo_manager->FindNode();
        if (!node) continue;
        if (node->IsOverlapping()) continue;

        std::string material_name = node->GetMedium()->GetMaterial()->GetName();
        std::string volume_name = node->GetVolume()->GetName();
        bool flag = std::regex_match(material_name, m_regex_material) && std::regex_match(volume_name, m_regex_volume);
        if (!flag) continue;
        nfound++;
      }
      
      if (nfound==0) {
        throw cet::exception("Decay0Gen") << "Didn't find the material " << m_material << " or the volume " << m_volume_gen << " in the specified volume " << m_volume_rand << ".\n";
      }
      
      double proportion = (double)nfound / ntries;
      std::cout << "There is " << proportion*100. << "% of " << m_material << " in the specified volume.\n";
      m_volume_cc *= proportion;
    }
    
    art::ServiceHandle<art::TFileService> tfs;
    m_timediff_TH1D = tfs->make<TH1D>("TimeDiff", ";Time Diff[ns];n particles" , (int)(m_T1/100), 0, m_T1);
    m_pdg_TH1D      = tfs->make<TH1D>("PDG"     , ";PDG;n particles"           ,  10,    0,  10);

    for (auto pdg : {1,2,3,4,5,6}) {
      std::string part="";
      switch (pdg) {
      case 1:
        part="alpha";
        break;
      case 2:
        part="gamma";
        break;
      case 3:
        part="positron";
        break;
      case 4:
        part="electron";
        break;
      case 5:
        part="neutron";
        break;
      case 6:
        part="proton";
        break;
      default:
        break;
      }
      
      m_pos_xy_TH2D[pdg]   = tfs->make<TH2D>(Form("posXY_%s"   , part.c_str()), ";X [cm];Y [cm]"             , 100, m_X0, m_X1, 100, m_Y0, m_Y1);
      m_pos_xz_TH2D[pdg]   = tfs->make<TH2D>(Form("posXZ_%s"   , part.c_str()), ";X [cm];Z [cm]"             , 100, m_X0, m_X1, 100, m_Z0, m_Z1);
      m_dir_x_TH1D [pdg]   = tfs->make<TH1D>(Form("dirX_%s"    , part.c_str()), ";X momentum projection"     , 100,   -1,   1);
      m_dir_y_TH1D [pdg]   = tfs->make<TH1D>(Form("dirY_%s"    , part.c_str()), ";Y momentum projection"     , 100,   -1,   1);
      m_dir_z_TH1D [pdg]   = tfs->make<TH1D>(Form("dirZ_%s"    , part.c_str()), ";Z momentum projection"     , 100,   -1,   1);
      m_mom_TH1D   [pdg]   = tfs->make<TH1D>(Form("Momentum_%s", part.c_str()), ";Momentum [MeV];n particles", 5000,   0, 500);
      m_ke_TH1D    [pdg]   = tfs->make<TH1D>(Form("KE_%s"      , part.c_str()), ";KE [MeV];n particles"      , 5000,   0, 500);
      m_time_TH1D  [pdg]   = tfs->make<TH1D>(Form("Time_%s"    , part.c_str()), ";Time[ns];n particles"      , 100, m_T0, m_T1);
          
      m_pdg_TH1D->GetXaxis()->SetBinLabel(pdg+1, part.c_str());
    }
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

    //unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);
    int track_id=-1;
    const std::string primary_str("primary");
    
    for (auto const& decay0_gen: m_decay0_generator) {
      int n_decay = GetNDecays();
      
      for (int iDecay=0; iDecay<n_decay; ++iDecay) {
        TLorentzVector position;
        if (GetGoodPositionTime(position)) {
        
          bxdecay0::event gendecay;     // Declare an empty decay event
          decay0_gen->shoot(*m_random_decay0, gendecay); // Randomize the decay event
          std::vector<bxdecay0::particle> part = gendecay.get_particles();
    
          double t_alpha=0;
          double t_electron=0;
          for (auto const& p: part) {
            // alpha particles need a little help since they're not in the TDatabasePDG table
            // so don't rely so heavily on default arguments to the MCParticle constructor
            simb::MCParticle part;
            if (not p.is_valid()) {
              p.print(std::cout);
              throw cet::exception("Decay0Gen") << "Invalid part generated by Decay0, printed above (no clue what that means so throw)";
            }
          
            double mass = bxdecay0::particle_mass_MeV(p.get_code());
            int simple_pdg=0;
          
            if      (p.is_alpha   ()) { simple_pdg = 1; part = simb::MCParticle(track_id, 1000020040, primary_str,-1,mass,1); }
            else if (p.is_gamma   ()) { simple_pdg = 2; part = simb::MCParticle(track_id,         22, primary_str); }
            else if (p.is_positron()) { simple_pdg = 3; part = simb::MCParticle(track_id,        -11, primary_str); }
            else if (p.is_electron()) { simple_pdg = 4; part = simb::MCParticle(track_id,         11, primary_str); }
            else if (p.is_neutron ()) { simple_pdg = 5; part = simb::MCParticle(track_id,       2112, primary_str); }
            else if (p.is_proton  ()) { simple_pdg = 6; part = simb::MCParticle(track_id,       2212, primary_str); }
            else {
              p.print(std::cout);
              throw cet::exception("Decay0Gen") << "Particle above is weird, cannot recognise it.";
            }
            if ((p.is_positron() or p.is_electron()) and t_electron == 0) {
              t_electron = p.get_time();
            }
            if (p.is_alpha() and t_alpha == 0) {
              t_alpha = p.get_time();
            }
            track_id--;
            
            TLorentzVector mom(p.get_px()/1000.,
                               p.get_py()/1000.,
                               p.get_pz()/1000.,
                               sqrt(p.get_p()*p.get_p() + mass*mass)/1000.);
            TLorentzVector this_part_position = position;
            double t = position.T();
            this_part_position.SetT(t+p.get_time()*1e9);

            part.AddTrajectoryPoint(this_part_position, mom);
            truth.Add(part);
          
            m_pos_xy_TH2D[simple_pdg]->Fill(position.X(), position.Y());
            m_pos_xz_TH2D[simple_pdg]->Fill(position.X(), position.Z());
            m_dir_x_TH1D [simple_pdg]->Fill(mom.Px()/mom.P());
            m_dir_y_TH1D [simple_pdg]->Fill(mom.Py()/mom.P());
            m_dir_z_TH1D [simple_pdg]->Fill(mom.Pz()/mom.P());
            double ke = (sqrt(mom.P()*mom.P()*1000.*1000.+mass*mass)-mass);
            m_ke_TH1D    [simple_pdg]->Fill(ke);
            m_mom_TH1D   [simple_pdg]->Fill(mom.P()*1000.);
            m_time_TH1D  [simple_pdg]->Fill(position.T());
            m_pdg_TH1D->Fill(simple_pdg);

          }
          m_timediff_TH1D->Fill(abs(t_alpha-t_electron)/1e9);
          gendecay.reset();
        } // GetGoodPosition
      } // idecay
    } // m_decay_generators
    
    MF_LOG_DEBUG("Decay0Gen") << truth;
    truthcol->push_back(truth);
    evt.put(std::move(truthcol));
  }

  bool Decay0Gen::GetGoodPositionTime(TLorentzVector& position) {

    double time = m_T0 + m_random_flat->fire()*(m_T1 - m_T0);
    int ntries=0;
    int nmaxtries=10000;
    bool flag=0;
    while (ntries++<nmaxtries) {
      position.SetXYZT(m_X0 + m_random_flat->fire()*(m_X1 - m_X0),
                       m_Y0 + m_random_flat->fire()*(m_Y1 - m_Y0),
                       m_Z0 + m_random_flat->fire()*(m_Z1 - m_Z0),
                       time);
      
      const TGeoNode* node = m_geo_manager->FindNode(position.X(),position.Y(),position.Z());
      std::string material_name = node->GetMedium()->GetMaterial()->GetName();
      std::string volume_name = node->GetVolume()->GetName();
      flag = std::regex_match(material_name, m_regex_material) && std::regex_match(volume_name, m_regex_volume);
      if (flag) return true;
    }
    return false;
  }


  //____________________________________________________________________________
  // Generate radioactive decays per isotope per volume according to the FCL parameters
  int Decay0Gen::GetNDecays() {

    if (m_rate_mode) {
      return m_random_poisson->fire(m_rate);
    } else {
      double rate = abs(m_Bq * (m_T1-m_T0) * m_volume_cc / 1.0E9);
      return m_random_poisson->fire(rate);
    }

  }
  
  bool Decay0Gen::findNode(const TGeoNode* curnode, std::string& tgtnname,
                           const TGeoNode* & targetnode)
  /// Shamelessly stolen from here: https://cdcvs.fnal.gov/redmine/attachments/6719/calc_bbox.C
  {
    std::string nname = curnode->GetName();
    std::string vname = curnode->GetVolume()->GetName();
    if ( nname == tgtnname || vname == tgtnname ) {
      targetnode = curnode;
      return true;
    }

    TObjArray* daunodes = curnode->GetNodes();
    if ( ! daunodes ) return false;
    TIter next(daunodes);
    const TGeoNode* anode = 0;
    while ( (anode = (const TGeoNode*)next()) ) {
      bool found = findNode(anode,tgtnname,targetnode);
      if ( found ) return true;
    }

    return false;
  }

  
  bool Decay0Gen::findMotherNode(const TGeoNode* cur_node, std::string& daughter_name,
                                 const TGeoNode* & mother_node)
  /// Adapted from above
  {
    TObjArray* daunodes = cur_node->GetNodes();
  
    if ( ! daunodes ) return false;
  
    TIter next(daunodes);

    const TGeoNode* anode = 0;
    bool found = 0;
    while ( (anode = (const TGeoNode*)next()) ) {
      std::string nname = anode->GetName();
      std::string vname = anode->GetVolume()->GetName();

      if ( nname == daughter_name || vname == daughter_name ) {
        mother_node = cur_node;
        return true;
      }
      found = findMotherNode(anode, daughter_name, mother_node);
    
    }

    return found;
  }

  
}//end namespace evgen

DEFINE_ART_MODULE(evgen::Decay0Gen)
