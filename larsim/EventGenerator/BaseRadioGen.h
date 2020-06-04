/// \file  BaseRadioGenHelper.cc
/// \brief Base module for radiogen/decay0gen for sampling etc
/// \author  plasorak@FNAL.GOV
///          June 2020 PLasorak
#pragma once

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

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"



namespace simb { class MCTruth; }
namespace evgen { class clhep_random; }

namespace evgen {
  class BaseRadioGen: public art::EDProducer {
  public:
    explicit BaseRadioGen(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);
    void beginRun(art::Run& run);
    void beginJob();
    void endJob();

    virtual void produce_radio(art::Event& evt) = 0;
    virtual void beginRun_radio(art::Run&) {};
    virtual void beginJob_radio() {}
    virtual void endJob_radio() {}

  protected:
    int GetNDecays();
    bool GetGoodPositionTime(TLorentzVector& position);
    TLorentzVector dirCalc(double p, double m);

    void FillHistos(simb::MCParticle& part);
    std::vector<std::string> m_isotope; ///< isotope to simulate.  Example:  "Ar39"

    CLHEP::HepRandomEngine& GetRandomEngine() { return m_engine; }
    int GetNEvents() { return m_nevent; }
    void GetTs(double& T0, double& T1) {T0=m_T0; T1=m_T1;}
    void GetXs(double& X0, double& X1) {X0=m_X0; X1=m_X1;}
    void GetYs(double& Y0, double& Y1) {Y0=m_Y0; Y1=m_Y1;}
    void GetZs(double& Z0, double& Z1) {Z0=m_Z0; Z1=m_Z1;}

  private:
    void SimplePDG(int pdg, int& simple, std::string& name);
    void DeclareOutputHistos();
    bool findNode(const TGeoNode* curnode, std::string& tgtnname,
                  const TGeoNode* & targetnode);
    bool findMotherNode(const TGeoNode* cur_node, std::string& daughter_name,
                        const TGeoNode* & mother_node);

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
    CLHEP::HepRandomEngine& m_engine;

    TGeoManager* m_geo_manager;
    double m_volume_cc;

    std::unique_ptr<CLHEP::RandFlat   > m_random_flat;
    std::unique_ptr<CLHEP::RandPoisson> m_random_poisson;

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
    TH1D* m_pdg_TH1D;


  };


  BaseRadioGen::BaseRadioGen(fhicl::ParameterSet const& pset):
    EDProducer(pset),
    m_engine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "BaseRadioGen", pset, "SeedBaseRadioGen")) {

    m_nevent=0;

    produces<std::vector<simb::MCTruth>>();
    produces<sumdata::RunData, art::InRun>();

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
        throw cet::exception("BaseRadioGen") << "Didn't find the node " << m_volume_rand << " exiting because I cannot generate events in this volume.";
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
          throw cet::exception("BaseRadioGen") << "Didn't find the mum of the following node: " << daughter_name;
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
            node_to_throw->LocalToMaster(local,master);
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


    m_random_flat    = std::make_unique<CLHEP::RandFlat   >(m_engine);
    m_random_poisson = std::make_unique<CLHEP::RandPoisson>(m_engine);


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
        throw cet::exception("BaseRadioGen") << "Didn't find the material " << m_material << " or the volume " << m_volume_gen << " in the specified volume " << m_volume_rand << ".\n";
      }

      double proportion = (double)nfound / ntries;
      std::cout << "There is " << proportion*100. << "% of " << m_material << " in the specified volume.\n";
      m_volume_cc *= proportion;
    }

  }

  void BaseRadioGen::produce(art::Event& evt) {
    m_nevent++;
    produce_radio(evt);
  }

  void BaseRadioGen::beginRun(art::Run& run) {
    art::ServiceHandle<geo::Geometry const> geo;
    run.put(std::make_unique<sumdata::RunData>(m_geo_service->DetectorName()));
    beginRun_radio(run);
  }

  void BaseRadioGen::beginJob() {
    if (m_isotope.empty()) {
      throw cet::exception("BaseRadioGen") << "m_isotope is empty, you need to fill it yourself in the constructor of your module\n";
    }

    DeclareOutputHistos();
    beginJob_radio();
    m_nevent=0;
  }

  void BaseRadioGen::endJob() {
    if (m_nevent) {
      m_pdg_TH1D->Scale(1./m_nevent);

      for (auto histo : m_pos_xy_TH2D) {
        m_pos_xy_TH2D[histo.first]->Scale(1./m_nevent);
        m_pos_xz_TH2D[histo.first]->Scale(1./m_nevent);
        m_dir_x_TH1D [histo.first]->Scale(1./m_nevent);
        m_dir_y_TH1D [histo.first]->Scale(1./m_nevent);
        m_dir_z_TH1D [histo.first]->Scale(1./m_nevent);
        m_mom_TH1D   [histo.first]->Scale(1./m_nevent);
        m_ke_TH1D    [histo.first]->Scale(1./m_nevent);
        m_time_TH1D  [histo.first]->Scale(1./m_nevent);
      }
    }
    endJob_radio();
  }

  bool BaseRadioGen::GetGoodPositionTime(TLorentzVector& position) {

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
  int BaseRadioGen::GetNDecays() {

    if (m_rate_mode) {
      return m_random_poisson->fire(m_rate);
    } else {
      double rate = abs(m_Bq * (m_T1-m_T0) * m_volume_cc / 1.0E9);
      return m_random_poisson->fire(rate);
    }

  }

  bool BaseRadioGen::findNode(const TGeoNode* curnode, std::string& tgtnname,
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


  bool BaseRadioGen::findMotherNode(const TGeoNode* cur_node, std::string& daughter_name,
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

  void BaseRadioGen::DeclareOutputHistos() {
    art::ServiceHandle<art::TFileService> tfs;
    auto pdgs = {1000020040,11,-11,22,2112,2212,9999};
    m_pdg_TH1D = tfs->make<TH1D>("PDG", ";PDG;n particles/event", pdgs.size(), 0, pdgs.size());

    for (auto pdg : pdgs) {

      std::string part="";
      int simple_pdg;
      SimplePDG(pdg, simple_pdg, part);
      art::TFileDirectory dir = tfs->mkdir(part.c_str());

      m_pos_xy_TH2D[simple_pdg] = dir.make<TH2D>("posXY"   , ";X [cm];Y [cm]"             , 100, m_X0, m_X1, 100, m_Y0, m_Y1);
      m_pos_xz_TH2D[simple_pdg] = dir.make<TH2D>("posXZ"   , ";X [cm];Z [cm]"             , 100, m_X0, m_X1, 100, m_Z0, m_Z1);
      m_dir_x_TH1D [simple_pdg] = dir.make<TH1D>("dirX"    , ";X momentum projection"     , 100,   -1,   1);
      m_dir_y_TH1D [simple_pdg] = dir.make<TH1D>("dirY"    , ";Y momentum projection"     , 100,   -1,   1);
      m_dir_z_TH1D [simple_pdg] = dir.make<TH1D>("dirZ"    , ";Z momentum projection"     , 100,   -1,   1);
      m_mom_TH1D   [simple_pdg] = dir.make<TH1D>("Momentum", ";Momentum [MeV];n particles/event", 5000,   0, 500);
      m_ke_TH1D    [simple_pdg] = dir.make<TH1D>("KE"      , ";KE [MeV];n particles/event"      , 5000,   0, 500);
      m_time_TH1D  [simple_pdg] = dir.make<TH1D>("Time"    , ";Time[ns];n particles/event"      , 100, m_T0, m_T1);

      m_pdg_TH1D->GetXaxis()->SetBinLabel(simple_pdg+1, part.c_str());
    }
  }

  void BaseRadioGen::SimplePDG(int pdg, int& simple, std::string& name) {
    switch (pdg) {
    case 1000020040:
      name = "alpha";
      simple = 0;
      return;
    case 22:
      name = "gamma";
      simple = 1;
      return;
    case -11:
      name = "positron";
      simple = 2;
      return;
    case 11:
      name = "electron";
      simple = 3;
      return;
    case 2112:
      name = "neutron";
      simple = 4;
      return;
    case 2212:
      name = "proton";
      simple = 5;
      return;
    default:
      name = "other";
      simple = 6;
      return;
    }
  }


  void BaseRadioGen::FillHistos(simb::MCParticle& part) {
    int pdg = part.PdgCode();
    int simple_pdg = 0;
    std::string dummy="";
    SimplePDG(pdg, simple_pdg, dummy);

    TLorentzVector position = part.Position();
    TLorentzVector mom = part.Momentum();
    double mass = mom.M();
    m_pos_xy_TH2D[simple_pdg]->Fill(position.X(), position.Y());
    m_pos_xz_TH2D[simple_pdg]->Fill(position.X(), position.Z());
    m_dir_x_TH1D [simple_pdg]->Fill(mom.Px()/mom.P());
    m_dir_y_TH1D [simple_pdg]->Fill(mom.Py()/mom.P());
    m_dir_z_TH1D [simple_pdg]->Fill(mom.Pz()/mom.P());
    double ke = (sqrt(mom.P()*mom.P()+mass*mass)-mass)*1000.;
    m_ke_TH1D    [simple_pdg]->Fill(ke);
    m_mom_TH1D   [simple_pdg]->Fill(mom.P()*1000.);
    m_time_TH1D  [simple_pdg]->Fill(position.T());
    m_pdg_TH1D->Fill(simple_pdg);
  }

  TLorentzVector BaseRadioGen::dirCalc(double p, double m) {
    // isotropic production angle for the decay product
    double costheta = (2.0*m_random_flat->fire() - 1.0);

    if (costheta < -1.0) costheta = -1.0;
    if (costheta > 1.0) costheta = 1.0;

    double const sintheta = sqrt(1.0-costheta*costheta);
    double const phi = 2.0*M_PI*m_random_flat->fire();

    return TLorentzVector{p*sintheta*std::cos(phi),
        p*sintheta*std::sin(phi),
        p*costheta,
        std::sqrt(p*p+m*m)};
  }

}//end namespace evgen
