////////////////////////////////////////////////////////////////////////
//
//  MCShowerRecoAlg source
//
//
////////////////////////////////////////////////////////////////////////

#include "MCShowerRecoAlg.h"

#include "RtypesCore.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/MCBase/MCLimits.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "larsim/MCSTReco/MCRecoEdep.h"
#include "larsim/MCSTReco/MCRecoPart.h"
#include "larsim/MCSTReco/MCShowerRecoPart.h"

#include <cmath> // std::hypot()
#include <memory>
#include <vector>

namespace sim {

  //##################################################################
  MCShowerRecoAlg::MCShowerRecoAlg(fhicl::ParameterSet const& pset)
    : fPartAlg(pset.get<fhicl::ParameterSet>("MCShowerRecoPart"))
    , fDebugMode(pset.get<bool>("DebugMode"))
    , fMinShowerEnergy(pset.get<double>("MinShowerEnergy"))
    , fMinNumDaughters(pset.get<unsigned int>("MinNumDaughters"))
  //##################################################################
  {}

  std::unique_ptr<std::vector<sim::MCShower>> MCShowerRecoAlg::Reconstruct(MCRecoPart& part_v,
                                                                           MCRecoEdep& edep_v)
  {

    art::ServiceHandle<geo::Geometry const> geo;

    auto pindex = details::createPlaneIndexMap();

    fPartAlg.ConstructShower(part_v);
    auto result = std::make_unique<std::vector<sim::MCShower>>();
    auto& mcshower = *result;
    //std::vector<sim::MCShower> mcshower;
    // Get shower info from grouped particles
    const std::vector<unsigned int> shower_index_v = fPartAlg.ShowerMothers();
    mcshower.reserve(shower_index_v.size());
    std::vector<size_t> mcs_to_spart_v;
    mcs_to_spart_v.reserve(shower_index_v.size());

    bool daughter_stored = false;
    for (size_t shower_index = 0; shower_index < shower_index_v.size(); ++shower_index) {

      unsigned int shower_candidate = shower_index_v.at(shower_index);
      auto const& shower_part = part_v.at(shower_candidate);

      unsigned int mother_track = part_v.MotherTrackID(shower_candidate);
      unsigned int ancestor_track = part_v.AncestorTrackID(shower_candidate);

      if (mother_track == kINVALID_UINT || ancestor_track == kINVALID_UINT)

        throw cet::exception(__FUNCTION__) << "LOGIC ERROR: mother/ancestor track ID is invalid!";

      MCMiniPart mother_part;
      MCMiniPart ancestor_part;

      unsigned int mother_index = part_v.TrackToParticleIndex(mother_track);
      unsigned int ancestor_index = part_v.TrackToParticleIndex(ancestor_track);

      if (mother_index != kINVALID_UINT)
        mother_part = part_v[mother_index];
      else
        mother_part._track_id = mother_track;

      if (ancestor_index != kINVALID_UINT)
        ancestor_part = part_v[ancestor_index];
      else
        ancestor_part._track_id = ancestor_track;

      double shower_g4_energy = shower_part._start_mom[3];

      if (fDebugMode)

        std::cout << "Found MCShower with mother energy: " << shower_g4_energy << " MeV";

      // Skip if mother energy is less than the enery threshold
      if (shower_g4_energy < fMinShowerEnergy) {
        if (fDebugMode) std::cout << " ... below energy threshold: skipping!" << std::endl;

        continue;
      }
      else if (shower_part._daughters.size() < fMinNumDaughters) {
        if (fDebugMode)
          std::cout << " ... below # daughter particle count threshold: skipping!" << std::endl;

        continue;
      }
      else if (fDebugMode) {
        std::cout << " ... condition matched. Storing this MCShower..." << std::endl;
      }

      // Record this MCShower
      mcs_to_spart_v.push_back(shower_index);

      if (fDebugMode)

        std::cout << " Storage index " << mcshower.size() << " => Shower index " << shower_index
                  << std::endl;

      ::sim::MCShower shower_prof;

      shower_prof.Origin(shower_part._origin);
      shower_prof.PdgCode(shower_part._pdgcode);
      shower_prof.TrackID(shower_part._track_id);
      shower_prof.Process(shower_part._process);

      shower_prof.MotherPdgCode(mother_part._pdgcode);
      shower_prof.MotherTrackID(mother_part._track_id);
      shower_prof.MotherProcess(mother_part._process);

      shower_prof.AncestorPdgCode(ancestor_part._pdgcode);
      shower_prof.AncestorTrackID(ancestor_part._track_id);
      shower_prof.AncestorProcess(ancestor_part._process);

      shower_prof.Start(MCStep(shower_part._start_vtx, shower_part._start_mom));
      shower_prof.End(MCStep(shower_part._end_vtx, shower_part._end_mom));
      shower_prof.MotherStart(MCStep(mother_part._start_vtx, mother_part._start_mom));
      shower_prof.MotherEnd(MCStep(mother_part._end_vtx, mother_part._end_mom));
      shower_prof.AncestorStart(MCStep(ancestor_part._start_vtx, ancestor_part._start_mom));
      shower_prof.AncestorEnd(MCStep(ancestor_part._end_vtx, ancestor_part._end_mom));

      // Daughter list
      std::vector<unsigned int> daughter_track_id;
      daughter_track_id.reserve(fPartAlg.ShowerDaughters(shower_index).size());

      for (auto const& index : fPartAlg.ShowerDaughters(shower_index))

        daughter_track_id.push_back(part_v.at(index)._track_id);

      shower_prof.DaughterTrackID(daughter_track_id);

      if (!daughter_stored && daughter_track_id.size() > 1) daughter_stored = true;

      mcshower.push_back(shower_prof);
    }

    if (fDebugMode)
      std::cout << " Found " << mcshower.size()
                << " MCShowers. Now computing DetProfile position..." << std::endl;

    //
    // Daughter vtx
    //
    std::vector<TLorentzVector> mcs_daughter_vtx_v(
      mcshower.size(),
      TLorentzVector(
        sim::kINVALID_DOUBLE, sim::kINVALID_DOUBLE, sim::kINVALID_DOUBLE, sim::kINVALID_DOUBLE));
    std::vector<TLorentzVector> mcs_daughter_mom_v(mcshower.size(), TLorentzVector());

    std::vector<std::vector<double>> plane_charge_v(mcshower.size(), std::vector<double>(3, 0));
    std::vector<std::vector<double>> plane_dqdx_v(mcshower.size(), std::vector<double>(3, 0));

    //For dEdx Calculation
    std::vector<double> mcs_daughter_dedx_v(mcshower.size(), 0);
    std::vector<double> mcs_daughter_dedxRAD_v(mcshower.size(), 0);
    std::vector<TVector3> mcs_daughter_dir_v(mcshower.size(), TVector3());

    for (size_t mcs_index = 0; mcs_index < mcshower.size(); ++mcs_index) {

      auto& mcs_daughter_vtx = mcs_daughter_vtx_v[mcs_index];
      auto& mcs_daughter_mom = mcs_daughter_mom_v[mcs_index];
      auto& mcs_daughter_dedx = mcs_daughter_dedx_v[mcs_index];
      auto& mcs_daughter_dedxRAD = mcs_daughter_dedxRAD_v[mcs_index];
      auto& mcs_daughter_dir = mcs_daughter_dir_v[mcs_index];
      auto& plane_charge = plane_charge_v[mcs_index];
      auto& plane_dqdx = plane_dqdx_v[mcs_index];

      for (auto const& daughter_trk_id : mcshower[mcs_index].DaughterTrackID()) {

        auto const daughter_part_index = part_v.TrackToParticleIndex(daughter_trk_id);

        auto const& daughter_part = part_v[daughter_part_index];

        auto const daughter_edep_index = edep_v.TrackToEdepIndex(daughter_trk_id);

        if (daughter_edep_index < 0) continue;

        auto const& daughter_edep = edep_v.GetEdepArrayAt(daughter_edep_index);

        if (!(daughter_edep.size())) continue;

        // Record first daughter's vtx point
        double min_dist = sim::kINVALID_DOUBLE;
        for (auto const& edep : daughter_edep) {

          double dist = std::hypot(edep.pos.X() - daughter_part._start_vtx[0],
                                   edep.pos.Y() - daughter_part._start_vtx[1],
                                   edep.pos.Z() - daughter_part._start_vtx[2]);

          if (dist < min_dist) {
            min_dist = dist;
            mcs_daughter_vtx[0] = edep.pos.X();
            mcs_daughter_vtx[1] = edep.pos.Y();
            mcs_daughter_vtx[2] = edep.pos.Z();
            mcs_daughter_vtx[3] = (dist / 100. / 2.998e8) * 1.e9 + daughter_part._start_vtx[3];
          }
        }
        if (!daughter_stored) {
          // If daughter is not stored, and shower id energetic enough, attempt to include angle info
          std::vector<double> shower_dir(3, 0);
          shower_dir[0] = mcshower[mcs_index].Start().Px();
          shower_dir[1] = mcshower[mcs_index].Start().Py();
          shower_dir[2] = mcshower[mcs_index].Start().Pz();
          double magnitude = 0;
          for (size_t i = 0; i < 3; ++i)
            magnitude += shower_dir[i] * shower_dir[i];

          magnitude = sqrt(magnitude);

          if (magnitude > 1.e-10) {
            // If enough momentum, include angle info
            min_dist = sim::kINVALID_DOUBLE;

            for (auto& v : shower_dir)
              v /= magnitude;

            for (auto const& edep : daughter_edep) {
              std::vector<double> shower_dep_dir(3, 0);
              shower_dep_dir[0] = edep.pos.X() - mcshower[mcs_index].Start().X();
              shower_dep_dir[1] = edep.pos.Y() - mcshower[mcs_index].Start().Y();
              shower_dep_dir[2] = edep.pos.Z() - mcshower[mcs_index].Start().Z();

              double const dist =
                std::hypot(shower_dep_dir[0], shower_dep_dir[1], shower_dep_dir[2]);
              if (dist == 0) continue; // it would yield nonsense anyway
              for (auto& v : shower_dep_dir)
                v /= dist;

              double angle =
                acos(shower_dep_dir[0] * shower_dir[0] + shower_dep_dir[1] * shower_dir[1] +
                     shower_dep_dir[2] * shower_dir[2]) /
                TMath::Pi() * 180.;

              if (dist < min_dist && angle < 10) {

                min_dist = dist;
                mcs_daughter_vtx[0] = edep.pos.X();
                mcs_daughter_vtx[1] = edep.pos.Y();
                mcs_daughter_vtx[2] = edep.pos.Z();
                mcs_daughter_vtx[3] =
                  (dist / 100. / 2.998e8) * 1.e9 + mcshower[mcs_index].Start().T();
              }
            }
          }
        }
        break;
      }
      // Now take care of momentum & plane charge

      std::vector<double> mom(3, 0);
      for (auto const& daughter_trk_id : mcshower[mcs_index].DaughterTrackID()) {

        //auto const daughter_part_index = part_v.TrackToParticleIndex(daughter_trk_id);

        // for c2: daughter_part is unused
        //auto const& daughter_part = part_v[daughter_part_index];

        auto const daughter_edep_index = edep_v.TrackToEdepIndex(daughter_trk_id);

        if (daughter_edep_index < 0) continue;

        auto const& daughter_edep = edep_v.GetEdepArrayAt(daughter_edep_index);

        if (!(daughter_edep.size())) continue;

        //bool first=true;  // unused
        for (auto const& edep : daughter_edep) {

          // Compute unit vector to this energy deposition
          mom[0] = edep.pos.X() - mcs_daughter_vtx[0];
          mom[1] = edep.pos.Y() - mcs_daughter_vtx[1];
          mom[2] = edep.pos.Z() - mcs_daughter_vtx[2];

          // Weight by energy (momentum)
          double magnitude = std::hypot(mom[0], mom[1], mom[2]);

          double energy = 0;
          double npid = 0;
          for (auto const& pid_energy : edep.deps) {
            npid++;
            energy += pid_energy.energy;
          }
          energy /= npid;
          if (magnitude > 1.e-10) {
            mom.at(0) = mom.at(0) * energy / magnitude;
            mom.at(1) = mom.at(1) * energy / magnitude;
            mom.at(2) = mom.at(2) * energy / magnitude;
            mcs_daughter_mom[0] += mom.at(0);
            mcs_daughter_mom[1] += mom.at(1);
            mcs_daughter_mom[2] += mom.at(2);
          }

          //Determine the direction of the shower right at the start point
          double E = 0;
          double N = 0;
          if (std::hypot(edep.pos.X() - mcs_daughter_vtx[0],
                         edep.pos.Y() - mcs_daughter_vtx[1],
                         edep.pos.Z() - mcs_daughter_vtx[2]) < 2.4 &&
              magnitude > 1.e-10) {

            mcs_daughter_dir[0] += mom.at(0);
            mcs_daughter_dir[1] += mom.at(1);
            mcs_daughter_dir[2] += mom.at(2);
            E += energy;
            N += 1;
          }

          if (E > 0) E /= N;
          mcs_daughter_dedxRAD += E;

          mcs_daughter_mom[3] += energy;

          // Charge
          auto const pid = edep.pid;
          if (pindex.hasPlane(pid))
            plane_charge[pid.Plane] += (double)(edep.deps[pindex[pid]].charge);

        } ///Looping through the MCShower daughter's energy depositions

      } ///Looping through MCShower daughters
      mcs_daughter_dedxRAD /= 2.4;

      for (auto const& daughter_trk_id : mcshower[mcs_index].DaughterTrackID()) {

        //auto const daughter_part_index = part_v.TrackToParticleIndex(daughter_trk_id);

        // for c2: daughter_part is unused
        //auto const& daughter_part = part_v[daughter_part_index];

        auto const daughter_edep_index = edep_v.TrackToEdepIndex(daughter_trk_id);

        if (daughter_edep_index < 0) continue;

        auto const& daughter_edep = edep_v.GetEdepArrayAt(daughter_edep_index);

        if (!(daughter_edep.size())) continue;

        for (auto const& edep : daughter_edep) {

          //Defining dEdx
          //Need to define a plane through the shower start point (x_0, y_0, z_0) with a normal along the momentum vector of the shower
          //The plane will be defined in the typical way:
          // a*x + b*y + c*z + d = 0
          // where, a = dir_x, b = dir_y, c = dir_z, d = - (a*x_0+b*y_0+c*z_0)
          // then the *signed* distance of any point (x_1, y_1, z_1) from this plane is:
          // D = (a*x_1 + b*y_1 + c*z_1 + d )/hypot(a,b,c)

          double p_mag = std::hypot(mcs_daughter_dir[0], mcs_daughter_dir[1], mcs_daughter_dir[2]);
          double a = 0, b = 0, c = 0, d = 0;
          if (p_mag > 1.e-10) {
            a = mcs_daughter_dir[0] / p_mag;
            b = mcs_daughter_dir[1] / p_mag;
            c = mcs_daughter_dir[2] / p_mag;
            d = -1 * (a * mcs_daughter_vtx[0] + b * mcs_daughter_vtx[1] + c * mcs_daughter_vtx[2]);
          }
          else {
            mcs_daughter_dedx += 0;
            continue;
          }
          //Radial Distance
          assert(std::hypot(a, b, c) > 0.0); // guaranteed by p_mag > 0
          double const radialDistance =
            (a * edep.pos.X() + b * edep.pos.Y() + c * edep.pos.Z() + d) / std::hypot(a, b, c);
          if (radialDistance < 2.4 && radialDistance > 0) {

            double E = 0;
            double N = 0;

            for (auto const& pid_energy : edep.deps) {
              N += 1;
              E += pid_energy.energy;
            }

            if (N > 0) { E /= N; }
            else {
              E = 0;
            }

            mcs_daughter_dedx += E;

            // Charge
            auto const pid = edep.pid;
            if (pindex.hasPlane(pid))
              plane_dqdx[pid.Plane] += (double)(edep.deps[pindex[pid]].charge);
          }
        }
      }
      mcs_daughter_dedx /= 2.4;
      plane_dqdx.at(0) /= 2.4;
      plane_dqdx.at(1) /= 2.4;
      plane_dqdx.at(2) /= 2.4;

    } ///Looping through MCShowers

    if (fDebugMode)
      std::cout << " Found " << mcshower.size() << " MCShowers. Now storing..." << std::endl;

    // Store plane charge & daughter momentum
    for (size_t mcs_index = 0; mcs_index < mcshower.size(); ++mcs_index) {

      auto& daughter_vtx = mcs_daughter_vtx_v[mcs_index];
      auto& daughter_mom = mcs_daughter_mom_v[mcs_index];
      auto& daughter_dedx = mcs_daughter_dedx_v[mcs_index];
      auto& daughter_dedxRAD = mcs_daughter_dedxRAD_v[mcs_index];
      auto& daughter_dir = mcs_daughter_dir_v[mcs_index];
      auto& plane_charge = plane_charge_v[mcs_index];
      auto& plane_dqdx = plane_dqdx_v[mcs_index];

      double magnitude = daughter_mom.P();

      if (daughter_mom[3] > 1.e-10) {
        daughter_mom[0] *= daughter_mom[3] / magnitude;
        daughter_mom[1] *= daughter_mom[3] / magnitude;
        daughter_mom[2] *= daughter_mom[3] / magnitude;
      }
      else
        for (size_t i = 0; i < 4; ++i)
          daughter_mom[i] = 0;

      mcshower.at(mcs_index).DetProfile(MCStep(daughter_vtx, daughter_mom));
      mcshower.at(mcs_index).Charge(plane_charge);
      mcshower.at(mcs_index).dQdx(plane_dqdx);
      mcshower.at(mcs_index).dEdx(daughter_dedx);
      mcshower.at(mcs_index).dEdxRAD(daughter_dedxRAD);
      mcshower.at(mcs_index).StartDir(daughter_dir);
    }

    if (fDebugMode) {

      for (auto const& prof : mcshower) {

        std::cout

          << Form("  Shower particle:     PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum "
                  "(%g,%g,%g,%g)",
                  prof.PdgCode(),
                  prof.TrackID(),
                  prof.Start().X(),
                  prof.Start().Y(),
                  prof.Start().Z(),
                  prof.Start().T(),
                  prof.Start().Px(),
                  prof.Start().Py(),
                  prof.Start().Pz(),
                  prof.Start().E())
          << std::endl
          << Form("    Mother particle:   PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum "
                  "(%g,%g,%g,%g)",
                  prof.MotherPdgCode(),
                  prof.MotherTrackID(),
                  prof.MotherStart().X(),
                  prof.MotherStart().Y(),
                  prof.MotherStart().Z(),
                  prof.MotherStart().T(),
                  prof.MotherStart().Px(),
                  prof.MotherStart().Py(),
                  prof.MotherStart().Pz(),
                  prof.MotherStart().E())
          << std::endl
          << Form("    Ancestor particle: PDG=%d : Track ID=%d Start @ (%g,%g,%g,%g) with Momentum "
                  "(%g,%g,%g,%g)",
                  prof.AncestorPdgCode(),
                  prof.AncestorTrackID(),
                  prof.AncestorStart().X(),
                  prof.AncestorStart().Y(),
                  prof.AncestorStart().Z(),
                  prof.AncestorStart().T(),
                  prof.AncestorStart().Px(),
                  prof.AncestorStart().Py(),
                  prof.AncestorStart().Pz(),
                  prof.AncestorStart().E())
          << std::endl
          << Form("    ... with %zu daughters: Start @ (%g,%g,%g,%g) with Momentum (%g,%g,%g,%g)",
                  prof.DaughterTrackID().size(),
                  prof.DetProfile().X(),
                  prof.DetProfile().Y(),
                  prof.DetProfile().Z(),
                  prof.DetProfile().T(),
                  prof.DetProfile().Px(),
                  prof.DetProfile().Py(),
                  prof.DetProfile().Pz(),
                  prof.DetProfile().E())
          << std::endl
          << "    Charge per plane: ";
        size_t const nPlanes = prof.Charge().size();
        for (size_t i = 0; i < nPlanes; ++i) {

          std::cout << " | Plane " << i << std::flush;
          std::cout << " ... Q = " << prof.Charge(i) << std::flush;
        }
        std::cout << std::endl << std::endl;
      }
    }
    return result;
  }
}
