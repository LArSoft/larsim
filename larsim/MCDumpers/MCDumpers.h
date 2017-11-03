/**
 * @file   larsim/MCDumpers/MCDumpers.h
 * @brief  Utility functions to print MC truth information.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 2, 2017
 * 
 * Functions dumping Monte Carlo data product objects (`sim::dump` namespace).
 * 
 */

#ifndef LARSIM_MCDUMPERS_MCDUMPERS_H
#define LARSIM_MCDUMPERS_MCDUMPERS_H

// LArSoft libraries
#include "larsim/MCDumpers/MCDumperUtils.h"
#include "larcorealg/CoreUtils/DumpUtils.h" // lar::dump namespace

// nutools libraries
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// ROOT libraries
#include "TLorentzVector.h"

// C/C++ standard libraries
#include <string>
#include <utility> // std::forward()


namespace sim {
  
  /// Functions to dump Monte Carlo object information into a stream.
  namespace dump {
    
    //--------------------------------------------------------------------------
    //@{
    /**
     * @brief Dumps the content of the specified particle in the output stream.
     * @tparam Stream the type of output stream
     * @param out the output stream
     * @param particle the particle to be dumped
     * @param indent base indentation string (default: none)
     * @param firstIndent indentation of first line (default: as `indent`)
     * 
     * The `indent` string is prepended to every line of output, except for
     * the first one, where `firstIndent` is used.
     * 
     * The output starts on the current line, and the last line is NOT broken.
     */
    template <typename Stream>
    void DumpMCParticle(
      Stream&& out, simb::MCParticle const& particle,
      std::string indent, std::string firstIndent
      );
    
    template <typename Stream>
    void DumpMCParticle
      (Stream&& out, simb::MCParticle const& particle, std::string indent = "")
      { DumpMCParticle(std::forward<Stream>(out), particle, indent, indent); }
    
    //@}
    
    
    //--------------------------------------------------------------------------
    //@{
    /**
     * @brief Dumps the specified particle trajectory into the output stream.
     * @tparam Stream the type of output stream
     * @param out the output stream
     * @param trajectory the particle trajectory to be dumped
     * @param pointsPerLine number of points dumped per line (default: all)
     * @param indent base indentation string (default: none)
     * 
     * All points of the specified Monte Carlo `trajectory` are printed
     * on screen, `pointsPerLine` on each line.
     * The points are printed starting on a new line, and each line is applied
     * the same indentation string (`indent`).
     * As an exception, if `pointsPerLine` is not specified, all points are
     * printed on the current line and `indent` is ignored.
     * 
     * The last line of the output is NOT broken.
     */
    template <typename Stream>
    void DumpMCParticleTrajectory(
      Stream&& out, simb::MCTrajectory const& trajectory,
      unsigned int pointsPerLine, std::string indent
      );
    
    template <typename Stream>
    void DumpMCParticleTrajectory
      (Stream&& out, simb::MCTrajectory const& trajectory)
      {
        DumpMCParticleTrajectory(std::forward<Stream>(out), trajectory, 0U, "");
      }
    
    // @}
    
    
    //--------------------------------------------------------------------------
    // @{
    /**
     * @brief Dumps the content of the specified neutrino in the output stream.
     * @tparam Stream the type of output stream
     * @param out the output stream
     * @param neutrino the neutrino to be dumped
     * @param indent base indentation string (default: none)
     * @param firstIndent string to be used for indentation of the first line
     *                    (default: as `indent`)
     * 
     * The `indent` string is prepended to every line of output, except for
     * the first one, where `firstIndent` is used.
     * 
     * The output starts on the current line, and the last line is NOT broken.
     */
    template <typename Stream>
    void DumpMCNeutrino(
      Stream&& out, simb::MCNeutrino const& neutrino,
      std::string indent, std::string firstIndent
      );
    
    template <typename Stream>
    void DumpMCNeutrino
      (Stream&& out, simb::MCNeutrino const& neutrino, std::string indent = "")
      { DumpMCNeutrino(std::forward<Stream>(out), neutrino, indent, indent); }
    
    // @}
    
    
    //--------------------------------------------------------------------------
    // @{
    /**
     * @brief Dumps the content of the specified MC truth in the output stream.
     * @tparam Stream the type of output stream
     * @param out the output stream
     * @param truth the truth information to be dumped
     * @param pointsPerLine (_optional_) number of trajectory points dumped per
     *                      line (default: `0`, no points shown)
     * @param indent base indentation string (default: none)
     * @param firstIndent string to be used for indentation of the first line
     *                    (default: same as `indent`)
     * 
     * The `indent` string is prepended to every line of output, except for
     * the first one, where `firstIndent` is used.
     * 
     * The argument `pointsPerLine` regulates the dump of trajectory points from
     * all the particles in the record (except the ones stored in the neutrino
     * object). Setting it to `0`, or leaving it out, will suppress the dump of
     * particle trajectories completely. There is no option to reduce the number
     * of printed trajectory points: it's just all or none.
     * 
     * The output starts on the current line, and the last line is NOT broken.
     */
    template <typename Stream>
    void DumpMCTruth(
      Stream&& out, simb::MCTruth const& truth, unsigned int pointsPerLine,
      std::string indent, std::string firstIndent
      );
    
    template <typename Stream>
    void DumpMCTruth(
      Stream&& out, simb::MCTruth const& truth, unsigned int pointsPerLine,
      std::string indent = ""
      )
      {
        DumpMCTruth
          (std::forward<Stream>(out), truth, pointsPerLine, indent, indent);
      }
    
    template <typename Stream>
    void DumpMCTruth(
      Stream&& out, simb::MCTruth const& truth,
      std::string indent, std::string firstIndent
      )
      { DumpMCTruth(std::forward<Stream>(out), truth, 0, indent, firstIndent); }
    
    template <typename Stream>
    void DumpMCTruth(
      Stream&& out, simb::MCTruth const& truth, std::string indent = ""
      )
      { DumpMCTruth(std::forward<Stream>(out), truth, indent, indent); }
    
    // @}
    
    
    //--------------------------------------------------------------------------
    // @{
    /**
     * @brief Dumps the content of the GENIE truth in the output stream.
     * @tparam Stream the type of output stream
     * @param out the output stream
     * @param truth the truth information to be dumped
     * @param indent base indentation string (default: none)
     * @param firstIndent string to be used for indentation of the first line
     *                    (default: as `indent`)
     * 
     * The `indent` string is prepended to every line of output, except for
     * the first one, where `firstIndent` is used.
     * 
     * The output starts on the current line, and the last line is NOT broken.
     */
    template <typename Stream>
    void DumpGTruth(
      Stream&& out, simb::GTruth const& truth,
      std::string indent, std::string firstIndent
      );
    
    template <typename Stream>
    void DumpGTruth
      (Stream&& out, simb::GTruth const& truth, std::string indent = "")
      { DumpGTruth(std::forward<Stream>(out), truth, indent, indent); }
    
    // @}
    
    //--------------------------------------------------------------------------
    template <typename Stream, typename Vector>
    void DumpLorentzVector(Stream&& out, Vector const& v);
    
    
    template <typename Stream>
    Stream& operator<< (Stream&& out, TLorentzVector const& v)
      { DumpLorentzVector(std::forward<Stream>(out), v); return out; }
    
    //--------------------------------------------------------------------------
    
  } // namespace dumpers
  
} // namespace sim



//------------------------------------------------------------------------------
//---  template implementation
//------------------------------------------------------------------------------
template <typename Stream, typename Vector>
void sim::dump::DumpLorentzVector(Stream&& out, Vector const& v) {
  out
    << "(" << v.X() << ", " << v.Y() << ", " << v.Z() << "; " << v.T() << ")";
} // sim::dump::DumpLorentzVector(Vector)


//------------------------------------------------------------------------------
template <typename Stream>
void sim::dump::DumpMCParticle(
  Stream&& out, simb::MCParticle const& particle,
  std::string indent, std::string firstIndent
) {
  out << firstIndent
    << "ID=" << particle.TrackId() << ": " << ParticleName(particle.PdgCode())
    << " mass=" << particle.Mass() << " GeV/c2 "
    << " status=" << particle.StatusCode()
      << " (" << ParticleStatusName(particle.StatusCode()) << ")"
    ;
  if (particle.Weight() != 1.0) out << " weight=" << particle.Weight();
  if (particle.Rescatter()) {
    out << " rescattered (" << particle.Rescatter()
      << ") at vertex " << particle.GetGvtx();
  }
  out << "\n" << indent << "created via "
    << (particle.Process().empty()? "magics": particle.Process());
  if (particle.Mother() == 0) out << " by the gods";
  else                        out << " from ID=" << particle.Mother();
  
  const unsigned int nDaughters = particle.NumberDaughters();
  const unsigned int nPoints = particle.NumberTrajectoryPoints();
  if (nPoints > 0) {
    TLorentzVector const& start = particle.Position();
    TLorentzVector const& start_mom = particle.Momentum();
    out << " at " << start << " cm with momentum " << start_mom << " GeV/c";
  }
  if (particle.Polarization().Mag2() != 0.) {
    out
      << " with polarization " << lar::dump::vector3D(particle.Polarization());
  }
  // does this particle seem to stop? (by decay, or because of extended traj.)
  if ((nPoints > 1) || (nDaughters > 0)) {
    out << " and " << ((nDaughters > 0)? "ends": "stops") << " by "
      << (particle.EndProcess().empty()? "magics": particle.EndProcess());
    if (nPoints > 1) {
      TLorentzVector const& stop = particle.EndPosition();
      TLorentzVector const& stop_mom = particle.EndMomentum();
      out << " at " << stop << " cm with momentum " << stop_mom << " GeV/c";
    }
    if (nDaughters > 0) {
      out << " into ";
      if (nDaughters == 1)
        out << "particle ID=" << particle.FirstDaughter();
      else {
        out << nDaughters << " particles from ID=" << particle.FirstDaughter()
          << " to ID=" << particle.LastDaughter();
      }
    } // if daughters
  } // if it stops
  if (nPoints > 1) {
    simb::MCTrajectory const& traj = particle.Trajectory();
    out << "\n" << indent << "comes with a trajectory " << traj.TotalLength()
      << " cm long in " << nPoints << " points";
  } // if has points
 
} // sim::dump::DumpMCParticle()


//------------------------------------------------------------------------------
template <typename Stream>
void sim::dump::DumpMCParticleTrajectory(
  Stream&& out, simb::MCTrajectory const& trajectory,
  unsigned int pointsPerLine, std::string indent
) {
  unsigned int page = 0;
  for (auto const& pair: trajectory) {
    if ((pointsPerLine > 0) && (page-- == 0)) {
      out << "\n" << indent << "  ";
      page = pointsPerLine - 1;
    }
    else out << " -- ";
    
    TLorentzVector const& pos = pair.first;
    out << pos;
  } // for trajectory points
  
} // sim::dump::DumpMCParticleTrajectory()


//------------------------------------------------------------------------------
template <typename Stream>
void sim::dump::DumpMCNeutrino(
  Stream&& out, simb::MCNeutrino const& nu,
  std::string indent, std::string firstIndent
) {

  out << firstIndent
    << "from " << TruthCCNCname(nu.CCNC())
      << ", " << TruthInteractionTypeName(nu.InteractionType())
      << ", mode: " << nu.Mode() << " (" << TruthReactionMode(nu.Mode()) << ")"
    << '\n' << indent
      << "target: " << nu.Target() << " (" << ParticleName(nu.Target()) << ")"
    ;
  if (nu.HitNuc() != 0) {
    out << ", hit nucleon: " << nu.HitNuc()
      << " (" << ParticleName(nu.HitNuc()) << ")";
  }
  if (nu.HitQuark() != 0) {
    out << ", hit quark: " << nu.HitQuark()
      << " (" << ParticleName(nu.HitQuark()) << ")";
  }
  out
    << '\n' << indent
      << "x=" << nu.X() << " y=" << nu.Y() << " w=" << nu.W()
      << " Q^2=" << nu.QSqr() << " GeV^2; theta=" << nu.Theta()
      << " rad pT=" << nu.Pt() << " GeV/c"
    ;
  out << '\n' << indent << "neutrino: ";
  DumpMCParticle(std::forward<Stream>(out), nu.Nu(), indent + "  ", "");
  out << '\n' << indent << "outgoing lepton: ";
  DumpMCParticle(std::forward<Stream>(out), nu.Lepton(), indent + "  ", "");
  
} // sim::dump::DumpMCNeutrino()


//------------------------------------------------------------------------------
template <typename Stream>
void sim::dump::DumpMCTruth(
  Stream&& out, simb::MCTruth const& truth, unsigned int pointsPerLine,
  std::string indent, std::string firstIndent
) {
  unsigned int const nParticles = truth.NParticles();
  out << firstIndent
    << nParticles << " particles from "
    << TruthOriginName(truth.Origin());
  if (truth.NeutrinoSet()) {
    out << '\n' << indent << "neutrino information: ";
    DumpMCNeutrino
      (std::forward<Stream>(out), truth.GetNeutrino(), indent + "  ", "");
  }
  for (unsigned int i = 0; i < nParticles; ++i) {
    out << '\n' << indent << "[#" << i << "] ";
    simb::MCParticle const& particle = truth.GetParticle(i);
    DumpMCParticle(std::forward<Stream>(out), particle, indent + "  ", "");
    
    const unsigned int nPoints = particle.NumberTrajectoryPoints();
    if ((nPoints > 0) && (pointsPerLine > 0)) {
      out << ":";
      DumpMCParticleTrajectory(
        std::forward<Stream>(out), particle.Trajectory(),
        pointsPerLine, indent + "    "
        );
    } // if has points
  } // for all particles
  
} // sim::dump::DumpMCTruth()


//------------------------------------------------------------------------------
template <typename Stream>
void sim::dump::DumpGTruth(
  Stream&& out, simb::GTruth const& truth,
  std::string indent, std::string firstIndent
) {
  
  unsigned int const nCharged
    = truth.fNumPiPlus + truth.fNumPiMinus + truth.fNumProton;
  unsigned int const nNeutral = truth.fNumPi0 + truth.fNumNeutron;
  unsigned int const nPions
    = truth.fNumPiPlus + truth.fNumPiMinus + truth.fNumPi0;
  unsigned int const nNucleons = truth.fNumProton + truth.fNumNeutron;
  unsigned int const nTotalParticles = nCharged + nNeutral;
  
  out << firstIndent
      << "interaction code: " << truth.fGint
      << ", neutrino scattering code: " << truth.fGscatter
      << " at " << truth.fVertex
    << "\n" << indent 
      << "probe: " << ParticleName(truth.fProbePDG)
      << " with cp=" << truth.fProbeP4
      << " hit nucleon with cp=" << truth.fHitNucP4 << " GeV"
      << " (" << (truth.fIsSeaQuark? "": "not a ") << "sea quark)"
      << " in target: " << ParticleName(truth.ftgtPDG)
      << " (Z: " << truth.ftgtZ << ", A: " << truth.ftgtA << ")"
    << "\n" << indent 
      << "event interaction weight (genie internal): " << truth.fweight
      << ", interaction probability: " << truth.fprobability
      << ", cross section: " << truth.fXsec
      << ", differential cross section: " << truth.fDiffXsec
    << "\n" << indent
      << "particles after reaction, before FSI: "
              << truth.fNumPiPlus  << " pi+"
      << ", " << truth.fNumPiMinus << " pi-"
      << ", " << truth.fNumPi0     << " pi0"
      << ", " << truth.fNumProton  << " p/pbar"
      << ", " << truth.fNumNeutron << " n/nbar"
    << "\n" << indent
      << "  total " << nTotalParticles << " particles after reaction before FSI"
        ": " << nCharged << "/" << nNeutral << " charged/neutral"
        ", " << nPions << " pions, " << nNucleons << " nucleons"
    << "\n" << indent << "process "
      << (truth.fIsCharm? "with": "without") << " charmed hadron";
  if (truth.fResNum == -1) out << ", no resonance";
  else                     out << ", resonance: #" << truth.fResNum;
  out
    << "\n" << indent
      << "internal (on shell) genie kinematics: Q^2: " << truth.fgQ2 << " GeV^2"
      << " q^2: " << truth.fgq2 << " GeV^2"
      << ", w: " << truth.fgW << " GeV^2"
      << ", t: " << truth.fgT << " GeV^2"
      << ", x: " << truth.fgX
      << ", y: " << truth.fgY
    << "\n" << indent
      << "FShadSyst: " << truth.fFShadSystP4
    ;
  
} // sim::DumpGTruth::DumpTruth()


//------------------------------------------------------------------------------


#endif // LARSIM_MCDUMPERS_MCDUMPERS_H
