/**
 * @file   larsim/MCDumpers/MCDumpersUtils.cxx
 * @brief  Utility functions to print MC truth information (implementation).
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   November 2, 2017
 * @see    larsim/MCDumpers/MCDumpersUtils.h
 * 
 */

// library header
#include "larsim/MCDumpers/MCDumperUtils.h"

// ROOT
#include "TDatabasePDG.h"
#include "TParticlePDG.h"


//------------------------------------------------------------------------------
std::string sim::TruthOriginName(simb::Origin_t origin) {
  switch (origin) {
    case simb::kUnknown          : return "unknown origin";
    case simb::kBeamNeutrino     : return "neutrinos from beam";
    case simb::kCosmicRay        : return "cosmic rays";
    case simb::kSuperNovaNeutrino: return "supernova neutrinos";
    case simb::kSingleParticle   : return "single particle";
    default:        return "unsupported (" + std::to_string((int)origin) + ")";
  } // switch
} // sim::TruthOriginName()


//------------------------------------------------------------------------------
std::string sim::TruthCCNCname(int ccnc) {
  switch (ccnc) {
    case simb::kCC: return "charged weak current";
    case simb::kNC: return "neutral weak current";
    default:        return "unsupported (" + std::to_string(ccnc) + ")";
  } // switch
} // sim::TruthCCNCname()


//------------------------------------------------------------------------------
std::string sim::TruthReactionMode(int mode) {
  
  switch (mode) {
    case  0: return "quasi-elastic";
    case  1: return "resonant";
    case  2: return "deep inelastic";
    case  3: return "coherent";
    default: return "unknown mode";
  } // switch
  
} // sim::TruthReactionMode()


//------------------------------------------------------------------------------
std::string sim::TruthInteractionTypeName(int type) {
  switch (type) {
    case simb::kUnknownInteraction            : return "unknown interaction";
    case simb::kQE                            : return "quasi-elastic scattering";
    case simb::kRes                           : return "resonant scattering";
    case simb::kDIS                           : return "deep inelastic scattering";
    case simb::kCoh                           : return "coherent scattering";
    case simb::kCohElastic                    : return "coherent elastic scattering";
    case simb::kElectronScattering            : return "electron scattering";
    case simb::kIMDAnnihilation               : return "inverse muon decay annihilation";
    case simb::kInverseBetaDecay              : return "inverse beta decay";
    case simb::kGlashowResonance              : return "Glashow resonance";
    case simb::kAMNuGamma                     : return "anomalous neutrino-photon interaction";
    case simb::kMEC                           : return "meson exchange current";
    case simb::kDiffractive                   : return "diffractive";
    case simb::kEM                            : return "electromagnetic";
    case simb::kWeakMix                       : return "weak mixing";
    case simb::kNuanceOffset                  : return "<nuance offset>";
    case simb::kCCQE                          : return "charged current quasi-elastic scattering";
    case simb::kNCQE                          : return "neutral current quasi-elastic scattering";
    case simb::kResCCNuProtonPiPlus           : return "resonant charged current neutrino proton pi+";
    case simb::kResCCNuNeutronPi0             : return "resonant charged current neutrino neutron pi0";
    case simb::kResCCNuNeutronPiPlus          : return "resonant charged current neutrino neutron pi+";
    case simb::kResNCNuProtonPi0              : return "resonant neutral current neutrino proton pi0";
    case simb::kResNCNuProtonPiPlus           : return "resonant neutral current neutrino proton pi+";
    case simb::kResNCNuNeutronPi0             : return "resonant neutral current neutrino neutron pi0";
    case simb::kResNCNuNeutronPiMinus         : return "resonant neutral current neutrino neutron pi-";
    case simb::kResCCNuBarNeutronPiMinus      : return "resonant charged current antineutrino neutron pi-";
    case simb::kResCCNuBarProtonPi0           : return "resonant charged current antineutrino proton pi0";
    case simb::kResCCNuBarProtonPiMinus       : return "resonant charged current antineutrino proton pi-";
    case simb::kResNCNuBarProtonPi0           : return "resonant neutral current antineutrino proton pi0";
    case simb::kResNCNuBarProtonPiPlus        : return "resonant neutral current antineutrino proton pi+";
    case simb::kResNCNuBarNeutronPi0          : return "resonant neutral current antineutrino neutron pi0";
    case simb::kResNCNuBarNeutronPiMinus      : return "resonant neutral current antineutrino neutron pi-";
    case simb::kResCCNuDeltaPlusPiPlus        : return "resonant charged current neutrino Delta+ pi+";
    case simb::kResCCNuDelta2PlusPiMinus      : return "resonant charged current neutrino Delta++ pi-";
    case simb::kResCCNuBarDelta0PiMinus       : return "resonant charged current antineutrino Delta0 pi-";
    case simb::kResCCNuBarDeltaMinusPiPlus    : return "resonant charged current antineutrino Delta- pi+";
    case simb::kResCCNuProtonRhoPlus          : return "resonant charged current neutrino proton rho+";
    case simb::kResCCNuNeutronRhoPlus         : return "resonant charged current neutrino neutron rho+";
    case simb::kResCCNuBarNeutronRhoMinus     : return "resonant charged current antineutrino neutron rho-";
    case simb::kResCCNuBarNeutronRho0         : return "resonant charged current antineutrino neutron rho0";
    case simb::kResCCNuSigmaPlusKaonPlus      : return "resonant charged current neutrino Sigma+ kaon+";
    case simb::kResCCNuSigmaPlusKaon0         : return "resonant charged current neutrino Sigma+ kaon0";
    case simb::kResCCNuBarSigmaMinusKaon0     : return "resonant charged current antineutrino Sigma- kaon0";
    case simb::kResCCNuBarSigma0Kaon0         : return "resonant charged current antineutrino Sigma0 kaon0";
    case simb::kResCCNuProtonEta              : return "resonant charged current neutrino proton eta";
    case simb::kResCCNuBarNeutronEta          : return "resonant charged current antineutrino neutron eta";
    case simb::kResCCNuKaonPlusLambda0        : return "resonant charged current neutrino Kaon+ Lambda0";
    case simb::kResCCNuBarKaon0Lambda0        : return "resonant charged current antineutrino kaon0 Lambda0";
    case simb::kResCCNuProtonPiPlusPiMinus    : return "resonant charged current neutrino proton pi+ pi-";
    case simb::kResCCNuProtonPi0Pi0           : return "resonant charged current neutrino proton pi0 pi0";
    case simb::kResCCNuBarNeutronPiPlusPiMinus: return "resonant charged current antineutrino neutron pi+ pi-";
    case simb::kResCCNuBarNeutronPi0Pi0       : return "resonant charged current antineutrino neutron pi0 pi0";
    case simb::kResCCNuBarProtonPi0Pi0        : return "resonant charged current antineutrino proton pi0 pi0";
    case simb::kCCDIS                         : return "charged current deep inelastic scattering";
    case simb::kNCDIS                         : return "neutral current deep inelastic scattering";
    case simb::kUnUsed1                       : return "unused (1)";
    case simb::kUnUsed2                       : return "unused (2)";
    case simb::kCCQEHyperon                   : return "charged current quasi-elastic scattering with hyperon";
    case simb::kNCCOH                         : return "neutral current coherent scattering";
    case simb::kCCCOH                         : return "charged current coherent scattering";
    case simb::kNuElectronElastic             : return "electron neutrino elastic";
    case simb::kInverseMuDecay                : return "inverse muon decay";
    default                                   : return "unsupported (" + std::to_string(type) + ")";
  } // switch
} // sim::TruthInteractionTypeName()


//------------------------------------------------------------------------------
std::string sim::ParticleStatusName(int code) {

  switch(code) {
    case -1: return "undefined";
    case  0: return "initial state";
    case  1: return "stable final state";
    case  2: return "intermediate";
    case  3: return "decayed";
    case 11: return "nucleon target";
    case 12: return "pre-fragmentation hadronic state";
    case 13: return "pre-decay resonant state";
    case 14: return "hadron in nucleus";
    case 15: return "final state nuclear remnant";
    case 16: return "nucleon cluster target";
    default: return "unknown";
  } // switch
  
} // sim::ParticleStatusName


//------------------------------------------------------------------------------
std::string sim::ParticleName(int pigid) {
  TParticlePDG const* PDGinfo = TDatabasePDG::Instance()->GetParticle(pigid);
  return PDGinfo? PDGinfo->GetTitle(): ("PDG ID " + std::to_string(pigid));
} // sim::ParticleName()

  
//------------------------------------------------------------------------------
