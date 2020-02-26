// OpParamAction.h - Ben Jones, MIT 2013
//
// This header file defines various optically parameterized volume actions.
//
// These objects are attached to OpParamSD objects in the OpDetReadoutGeometry class.
// When a photon steps into a volume with such a sensitive detector object attached,
// the GetAttenuationFraction method of one of the classes in this file
// is called to provide a survival probability for this photon based on its
// position and momentum
//
// In this way, generically opaque surfaces with different position and
// directionally dependent attenuation coefficients can be easily implemented
// into LArSoft and attached to physical volumes in the detector geometry
// by adding new implementations of the OpParamAction class.
//
// Two simple examples are provided:
//
// - SimpleWireplaneAction represents the angular tranission coefficient of
//     a single wireplane in 3D
//
// - OverlaidWireplaneAction represents the angular transmission coefficient
//    of multiple overlaid wireplanes in 3D.  This is the implementation used
//    to model the optical transmission of wireplanes in MicroBooNE.
//

#ifndef OPPARAMACTION_H
#define OPPARAMACTION_H

#include "TVector3.h"
#include "Geant4/G4ThreeVector.hh"

#include <vector>

namespace larg4
{



  //---------------------------------------------------
  // Abstract base class
  //---------------------------------------------------

  class OpParamAction
  {
  public:
    OpParamAction();
    virtual ~OpParamAction();
    virtual double GetAttenuationFraction(G4ThreeVector PhotonDirection, G4ThreeVector PhotonPosition);

  private:

  };


  //---------------------------------------------------
  // TransparentPlaneAction class
  //---------------------------------------------------

  class TransparentPlaneAction: public OpParamAction
  {
  public:
    TransparentPlaneAction() {};
    ~TransparentPlaneAction() {};
    double GetAttenuationFraction(G4ThreeVector /*PhotonDirection*/, G4ThreeVector /*PhotonPosition*/) {return 1;}

  private:

  };


  //---------------------------------------------------
  // SimpleWireplaneAction class
  //---------------------------------------------------

  class SimpleWireplaneAction: public OpParamAction
  {
  public:
    SimpleWireplaneAction(TVector3 WireDirection, TVector3 PlaneNormal, double WirePitch, double WireDiameter)  ;
    SimpleWireplaneAction(std::vector<std::vector<double> >, int);
    ~SimpleWireplaneAction() ;

    double GetAttenuationFraction(G4ThreeVector PhotonDirection, G4ThreeVector PhotonPosition);

  private:
    G4ThreeVector fWireDirection;
    G4ThreeVector fPlaneNormal;
    double   fDPRatio;
  };


  //---------------------------------------------------
  // OverlaidWireplanesAction class
  //---------------------------------------------------

  class OverlaidWireplanesAction: public OpParamAction
  {
  public:
    OverlaidWireplanesAction(std::vector<std::vector<double> >, int);
    ~OverlaidWireplanesAction() ;

    double GetAttenuationFraction(G4ThreeVector PhotonDirection, G4ThreeVector PhotonPosition);

  private:

    G4ThreeVector fPlaneNormal;
    std::vector<G4ThreeVector> fWireDirections;
    std::vector<double> fDPRatios;
  };


}

#endif
