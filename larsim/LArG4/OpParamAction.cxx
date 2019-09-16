// OpParamAction class implementation : Ben Jones, MIT 2013
// See header file for detailed description.

#include "larsim/LArG4/OpParamAction.h"

#include "cetlib_except/exception.h"

#include <cmath>

namespace larg4
{

  //-----------------------------------------------
  //  OpParamAction base class methods
  //-----------------------------------------------



  OpParamAction::OpParamAction()
  {
  }

  OpParamAction::~OpParamAction()
  {
  }

  double OpParamAction::GetAttenuationFraction(G4ThreeVector /*PhotonDirection*/, G4ThreeVector /*PhotonPosition*/)
  {
    return 0;
  }



  //-----------------------------------------------
  //  SimpleWireplaneAction Methods
  //-----------------------------------------------


  SimpleWireplaneAction::SimpleWireplaneAction(TVector3 WireDirection, TVector3 PlaneNormal, double WirePitch, double WireDiameter)
  {
    fWireDirection = G4ThreeVector(WireDirection.X(),
				   WireDirection.Y(),
				   WireDirection.Z()).unit();
    fPlaneNormal   = G4ThreeVector(PlaneNormal.X(),
				   PlaneNormal.Y(),
				   PlaneNormal.Z()).unit();
    fDPRatio       = WireDiameter/WirePitch;
  }

  //-----------------------------------------------

  SimpleWireplaneAction::~SimpleWireplaneAction()
  {
  }

  //-----------------------------------------------
  // An ideal simple wireplane attenuates the light by a fraction
  //  Wire_Diameter / (Pitch * cos theta) where theta is the angle
  //  of incident light projected into the plane perpendicular to the
  //  wires.  The photon position is not used.
  //
  double SimpleWireplaneAction::GetAttenuationFraction(G4ThreeVector PhotonDirection, G4ThreeVector /*PhotonPosition*/)
  {
    G4ThreeVector ProjDirection = PhotonDirection - fWireDirection*(fWireDirection.dot(PhotonDirection));
    double CosTheta = std::abs(fPlaneNormal.dot(ProjDirection));
    if(CosTheta < fDPRatio)
      return 0;
    else
      return (1.0 - fDPRatio /  CosTheta);
   }



  //-----------------------------------------------
  //  OverlaidWireplanesAction Methods
  //-----------------------------------------------


  OverlaidWireplanesAction::OverlaidWireplanesAction(std::vector<std::vector<double> > InputVectors, int Orientation)
  {

    G4ThreeVector WireBasis1, WireBasis2;

    if(Orientation==0)
      {
	fPlaneNormal=G4ThreeVector(1,0,0);
	WireBasis1 =G4ThreeVector(0,1,0);
	WireBasis2 =G4ThreeVector(0,0,1);
      }
    else if(Orientation==1)
      {
	fPlaneNormal=G4ThreeVector(0,1,0);
	WireBasis1 =G4ThreeVector(1,0,0);
	WireBasis2 =G4ThreeVector(0,0,1);
      }
    else if(Orientation==2)
      {
	fPlaneNormal=G4ThreeVector(0,0,1);
	WireBasis1 =G4ThreeVector(0,1,0);
	WireBasis2 =G4ThreeVector(0,0,1);
      }
    else
      {
	throw cet::exception("OpParamAction") << "Unrecognized wireplane orientation. Options are 1=Xdrift, 2=Ydrift, 3=Zdrift\n";
      }
    for(size_t i=0; i!=InputVectors.size(); ++i)
      {
	if(InputVectors.at(i).size()!=3)
	  {
	    throw cet::exception("OpParamAction") << "Unrecognized wireplane parameter format. Expected vector(3)'s with v[0] = wire angle, v[1] = wire pitch, v[2] = wire diameter\n";
	  }
	else
	  {
	    double theta = InputVectors[i][0]*3.142/180.;
	    fWireDirections.push_back(cos(theta)*WireBasis1 + sin(theta)*WireBasis2);
	    fDPRatios.push_back(InputVectors[i][2]/InputVectors[i][1]);
	  }
      }
  }


  //-----------------------------------------------

  OverlaidWireplanesAction::~OverlaidWireplanesAction()
  {
  }


  //-----------------------------------------------

  double OverlaidWireplanesAction::GetAttenuationFraction(G4ThreeVector PhotonDirection, G4ThreeVector /*PhotonPosition*/)
  {

    double AttenFraction=1.;

    for(size_t i=0; i!=fWireDirections.size(); ++i)
      {
	G4ThreeVector ProjDirection = PhotonDirection - fWireDirections.at(i) * (fWireDirections.at(i).dot(PhotonDirection.unit()));

	  // fWireDirections.at(i).cross(PhotonDirection.cross(fWireDirections.at(i)).unit());
	double CosTheta =  (ProjDirection.mag() > 0 ? std::abs(fPlaneNormal.dot(ProjDirection))/ProjDirection.mag() : 1.0);
	if(CosTheta < fDPRatios.at(i))
	  return 0;
	else
	  AttenFraction *= (1.0 - fDPRatios.at(i) / CosTheta);
      }
    return AttenFraction;
  }





}
