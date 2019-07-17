////////////////////////////////////////////////////////////////////////
/// \file OpDetReadoutGeometry.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
#include "nug4/G4Base/DetectorConstruction.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Types.hh"
#include "Geant4/G4VPhysicalVolume.hh"

#include "larsim/Simulation/LArG4Parameters.h"
#include "larsim/LArG4/OpDetReadoutGeometry.h"
#include "larsim/LArG4/OpDetLookup.h"
#include "larsim/LArG4/OpDetSensitiveDetector.h"
#include "larsim/LArG4/OpParamSD.h"

#include <sstream>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {

  OpDetReadoutGeometry::OpDetReadoutGeometry(G4String OpDetSensitiveName, const G4String name) :
    G4VUserParallelWorld(name)
  {
    fOpDetSensitiveName = OpDetSensitiveName;
  }

  OpDetReadoutGeometry::~OpDetReadoutGeometry()
  {}

  void OpDetReadoutGeometry::Construct()
  {
    mf::LogInfo("OpDetReadoutGeometry") << "constructing parallel world, looking for "
				      << fOpDetSensitiveName;

    // Get an empty parallel world volume
    G4VPhysicalVolume * ParallelWorld = GetWorld();

    // Start with empty vectors
    std::vector<G4LogicalVolume*>  OpDetVolumes;
    std::vector<G4Transform3D>     OpDetTransformations;

    // Get the primary world volume
    G4VPhysicalVolume * WorldPhysical = g4b::DetectorConstruction::GetWorld();

    // Find the OpDet volumes
    std::vector<G4Transform3D> EmptyVector;
    EmptyVector.clear();
    FindVolumes(WorldPhysical, fOpDetSensitiveName, EmptyVector, OpDetVolumes, OpDetTransformations );

    fOpDetVolumes          = OpDetVolumes;
    fOpDetTransformations  = OpDetTransformations;

    // Get the OpDet Lookup Table
    OpDetLookup * TheOpDetLookup = OpDetLookup::Instance();

    // Create sensitive detector
    OpDetSensitiveDetector * TheSD = new OpDetSensitiveDetector("OpDetSensitiveDetector");


    if(OpDetVolumes.size()>0)
      {



	// Make placements
	for(unsigned int i=0; i!=OpDetVolumes.size(); i++)
	{
	  std::stringstream VolumeName;
	  VolumeName.flush();
	  VolumeName.str("OpDetVolume_");
	  VolumeName<<i;

	  G4Transform3D      TheTransform = OpDetTransformations.at(i);

	  G4VSolid * TheSolid = OpDetVolumes.at(i)->GetSolid();
	  G4Material * TheMaterial = OpDetVolumes.at(i)->GetMaterial();
	  G4LogicalVolume * TheLogVolume = new G4LogicalVolume(TheSolid,TheMaterial,VolumeName.str().c_str());

	  TheLogVolume->SetSensitiveDetector(TheSD   );


	  G4PVPlacement * ThePlacement
	    = new G4PVPlacement( TheTransform,
				 VolumeName.str().c_str(),
				 TheLogVolume,
				 ParallelWorld,
				 false,
				 0);

	  //	  CLHEP::Hep3Vector trans = ThePlacement->GetTranslation();

	  TheOpDetLookup->AddPhysicalVolume(ThePlacement);

	}
      }


    // Now add any optically parameterized volumes

    std::vector<G4LogicalVolume*>  OpParamVolumesFound;
    std::vector<G4Transform3D>     OpParamTransformationsFound;

    art::ServiceHandle<sim::LArG4Parameters const> lgp;

    std::vector<std::string> OpParamModels       = lgp->OpticalParamModels();
    std::vector<std::string> OpParamVolumes      = lgp->OpticalParamVolumes();
    std::vector<int>         OpParamOrientations = lgp->OpticalParamOrientations();;
    std::vector<std::vector<std::vector<double > > > OpParamParameters    = lgp->OpticalParamParameters();

    if((OpParamModels.size()!=OpParamVolumes.size())||
       (OpParamModels.size()!=OpParamOrientations.size())||
       (OpParamModels.size()!=OpParamParameters.size()))
      {
        throw cet::exception("OpDetReadoutGeometry")<<"sizes of OpParam specification vectors do not match\n";
      }

    for(size_t imodel=0; imodel!=OpParamVolumes.size(); ++imodel)
      {
	EmptyVector.clear();
	FindVolumes(WorldPhysical, OpParamVolumes.at(imodel), EmptyVector, OpParamVolumesFound, OpParamTransformationsFound );
	mf::LogInfo("OpDetReadoutGeometry")<< "Found " << OpParamVolumesFound.size()<< " volumes of name " << OpParamVolumes.at(imodel)<< " to attach optical parameterization " << OpParamModels.at(imodel)<<std::endl;


	// Since the same named model may be instantiated more than once with
	//  different parameters,  create a unique sensitive detector name for this
	//  instance

	std::stringstream SDName("");
	SDName<<OpParamModels.at(imodel)<<"_"<<imodel;

	OpParamSD * ParamSD = new OpParamSD(SDName.str().c_str(), OpParamModels.at(imodel), OpParamOrientations.at(imodel), OpParamParameters.at(imodel));


	if(OpParamVolumesFound.size()>0)
	  {
	    for(unsigned int ivol=0; ivol!=OpParamVolumes.size(); ivol++)
	      {
		std::stringstream VolumeName;
		VolumeName.flush();
		VolumeName.str("OpParamVolume_");
		VolumeName<<ivol;

		G4Transform3D      TheTransform = OpParamTransformationsFound.at(ivol);

		G4VSolid * TheSolid = OpParamVolumesFound.at(ivol)->GetSolid();
		G4Material * TheMaterial = OpParamVolumesFound.at(ivol)->GetMaterial();
		G4LogicalVolume * TheLogVolume = new G4LogicalVolume(TheSolid,TheMaterial,VolumeName.str().c_str());


		G4PVPlacement * ThePlacement
		  = new G4PVPlacement( TheTransform,
				       VolumeName.str().c_str(),
				       TheLogVolume,
				       ParallelWorld,
				       false,
				       0);

		TheLogVolume->SetSensitiveDetector(ParamSD   );


		// This line just suppressed a compiler warning
		ThePlacement->GetTranslation();


	      }

	  }
      }
  }



  void OpDetReadoutGeometry::FindVolumes(G4VPhysicalVolume * PhysicalVolume, G4String OpDetName, std::vector<G4Transform3D> TransformSoFar, std::vector<G4LogicalVolume*>& OpDetVolumes, std::vector<G4Transform3D>& OpDetTransformations)
  {

    // Add the next layer of transformation to the vector
    G4ThreeVector Translation = PhysicalVolume->GetObjectTranslation();
    G4RotationMatrix Rotation = PhysicalVolume->GetObjectRotationValue();
    G4Transform3D NextTransform( Rotation, Translation );

    TransformSoFar.push_back(NextTransform);


    // Check if this volume is a OpDet
    G4String OpDetNameUnderscore = OpDetName+"_";
    G4String VolumeName = PhysicalVolume->GetName();
    if( ( VolumeName == OpDetName ) ||
	( VolumeName.find( OpDetNameUnderscore,0,OpDetNameUnderscore.length() )==0 )
	)
      {

	// We found a OpDet! Store its volume and global transformation
	G4ThreeVector     Trans(0,0,0);
	G4RotationMatrix  Rot(0,0,0);
	G4Transform3D TotalTransform(Rot,Trans);
	//for ( std::vector<G4Transform3D>::reverse_iterator it = TransformSoFar.rbegin();
	//     it != TransformSoFar.rend(); ++it )



	for ( std::vector<G4Transform3D>::iterator it = TransformSoFar.begin();
	      it != TransformSoFar.end(); ++it )
	  {
	    CLHEP::Hep3Vector trans = (*it).getTranslation();
	    CLHEP::HepRotation rot =  (*it).getRotation();
	    TotalTransform =  TotalTransform * (*it);
	  }

	OpDetVolumes.push_back(PhysicalVolume->GetLogicalVolume());
	OpDetTransformations.push_back(TotalTransform);

      }
    else
      {
	// We did not find a OpDet.  Keep looking through daughters.
	G4LogicalVolume * LogicalVolume = PhysicalVolume->GetLogicalVolume();

	// Loop through the daughters of the volume
	G4int NumberDaughters = LogicalVolume->GetNoDaughters();
	for(G4int i=0; i!=NumberDaughters; ++i)
	  {
	    // Get the ith daughter volume
	    G4VPhysicalVolume * Daughter = LogicalVolume->GetDaughter(i);

	    // Recursively step into this volume
	    FindVolumes(Daughter, OpDetName, TransformSoFar, OpDetVolumes, OpDetTransformations);
	  }
      }
  }

}
