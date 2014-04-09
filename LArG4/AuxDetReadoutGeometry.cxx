////////////////////////////////////////////////////////////////////////
/// \file  AuxDetReadoutGeometry.cxx
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
/// \author miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "G4Base/DetectorConstruction.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "LArG4/AuxDetReadoutGeometry.h"
#include "LArG4/AuxDetReadout.h"

// G4 includes
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4PVReplica.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4VisAttributes.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Point3D.hh"
#include "Geant4/globals.hh"

#include <vector>
#include <cmath>

namespace larg4 {

  // Constructor and destructor.
  AuxDetReadoutGeometry::AuxDetReadoutGeometry(const G4String name)
    : G4VUserParallelWorld(name)
  {}

  ////////////////////////////////////////////////////////////////////
  AuxDetReadoutGeometry::~AuxDetReadoutGeometry() 
  {}

  ////////////////////////////////////////////////////////////////////
  void AuxDetReadoutGeometry::Construct()

  {
    // With a "parallel geometry", Geant4 has already created a clone
    // of the world physical and logical volumes.  We want to place
    // the Auxiliary Detector, and only the Auxiliary Detector, within
    // this cloned world.

    // Get the parallel world physical volume.
    //G4VPhysicalVolume* parallelPhysical = GetWorld();

    // Now we want to place a parallel AuxDet volume within this
    // parallel world volume.  We only want to duplicate the AuxDet
    // volume; any other volumes in the "official" geometry are going
    // to be ignored.  Our parallel world will consist only of the
    // world volume and the AuxDet volume.

    G4Transform3D      worldTransform;

    // first get the volDetEnclosure
    std::string daughterName("volDetEnclosure");
    G4Transform3D      detEnclosureTransform;
    G4VPhysicalVolume* detEnclosureVolume = this->FindNestedVolume(g4b::DetectorConstruction::GetWorld(),
								                                                   worldTransform,
								                                                   detEnclosureTransform,
								                                                   daughterName,
								                                                   0);
    for(unsigned int c = 0; c < fGeo->NAuxDets(); ++c){
      // next get the Auxiliary Detector
      std::string volName = "volAuxDet";
		
      G4Transform3D       auxDetTransform;
      G4VPhysicalVolume*  auxDetVolume = this->FindNestedVolume(detEnclosureVolume,
								                                                detEnclosureTransform,
								                                                auxDetTransform,
								                                                volName,
								                                                c);
      
        // Get the AuxDet volume, and its shape.
        G4LogicalVolume* auxDetLogical = auxDetVolume->GetLogicalVolume();
        //G4VSolid* auxDetShape = auxDetLogical->GetSolid();
        
        // Define the sensitive detector for the auxiliary detector. This
        // routine will be called every time a particle deposits energy
        // in the AuxDet.
        std::string name("AuxDetSD");
        char nums[32];
        sprintf(nums, "_AuxDet%u", c);
        name += nums;
        AuxDetReadout* auxDetReadout = new AuxDetReadout(name);
        
        // Tell Geant4's sensitive-detector manager that the AuxDet
        // class exists.
        G4SDManager* sdManager = G4SDManager::GetSDMpointer();
        sdManager->AddNewDetector(auxDetReadout);
        
        // Set the sensitive detector of the AuxDet to be the AuxDet readout.
        auxDetLogical->SetSensitiveDetector(auxDetReadout);
        
    } // end loop over Auxiliary Detectors
    return;
  }// end Construct

  //---------------------------------------------------------------
  ////////////////////////////////////////////////////////////////////
  // We know the ordering of the volumes in the Geometry, 
  // see https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Geometry
  // Make use of that knowledge to efficiently get the desired volumes and 
  // their total transforms.  
  // the daughterTransform is the total transform to the world coordinate system
  G4VPhysicalVolume* AuxDetReadoutGeometry::FindNestedVolume(G4VPhysicalVolume* mother,
							                                               G4Transform3D&     motherTransform,
							                                               G4Transform3D&     daughterTransform,
							                                               std::string&       daughterName,
							                                               unsigned int       expectedNum)
  {
    G4LogicalVolume* logicalVolume = mother->GetLogicalVolume();
    G4int numberDaughters = logicalVolume->GetNoDaughters();
    for ( G4int i = 0; i != numberDaughters; ++i ){
      G4VPhysicalVolume* d = logicalVolume->GetDaughter(i);

      LOG_DEBUG("AuxDetReadoutGeometry") << d->GetName() << ":" << mother->GetName();
	  std::string volAuxDetNumString("volAuxDet");
	  char n[32];
      sprintf(n, "%u", expectedNum);
	  volAuxDetNumString += n;
	  volAuxDetNumString += "_PV";

      if(d->GetName().contains(daughterName)){

	      // check that this cryostat is the requested one using fCryostat
	      G4ThreeVector translation = d->GetObjectTranslation();
	      G4RotationMatrix rotation = d->GetObjectRotationValue();
	      G4Transform3D transform(rotation, translation);
	      daughterTransform = motherTransform * transform;


	      // take the origin of the volume and transform it to 
	      // world coordinated
	      G4Point3D local(0., 0., 0.);
	      G4Point3D world = daughterTransform * local;


	      LOG_DEBUG("AuxDetReadoutGeometry") << "current daughter=" << daughterName 
					     << " origin is at (" 
					     << world.x() / cm << ","
					     << world.y() / cm << ","
					     << world.z() / cm << ")";

	      // we don't bother with the cryostat number when calling Geometry::PositionToTPC
	      // because we know we have already started off with the correct cryostat volume
	      // G4 uses mm, we want cm
	      double worldPos[3] = { world.x() / cm, world.y() / cm, world.z() / cm };

        if(d->GetName().compare(volAuxDetNumString)==0){
		      fGeo->PositionToAuxDet(worldPos, expectedNum);
		      LOG_DEBUG("AuxDetReadoutGeometry") << "found the desired " << daughterName;
		      return d;
	      }else if(d->GetName().compare("volDetEnclosure_PV") == 0){
	        // for either of these volumes, we know there is only 1 in the mother volume
	        LOG_DEBUG("AuxDetReadoutGeometry") << "found the desired " << daughterName;
	        return d;
	      }
      }// end if the volume contains the right name
    }// end loop over daughter volumes

    throw cet::exception("AuxDetReadoutGeometry") << "could not find the desired "
						    << daughterName
						    << " to make AuxDetReadoutGeometry\n";
    
    return 0;
  }// end FindNestedVolume


} // namespace larg4
