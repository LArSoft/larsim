////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadoutGeometry.cxx
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
///
/// \version $Id: LArVoxelReadoutGeometry.cxx,v 1.3 2009/03/31 17:58:39 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

#include "G4Base/DetectorConstruction.h"

// C/C++ libraries
#include <vector>
#include <cmath>
#include <map>
#include <memory> // std::unique_ptr()
#include <sstream> // std::ostringstream
#include <iostream> // std::endl

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larsim/LArG4/LArVoxelReadoutGeometry.h"
#include "larsim/LArG4/LArVoxelReadout.h"
#include "larsim/Simulation/LArVoxelCalculator.h"

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


constexpr bool DisableVoxelCaching = false;

#include <typeinfo>

#ifdef __GNUG__
#include <cxxabi.h>

template <typename T>
std::string demangle_cxx_symbol(const T& obj) {
	int status = -4; // some arbitrary value to eliminate the compiler warning
	
	std::string name = typeid(obj).name();
	
	// __cxa_demangle() allocates with malloc(), so we ask unique_ptr to
	// deallocate by free()
	std::unique_ptr<char, void(*)(void*)> res
		{ abi::__cxa_demangle(name.c_str(), NULL, NULL, &status), std::free };
	
	return (status == 0)? res.get(): name;
} // demangle_cxx_symbol()

#else

template <typename T>
std::string demangle_cxx_symbol(const T& obj) { return typeid(obj).name(); }

#endif


template <class STREAM>
int DumpPhysicalVolume
	(STREAM& out, const G4VPhysicalVolume& PV, std::string indentstr = "")
{
	
	const G4ThreeVector& pos = PV.GetTranslation();
	const G4LogicalVolume* LV = PV.GetLogicalVolume();
	
	int count = 1;
	out << indentstr << PV.GetName() << " [" << demangle_cxx_symbol(PV) << "]"
		<< " at (" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")";
	if (LV) {
		out << ", a " << LV->GetName();
		
		const G4VSolid* Solid = LV->GetSolid();
		if (Solid) {
			out << " shaped as " << Solid->GetName();
			const G4Box* pBox;
			
			try { pBox = dynamic_cast<const G4Box*>(Solid); }
			catch (std::bad_cast&) { pBox = nullptr; }
			
			if (pBox) {
				out << ", a (" << (2.*pBox->GetXHalfLength())
					<< " x " << (2.*pBox->GetYHalfLength())
					<< " x " << (2.*pBox->GetZHalfLength()) << ") cm box";
			} // if box
			else {
				out << ", a " << demangle_cxx_symbol(*Solid);
			}
		}
		else {
			out << " with no shape (!?!)";
		}
		
		G4int nDaughters = LV->GetNoDaughters();
		if (nDaughters > 0) {
			out << " with " << nDaughters << " subvolumes:\n";
			for (G4int i = 0; i < nDaughters; ++i) {
				count += DumpPhysicalVolume(out, *LV->GetDaughter(i), indentstr + "  ");
			}
		}
		else out << '\n';
	}
	else out << " with no logical volume (!?)\n";
	
	return count;
} // DumpPhysicalVolume()


namespace larg4 {

  // Constructor and destructor.
  LArVoxelReadoutGeometry::LArVoxelReadoutGeometry(const G4String name)
    : G4VUserParallelWorld(name)
  {
    larg4::IonizationAndScintillation *ios = larg4::IonizationAndScintillation::Instance();
    std::unique_ptr<G4UserLimits> fStepLimit(new G4UserLimits(ios->StepSizeLimit()));
  }

  ////////////////////////////////////////////////////////////////////
  LArVoxelReadoutGeometry::~LArVoxelReadoutGeometry() 
  {
  }

  ////////////////////////////////////////////////////////////////////
  void LArVoxelReadoutGeometry::Construct()

  {
    // With a "parallel geometry", Geant4 has already created a clone
    // of the world physical and logical volumes.  We want to place
    // the LAr TPC, and only the LAr TPC, within this cloned world.

    // Get the parallel world physical volume.
    G4VPhysicalVolume* parallelPhysical = GetWorld();

    // Now we want to place a parallel LAr TPC volume within this
    // parallel world volume.  We only want to duplicate the LAr TPC
    // volume; any other volumes in the "official" geometry are going
    // to be ignored.  Our parallel world will consist only of the
    // world volume and the LAr TPC volume.

    G4Transform3D      worldTransform;

    // first get the volDetEnclosure
    std::string        daughterName("volDetEnclosure");
    G4Transform3D      detEnclosureTransform;
    G4VPhysicalVolume* detEnclosureVolume = this->FindNestedVolume(g4b::DetectorConstruction::GetWorld(),
                                                                   worldTransform,
                                                                   detEnclosureTransform,
                                                                   daughterName,
                                                                   0);



    // we keep track of all the voxel boxes with different geometries we create
    struct VoxelSpecs_t {
      double w, h, d; ///< sizes of the volume
      unsigned int nw, nh, nd; ///< divisions in each volume
      
      /// have some sorting...
      bool operator< (const VoxelSpecs_t& vs) const
        {
          if (w < vs.w) return true;
          if (w > vs.w) return false;
          if (h < vs.h) return true;
          if (h > vs.h) return false;
          if (d < vs.d) return true;
          if (d > vs.d) return false;
          if (nw < vs.nw) return true;
          if (nw > vs.nw) return false;
          if (nh < vs.nh) return true;
          if (nh > vs.nh) return false;
          if (nd < vs.nd) return true;
          return false;
        } // operator<
    }; // VoxelSpecs_t
    
    // TODO we don't need to cache the cell any more: a simple map will suffice
    struct VoxelVolumes_t {
      G4LogicalVolume* pBox;
      G4LogicalVolume* pCell;
      
      VoxelVolumes_t
        (G4LogicalVolume* box = nullptr, G4LogicalVolume* cell = nullptr):
        pBox(box), pCell(cell) {}
    }; // VoxelVolumes_t
    
    
    // Define the sensitive detector for the voxel readout.  This class
    // routines will be called every time a particle deposits energy in
    // a voxel that overlaps the LAr TPC.
    LArVoxelReadout* larVoxelReadout = new LArVoxelReadout("LArVoxelSD");
    if ((fGeo->Ncryostats() == 1) && (fGeo->Cryostat(0).NTPC() == 1))
      larVoxelReadout->SetSingleTPC(0, 0); // just one TPC in the detector...

    // Tell Geant4's sensitive-detector manager that the voxel SD
    // class exists.
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    sdManager->AddNewDetector(larVoxelReadout);

    // hope hashing doubles is reliable...
    typedef std::map<VoxelSpecs_t, VoxelVolumes_t> VoxelCache_t;
    VoxelCache_t VoxelCache;
    
    for(unsigned int c = 0; c < fGeo->Ncryostats(); ++c){

      // next get the cryostat
      daughterName = "volCryostat";
      G4Transform3D       cryostatTransform;
      G4VPhysicalVolume*  cryostatVolume = this->FindNestedVolume(detEnclosureVolume,
                                                                  detEnclosureTransform,
                                                                  cryostatTransform,
                                                                  daughterName,
                                                                  c);

      for(unsigned int t = 0; t < fGeo->Cryostat(c).NTPC(); ++t){

        // now for the TPC
        daughterName = "volTPC";    
        G4Transform3D       tpcTransform;
        G4VPhysicalVolume*  tpcVolume = this->FindNestedVolume(cryostatVolume, 
                                                               cryostatTransform,
                                                               tpcTransform,
                                                               daughterName,
                                                               t);

        daughterName = "volTPCActive";    
        G4Transform3D       transform;
        G4VPhysicalVolume*  larTPCPhysical = this->FindNestedVolume(tpcVolume, 
                                                                    tpcTransform,
                                                                    transform,
                                                                    daughterName,
                                                                    t);

        // Get the LAr TPC volume, and its shape.
        G4LogicalVolume* larTPCLogical = larTPCPhysical->GetLogicalVolume();
        larTPCLogical->SetUserLimits(fStepLimit.get());

        G4VSolid* larTPCShape = larTPCLogical->GetSolid();
        
        // We're not going to exactly duplicate the LAr TPC in our
        // parallel world.  We're going to construct a box of voxels.
        // What should the size and position of that box be?

        // To get our first hints, we need the overall dimensions of a
        // "bounding box" that contains the shape.  For now, we'll allow
        // two possible shapes: a box (by the far the most likely) and a
        // cylinder (for bizarre future detectors that I know nothing
        // about).

        G4double larTPCHalfXLength = 0;
        G4double larTPCHalfYLength = 0;
        G4double larTPCHalfZLength = 0;
        G4Box* tpcBox = dynamic_cast< G4Box* >( larTPCShape );
        if ( tpcBox != 0 ){
          larTPCHalfXLength = tpcBox->GetXHalfLength();
          larTPCHalfYLength = tpcBox->GetYHalfLength();
          larTPCHalfZLength = tpcBox->GetZHalfLength();
        }
        else{
          // It's not a box.  Try a cylinder.
          G4Tubs* tube = dynamic_cast< G4Tubs* >( larTPCShape );
          if ( tube != 0 ){
            larTPCHalfXLength = tube->GetOuterRadius();
            larTPCHalfYLength = tube->GetOuterRadius();
            larTPCHalfZLength = tube->GetZHalfLength();
          }
          else{
            throw cet::exception("LArVoxelReadoutGeometry") << "Unknown shape in readout geometry"
                                                            << "The LAr TPC volume is not a box or a tube. "
                                                            << "This routine can't convert any other shapes.\n";
          }
        }

        LOG_DEBUG("LArVoxelReadoutGeometry") << ": larTPCHalfXLength=" << larTPCHalfXLength
                                             << ": larTPCHalfYLength=" << larTPCHalfYLength
                                             << ": larTPCHalfZLength=" << larTPCHalfZLength;

        // Get some constants from the LAr voxel information object.
        // Remember, ROOT uses cm.
        art::ServiceHandle<sim::LArVoxelCalculator> lvc;
        G4double voxelSizeX   = lvc->VoxelSizeX() * CLHEP::cm;
        G4double voxelSizeY   = lvc->VoxelSizeY() * CLHEP::cm;
        G4double voxelSizeZ   = lvc->VoxelSizeZ() * CLHEP::cm;
        G4double voxelOffsetX = lvc->VoxelOffsetX() * CLHEP::cm;
        G4double voxelOffsetY = lvc->VoxelOffsetY() * CLHEP::cm;
        G4double voxelOffsetZ = lvc->VoxelOffsetZ() * CLHEP::cm;

        LOG_DEBUG("LArVoxelReadoutGeometry") << ": voxelSizeX=" << voxelSizeX
                                             << ", voxelSizeY=" << voxelSizeY
                                             << ", voxelSizeZ=" << voxelSizeZ;

        LOG_DEBUG("LArVoxelReadoutGeometry") << ": voxelOffsetX=" << voxelOffsetX
                                             << ", voxelOffsetY=" << voxelOffsetY
                                             << ", voxelOffsetZ=" << voxelOffsetZ;
        
        // We want our voxelization region to be an integer multiple of
        // the voxel sizes in all directions; if we didn't do this, we
        // might get into trouble when we start playing with replicas.
        // Compute the the dimensions of our voxelization to be about the
        // size of the LAr TPC region, adjusted to be an integer number of
        // voxels in all directions.

        G4double numberXvoxels = 2.*larTPCHalfXLength / voxelSizeX;
        G4double numberYvoxels = 2.*larTPCHalfYLength / voxelSizeY;
        G4double numberZvoxels = 2.*larTPCHalfZLength / voxelSizeZ;
        numberXvoxels = trunc(numberXvoxels) + 1.;
        numberYvoxels = trunc(numberYvoxels) + 1.;
        numberZvoxels = trunc(numberZvoxels) + 1.;
        
        LOG_DEBUG("LArVoxelReadoutGeometry") << "Active volume of cryo #" << c << " TPC #" << t
          << " will be split in " << numberXvoxels << " x " << numberYvoxels << " x " << numberZvoxels
          << " = " << (numberXvoxels * numberYvoxels * numberZvoxels)
          << " voxels of size " << voxelSizeX << " x " << voxelSizeY << " x " << voxelSizeZ << " cm"
          ;

        VoxelSpecs_t VoxelSpecs{
          /* w */  2.*larTPCHalfXLength,
          /* h */  2.*larTPCHalfYLength,
          /* d */  2.*larTPCHalfZLength,
          /* nw */ (unsigned int) numberXvoxels,
          /* nh */ (unsigned int) numberYvoxels,
          /* nd */ (unsigned int) numberZvoxels
        };
        
        VoxelCache_t::iterator iVoxelVol = DisableVoxelCaching?
          VoxelCache.end(): VoxelCache.find(VoxelSpecs);
        if (iVoxelVol == VoxelCache.end()) {
          LOG_DEBUG("LArVoxelReadoutGeometry") << "Creating a new voxel volume "
            << VoxelSpecs.w << " x " << VoxelSpecs.h << " x " << VoxelSpecs.d
            << " cm, "
            << VoxelSpecs.nw << " x " << VoxelSpecs.nh << " x " << VoxelSpecs.nd
            << " voxels";
          
          G4double voxelBoxHalfX = numberXvoxels * voxelSizeX / 2.;
          G4double voxelBoxHalfY = numberYvoxels * voxelSizeY / 2.;
          G4double voxelBoxHalfZ = numberZvoxels * voxelSizeZ / 2.;

          LOG_DEBUG("LArVoxelReadoutGeometry") << ": voxelBoxHalfX=" << voxelBoxHalfX
                                               << ", voxelBoxHalfY=" << voxelBoxHalfY
                                               << ", voxelBoxHalfZ=" << voxelBoxHalfZ;
          
          // Now we have a box that will include an integer number of voxels
          // in each direction.  Note that the material is irrelevant for a
          // "parallel world."
          G4Box* voxelBox = new G4Box("VoxelBox",voxelBoxHalfX,voxelBoxHalfY,voxelBoxHalfZ);
          G4LogicalVolume* voxelBoxLogical = new G4LogicalVolume(voxelBox,
                                                                 0,
                                                                 "VoxelizationLogicalVolume" );
          
          // If we generate an event display within Geant4, we won't want to
          // see this box.
          G4VisAttributes* invisible = new G4VisAttributes();
          invisible->SetVisibility(false);
          voxelBoxLogical->SetVisAttributes(invisible);
          
          //LOG_DEBUG("LArVoxelReadoutGeometry") << ": transform = \n";
          //for ( G4int i = 0; i < 3; ++i ){
          //  for ( G4int j = 0; j < 4; ++j ){ 
          //    LOG_DEBUG("LArVoxelReadoutGeometry") << transform[i][j] << " ";
          //  }
          //  LOG_DEBUG("LArVoxelReadoutGeometry") << "\n";
          //}
          
          // Now we've fill our "box of voxels" with the voxels themselves.
          // We'll do this by sub-dividing the volume in x, then y, then z.
          
          // Create an "x-slice".
          G4Box* xSlice = new G4Box("xSlice",voxelSizeX/2.,voxelBoxHalfY,voxelBoxHalfZ);
          G4LogicalVolume* xSliceLogical = new G4LogicalVolume( xSlice, 0, "xLArVoxelSlice" );
          xSliceLogical->SetVisAttributes(invisible);
          
          // Use replication to slice up the "box of voxels" along the x-axis.
          new G4PVReplica( "VoxelSlicesInX",
                           xSliceLogical,
                           voxelBoxLogical,
                           kXAxis,
                           G4int( numberXvoxels ),
                           voxelSizeX );
          
          // Now do the same thing, dividing that x-slice along the y-axis.
          G4Box* ySlice = new G4Box("ySlice",voxelSizeX/2.,voxelSizeY/2., voxelBoxHalfZ);
          G4LogicalVolume* ySliceLogical = new G4LogicalVolume( ySlice, 0, "yLArVoxelSlice" );
          ySliceLogical->SetVisAttributes(invisible);
          new G4PVReplica( "VoxelSlicesInY",
                           ySliceLogical,
                           xSliceLogical,
                           kYAxis,
                           G4int( numberYvoxels ),
                           voxelSizeY );
          
          // Now divide the y-slice along the z-axis, giving us our actual voxels.
          G4Box* zSlice = new G4Box("zSlice",voxelSizeX/2.,voxelSizeY/2., voxelSizeZ/2.);
          G4LogicalVolume* voxelLogical = new G4LogicalVolume( zSlice, 0, "LArVoxel" );
          voxelLogical->SetVisAttributes(invisible);
          new G4PVReplica( "LArVoxel",
                           voxelLogical,
                           ySliceLogical,
                           kZAxis,
                           G4int( numberZvoxels ),
                           voxelSizeZ );
          
          iVoxelVol = VoxelCache.emplace(VoxelSpecs,
            VoxelVolumes_t(voxelBoxLogical, voxelLogical)
            ).first; // iterator to the inserted element

          // Set the sensitive detector of this LAr TPC (logical) volume
          // to be the voxel readout.
          voxelLogical->SetSensitiveDetector(larVoxelReadout);
          
        } // if not cached yet
        G4LogicalVolume* voxelBoxLogical = iVoxelVol->second.pBox;
        
        // We have a "box of voxels" that's the right dimensions, but we
        // have to know exactly where to put it.  The user has the option
        // to offset the voxel co-ordinate system.  We want to place our
        // box so the edges of our voxels align with that co-ordinate
        // system.  In effect, we want to offset our "box of voxels" by
        // the user's offsets, modulo the size of the voxel in each
        // direction.
        
        G4double offsetInVoxelsX = voxelOffsetX / voxelSizeX;
        G4double offsetInVoxelsY = voxelOffsetY / voxelSizeY;
        G4double offsetInVoxelsZ = voxelOffsetZ / voxelSizeZ;
        G4double fractionOffsetX = offsetInVoxelsX - trunc(offsetInVoxelsX);
        G4double fractionOffsetY = offsetInVoxelsY - trunc(offsetInVoxelsY);
        G4double fractionOffsetZ = offsetInVoxelsZ - trunc(offsetInVoxelsZ);
        G4double offsetX         = fractionOffsetX * voxelSizeX;
        G4double offsetY         = fractionOffsetY * voxelSizeY;
        G4double offsetZ         = fractionOffsetZ * voxelSizeZ;
        
        // Now we know how much to offset the "box of voxels".  Include
        // that in the transformation of the co-ordinates from world
        // volume to LAr TPC volume.
        transform = G4Translate3D( offsetX, offsetY, offsetZ ) * transform;
        
        LOG_DEBUG("LArVoxelReadoutGeometry") << ": offsetX=" << offsetX
                                             << ", offsetY=" << offsetY
                                             << ", offsetZ=" << offsetZ;
        
        // suffix representing this TPC
        std::ostringstream nums;
        nums << "_Cryostat" << c << "_TPC" << t;
        
        // Place the box of voxels, with the accumulated transformations
        // computed above.
        new G4PVPlacementInTPC( transform,
          "VoxelizationPhysicalVolume" + nums.str(),
          voxelBoxLogical,
          parallelPhysical,
          false,                                             // Only one volume
          0,                                                 // Copy number
          false,                                             // No surface check
          { (unsigned short int) c, (unsigned short int) t } // TPC ID
          );
        
      } // end loop over tpcs
    } // end loop over cryostats

    if (mf::isDebugEnabled()) {
      LOG_DEBUG("LArVoxelReadoutGeometryDump") << "Dump of voxelized volume";
      {
        mf::LogDebug log("LArVoxelReadoutGeometryDump");
        DumpPhysicalVolume(log, *parallelPhysical);
      }
      LOG_DEBUG("LArVoxelReadoutGeometryDump")
        << "End of dump of voxelized volume";
    }
    return;
  }

  //---------------------------------------------------------------
  ////////////////////////////////////////////////////////////////////
  // We know the ordering of the volumes in the Geometry, 
  // see https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Geometry
  // Make use of that knowledge to efficiently get the desired volumes and 
  // their total transforms.  
  // the daughterTransform is the total transform to the world coordinate system
  G4VPhysicalVolume* LArVoxelReadoutGeometry::FindNestedVolume(G4VPhysicalVolume* mother,
                                                               G4Transform3D&     motherTransform,
                                                               G4Transform3D&     daughterTransform,
                                                               std::string&       daughterName,
                                                               unsigned int       expectedNum)
  {
    G4LogicalVolume* logicalVolume = mother->GetLogicalVolume();
    G4int numberDaughters = logicalVolume->GetNoDaughters();
    for ( G4int i = 0; i != numberDaughters; ++i ){
      G4VPhysicalVolume* d = logicalVolume->GetDaughter(i);

    //  LOG_DEBUG("LArVoxelReadoutGeometry") << d->GetName() << ":" << mother->GetName();

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


        LOG_DEBUG("LArVoxelReadoutGeometry") << "current " << daughterName 
                                             << " origin is at (" 
                                             << world.x() / CLHEP::cm << ","
                                             << world.y() / CLHEP::cm << ","
                                             << world.z() / CLHEP::cm << ")";

        // we don't bother with the cryostat number when calling Geometry::PositionToTPC
        // because we know we have already started off with the correct cryostat volume
        // G4 uses mm, we want cm
        double worldPos[3] = { world.x() / CLHEP::cm, world.y() / CLHEP::cm, world.z() / CLHEP::cm };
        unsigned int daughterNum = 0;
        unsigned int extra       = 0;
        if(daughterName.compare("volCryostat") == 0)
          fGeo->PositionToCryostat(worldPos, daughterNum);
        else if(daughterName.compare("volTPC") == 0)
          fGeo->PositionToTPC(worldPos, daughterNum, extra);
        else if(daughterName.compare("volTPCActive") == 0 || 
                daughterName.compare("volDetEnclosure") == 0){
          // for either of these volumes, we know there is only 1 in the mother volume
          LOG_DEBUG("LArVoxelReadoutGeometry") << "found the desired " << daughterName;
          return d;
        }

        // if we found the desired volume, stop looking
        if(daughterNum == expectedNum){
          LOG_DEBUG("LArVoxelReadoutGeometry") << "found the desired " << daughterName;
          return d;
        }
      }// end if the volume has the right name
    }// end loop over volumes

    throw cet::exception("LArVoxelReadoutGeometry") << "could not find the desired "
                                                    << daughterName
                                                    << " to make LArVoxelReadoutGeometry\n";
    
    return 0;
  }


} // namespace larg4
