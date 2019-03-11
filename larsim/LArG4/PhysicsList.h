////////////////////////////////////////////////////////////////////////
/// \file  PhysicsList.h
/// \brief Create the physics lists to be used by Geant4.
///
/// \author  seligman@nevis.columbia.edu
//  \modified by: drivera@fnal.gov 
//  \changes: removed QGSP_BERT.h include because it is not pertinent 
//            here nor where this header is included 
//      
////////////////////////////////////////////////////////////////////////
///
/// Without a physics list, Geant4 won't do anything.  G4 comes with a
/// number of pre-constructed lists, and for now I plan to use
/// "QGSP_BERT".  It has the following properties:
///
/// - Standard EM physics processes.
/// - Quark-gluon string model for high energies (> 20GeV)
/// - Low Energy Parameterized (LEP) for medium energies (10<E<20GeV)
/// - Gertini-style cascade for low energies (< 10GeV)
/// - LEP, HEP for all anti-baryons (LEP,HEP = low/high energy parameterized, from GHEISHA)
/// - Gamma-nuclear model added for E<3.5 GeV
/// (comments from "Guided Tour of Geant4 Physics List II",
/// talk given at JPL by Dennis Wright)
///
/// IMPORTANT: For now, I'm just copying this physics list from the
/// work I did for NuSOnG, which in turn I copied from ATLAS.  More
/// thought is needed for the physics list for MicroBooNE.  
///
/// If you decide to replace QGSP_BERT with another of G4's
/// pre-supplied physics lists, you can just do a global replace on
/// that name in this header; nothing else need change (except the
/// comments, of course!).
///
/// Things become more complex when you're using parallel geometries.
/// The physics processes have to be duplicated in the parallel world.
/// There's a G4 facility for doing this, but it's not part of the
/// pre-supplied physics lists.
///
/// Hence this class.  It takes the physics from QGSP_BERT, but it
/// extends one of the methods to include recognition of physics
/// processes in the parallel world.
/// -------------------------------------------------------------------------------------------
/// 03/04/19 - D. Rivera
/// The physics list is not chosen in this header. Instead this header file defines the Modular
/// Physics list object. The actual available physics lists are limited to what is defined in the
/// CustomPhysicsBuiltIns.hh which is where various other physics lists can be included. 
/// 

#ifndef LArG4_PhysicsList_h
#define LArG4_PhysicsList_h

#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4VModularPhysicsList.hh"
#include "Geant4/G4VPhysicsConstructor.hh"
#include "Geant4/G4String.hh"
#include "Geant4/globals.hh"
#include "larsim/LArG4/ConfigurablePhysicsList.hh"

namespace larg4 {

  // Under normal circumstances, there is no need to inherit a class
  // from G4VModularPhysicsList.  It's necessary here because we have
  // to modify one of its routines.
  class ModularPhysicsList : public G4VModularPhysicsList
  {
  public:
    // Constructor.
    explicit ModularPhysicsList();
    virtual ~ModularPhysicsList();

    // This is the one method we're overridding to include the code
    // for the parallel geometry.
    virtual void ConstructProcess();

    // Non-virtual methods in G4VModularPhysicsList.  Just call the
    // method with the same name in G4VModularPhysicsList.
    void RegisterPhysics(G4VPhysicsConstructor* g) 
    { 
      G4VModularPhysicsList::RegisterPhysics(g); 
    }
    const G4VPhysicsConstructor* GetPhysics(G4int index) const 
    { 
      return G4VModularPhysicsList::GetPhysics(index); 
    }
    const G4VPhysicsConstructor* GetPhysics(const G4String& name) const 
    { 
      return G4VModularPhysicsList::GetPhysics(name); 
    }

  };

  /// This typedef is what defines the name "larg4::PhysicsList" in
  /// any class that includes this header. 
  typedef TConfigurablePhysicsList<ModularPhysicsList> PhysicsList;

} // namespace larg4

#endif // LArG4_PhysicsList_h
