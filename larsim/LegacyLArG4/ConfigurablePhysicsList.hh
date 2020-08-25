////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsList.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Configurable Physics List Class
//
// Ben Jones, MIT, 24/06/09
//
// Based on the QGSP_BERT physics list supplied with geant 4, but with
// options to switch on and off different physics processed from the config
// file.
//
// To include a new physics process, create a new physics builder registering
// the relevant particles and interactions.  Then add the builder to the
// GetPhysicsBuilders and GetDefaultSettings functions in this class,
// using the name by which the builder will be referenced in the config
// file.  Physics builders to be enabled are specified in the LArG4 config.
//
//

#ifndef TConfigurablePhysicsList_h
#define TConfigurablePhysicsList_h 1

#include "Geant4/CompileTimeConstraints.hh"
#include "Geant4/G4VModularPhysicsList.hh"
#include "Geant4/globals.hh"

namespace larg4 {

  template <class T>
  class TConfigurablePhysicsList : public T {
  public:
    TConfigurablePhysicsList(G4int ver = 1);
    virtual ~TConfigurablePhysicsList();
    bool UseCustomPhysics;
    std::vector<std::string> EnabledPhysics;
    std::vector<std::string> GetDefaultSettings();
    virtual void SetCuts();
  };

}

#include "larsim/LegacyLArG4/ConfigurablePhysicsList.icc"

namespace larg4 {
  using ConfigurablePhysicsList = TConfigurablePhysicsList<G4VModularPhysicsList>;
}

#endif
