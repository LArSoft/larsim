////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsBuiltIns.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////


// The parameters of GEANT4 built in physics modules must be supplied to their constructors.
// This code initiates CustomPhysicsFactories to pass the default QGSP_BERT parameters
// and add built in physics to the CustomPhysicsFactory.  Parameters are passed exactly
// as they were in QGSP_BERT (hence apparent inconsistencies are not the fault of the current
// author!
//
// For new modules, parameters for new physics should be read from xml files within the modules
// themselves rather than passing parameters directly to the constructor. See OpticalPhysics as
// an example.
//
// Therefore no new modules should be initialised as those below (this is a somewhat
// hacky patch to reconcile the Geant4 way of declaring physics with the physics table
// scheme)


#include "Geant4/G4DecayPhysics.hh"
#include "Geant4/G4EmExtraPhysics.hh"
#include "Geant4/G4IonPhysics.hh"
#include "Geant4/G4StoppingPhysics.hh"
#include "Geant4/G4HadronElasticPhysics.hh"
#include "Geant4/G4NeutronTrackingCut.hh"
#include "Geant4/G4HadronPhysicsQGSP_BERT.hh"
#include "Geant4/G4EmStandardPhysics.hh"
#include "Geant4/G4EmLivermorePhysics.hh"
#include "Geant4/G4ChargeExchangePhysics.hh"

#include "larsim/LArG4/CustomPhysicsFactory.hh"


namespace larg4 {

  class EmPhysicsFactory : public CustomPhysicsFactory<G4EmStandardPhysics>
  {
  public:
    G4VPhysicsConstructor * Build()  {return new G4EmStandardPhysics(0);}
    EmPhysicsFactory() : CustomPhysicsFactory<G4EmStandardPhysics>("Em") {}
    virtual ~EmPhysicsFactory() {};
  };

  class SynchrotronAndGN : public CustomPhysicsFactory<G4EmExtraPhysics>
  {
  public:
    G4VPhysicsConstructor * Build() {return new G4EmExtraPhysics("extra EM");}
    SynchrotronAndGN() : CustomPhysicsFactory<G4EmExtraPhysics>("SynchrotronAndGN") {}
    virtual ~SynchrotronAndGN() {}
  };

  class ChargeExchange : public CustomPhysicsFactory<G4ChargeExchangePhysics>
  {
  public:
    G4VPhysicsConstructor * Build() {return new G4ChargeExchangePhysics(2);}
    ChargeExchange() : CustomPhysicsFactory<G4ChargeExchangePhysics>("ChargeExchange") {}
    virtual ~ChargeExchange() {}
  };


  class DecayPhysicsFactory : public CustomPhysicsFactory<G4DecayPhysics>
  {
  public:
    G4VPhysicsConstructor * Build() {return new G4DecayPhysics("Decay",0);}
    DecayPhysicsFactory() : CustomPhysicsFactory<G4DecayPhysics>("Decay") {}
    virtual ~DecayPhysicsFactory() {}
  };

  class HadronElasticPhysicsFactory : public CustomPhysicsFactory<G4HadronElasticPhysics>
  {
  public:
    G4VPhysicsConstructor * Build()  {return new G4HadronElasticPhysics(0);}
    HadronElasticPhysicsFactory() : CustomPhysicsFactory<G4HadronElasticPhysics>("HadronElastic") {}
   virtual ~HadronElasticPhysicsFactory() {}
  };

  class HadronPhysicsFactory : public CustomPhysicsFactory<G4HadronPhysicsQGSP_BERT>
  {
  public:
    G4bool quasiElastic;
    G4VPhysicsConstructor * Build() {return new G4HadronPhysicsQGSP_BERT("hadron",quasiElastic=true);}
    HadronPhysicsFactory() : CustomPhysicsFactory<G4HadronPhysicsQGSP_BERT>("Hadron") {}
   virtual ~HadronPhysicsFactory() {}
  };

  class StoppingPhysicsFactory : public CustomPhysicsFactory<G4StoppingPhysics>
  {
  public:
    G4VPhysicsConstructor * Build()  {return new G4StoppingPhysics("stopping");}
    StoppingPhysicsFactory() : CustomPhysicsFactory<G4StoppingPhysics>("Stopping") {}
    virtual ~StoppingPhysicsFactory() {};
  };

  class IonPhysicsFactory : public CustomPhysicsFactory<G4IonPhysics>
  {
  public:
    G4VPhysicsConstructor * Build()  {return new G4IonPhysics("ion");}
    IonPhysicsFactory() : CustomPhysicsFactory<G4IonPhysics>("Ion") {}
    virtual ~IonPhysicsFactory() {};
  };

  class NeutronTrackingCutFactory : public CustomPhysicsFactory<G4NeutronTrackingCut>
  {
  public:
    G4VPhysicsConstructor * Build()  {return new G4NeutronTrackingCut("Neutron tracking cut",0);}
    NeutronTrackingCutFactory() : CustomPhysicsFactory<G4NeutronTrackingCut>("NeutronTrackingCut") {}
    virtual ~NeutronTrackingCutFactory() {};
  };

    class LowEnergyEmFactory : public CustomPhysicsFactory<G4EmLivermorePhysics>
    {
    public:
      G4VPhysicsConstructor * Build() {return new G4EmLivermorePhysics();}
      LowEnergyEmFactory() : CustomPhysicsFactory<G4EmLivermorePhysics>("LowEnergyEm"){}
      virtual ~LowEnergyEmFactory() {}
    };


}





// Sept 2009 - Ben Jones, MIT
