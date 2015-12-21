////////////////////////////////////////////////////////////////////////
/// \file CustomPhysicsFactory.h
//
/// \version $Id: PrimaryParticleInformation.cxx,v 1.3 2009/10/05 23:21:51 t962cvs Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
#ifndef CUSTOMPHYSICSFACTORY_hh
#define CUSTOMPHYSICSFACTORY_hh

#include "Geant4/G4VPhysicsConstructor.hh"

#ifndef CUSTOMPHYSICSTABLE_hh
namespace larg4 { class CustomPhysicsTable; }
#endif

namespace larg4 {
  class CustomPhysicsFactoryBase
  {
  public:
    CustomPhysicsFactoryBase() {}
    virtual ~CustomPhysicsFactoryBase() {}
    
    virtual bool Registered() = 0;
    virtual std::string GetName() = 0;
    virtual G4VPhysicsConstructor * Build() = 0;
  };
  
  template<class T>
  class CustomPhysicsFactory : public CustomPhysicsFactoryBase
  {
  public:
    CustomPhysicsFactory();
    CustomPhysicsFactory(std::string);
    virtual ~CustomPhysicsFactory() {};
    bool Registered () 
    {return registered;}
    std::string GetName() 
    { return myName;}
    virtual G4VPhysicsConstructor * Build();
    
    
  private:
    std::string myName;
    bool registered;
    bool verbose;
  };
}  
#include "larsim/LArG4/CustomPhysicsTable.hh"
 
namespace larg4 {
 
  template<class T> G4VPhysicsConstructor * CustomPhysicsFactory<T>::Build()
  {
    return new T();
  }
  
  template<class T> CustomPhysicsFactory<T>::CustomPhysicsFactory(std::string Name)
  {
    
    // For debugging. 
    verbose=true;
    
    if(Name!="")
      myName=Name;
    else
      std::cerr<<"CustomPhysicsFactory Error : Physics registered with no name!"<<std::endl;
    
    // register self in physics table - note, factory is actually registered
    // in static TheCustomPhysicsTable, not the instance created below
    // which just acts to pass information along
    new CustomPhysicsTable(this);
    registered=true;
  }
  
  
  
  template<class T> CustomPhysicsFactory<T>::CustomPhysicsFactory()
  {
    registered=false;
  }
  
}
  
#endif
  

// Sept 2009 - Ben Jones, MIT
