////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently used to supply optical properties to LAr
// and other optical elements of the detector.
//

#ifndef LArG4_MaterialPropertyLoader_h
#define LArG4_MaterialPropertyLoader_h

#include "Geant4/G4LogicalVolumeStore.hh"
#include <map>

namespace larg4 {
  
  class MaterialPropertyLoader
  {
  public:
    
    MaterialPropertyLoader() {}
    ~MaterialPropertyLoader() {}
    
    
    
  public:
    
    //Accessors
    std::map<double,double> GetMaterialProperty(std::string Material,std::string Property) 
      {return fPropertyList[Material][Property];}
    
    double GetMaterialConstProperty(std::string Material, std::string Property)
      {return fConstPropertyList[Material][Property];}
    
    std::map<std::string,double> GetMaterialConstProperties(std::string Material)
      {return fConstPropertyList[Material];}
    
    std::map<std::string,std::map<double,double> >  GetMaterialProperties(std::string Material)
      {return fPropertyList[Material];}
    
    
    // Methods to set material properties
    void SetMaterialProperty(       std::string Material, std::string Property, std::map<double,double> Values, double Unit);
    void SetMaterialConstProperty(  std::string Material, std::string Property, double Value,                   double Unit);

    // Method to set LArG4 Birks constant
    void SetBirksConstant( std::string, double, double );
    
    void SetReflectances( std::string, std::map<std::string, std::map<double,double> >, std::map<std::string, std::map<double, double> >);
    
    // Set material properties supplied from services
    void GetPropertiesFromServices();

    // Geometry updater.  Generate Geant4 material property tables and attach to detector geometry
    void UpdateGeometry( G4LogicalVolumeStore* );
    


  private:

    //         materials                properties            values
    std::map < std::string , std::map < std::string,double> > fConstPropertyList;

    //         materials                properties               energies  values
    std::map < std::string , std::map < std::string , std::map < double ,  double > > > fPropertyList;

    std::map<std::string, double> fBirksConstants;
    
    

    
  };  
  
}
#endif
