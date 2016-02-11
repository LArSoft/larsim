////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently mainly used to set optical properties
// for LAr and other optical components
//



#include "larsim/LArG4/MaterialPropertyLoader.h"
#include "lardata/Utilities/LArProperties.h"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace larg4 {

  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialProperty(std::string Material,
						   std::string Property, 
						   std::map<double, double> PropertyVector,
						   double Unit)
  {
    std::map<double,double> PropVectorWithUnit;
    for(std::map<double,double>::const_iterator it=PropertyVector.begin();
	it!=PropertyVector.end();
	it++)
      {
	PropVectorWithUnit[it->first*eV]=it->second*Unit;
      }
    fPropertyList[Material][Property]=PropVectorWithUnit;
    mf::LogInfo("MaterialPropertyLoader")<<"Added property " 
					 << Material<< "  " 
					 << Property;
  }
  
  //----------------------------------------------
  void MaterialPropertyLoader::SetMaterialConstProperty(std::string Material, 
							std::string Property, 
							double PropertyValue,
							double Unit)
  {
    fConstPropertyList[Material][Property]=PropertyValue*Unit;
    mf::LogInfo("MaterialPropertyLoader") << "Added const property " 
					  << Material << "  " 
					  << Property << " = " << PropertyValue;
  }
  
  //----------------------------------------------
  void MaterialPropertyLoader::SetBirksConstant(std::string Material, 
						double PropertyValue,
						double Unit)
  {
    fBirksConstants[Material]=PropertyValue*Unit;	
    mf::LogInfo("MaterialPropertyLoader") << "Set Birks constant " 
					  << Material;
  }

  //----------------------------------------------  
  void MaterialPropertyLoader::UpdateGeometry(G4LogicalVolumeStore * lvs)
  {
    std::map<std::string,G4MaterialPropertiesTable*> MaterialTables;
    std::map<std::string,bool> MaterialsSet;
    
    mf::LogInfo("MaterialPropertyLoader") << "UPDATING GEOMETRY";
    
    // Loop over each material with a property vector and create a new material table for it
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=fPropertyList.begin(); i!=fPropertyList.end(); i++){
      std::string Material=i->first;
      MaterialsSet[Material]=true;
      MaterialTables[Material]=new G4MaterialPropertiesTable;
    }
    
    // Loop over each material with a const property, 
    // if material table does not exist, create one
    for(std::map<std::string,std::map<std::string,double> >::const_iterator i=fConstPropertyList.begin(); i!=fConstPropertyList.end(); i++){
      std::string Material=i->first;
      if(!MaterialsSet[Material]){
	MaterialsSet[Material]=true;
	MaterialTables[Material]=new G4MaterialPropertiesTable;
      }
    }
    
    // For each property vector, convert to an array of g4doubles and 
    // feed to materials table Lots of firsts and seconds!  See annotation 
    // in MaterialPropertyLoader.h to follow what each element is
    
    for(std::map<std::string,std::map<std::string,std::map<double,double> > >::const_iterator i=fPropertyList.begin(); i!=fPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,std::map<double,double> >::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
	std::string Property=j->first;
	std::vector<G4double> g4MomentumVector;
	std::vector<G4double> g4PropertyVector;
	
	for(std::map<double,double>::const_iterator k=j->second.begin(); k!=j->second.end(); k++){
	  g4MomentumVector.push_back(k->first);
	  g4PropertyVector.push_back(k->second);
	}
	int NoOfElements=g4MomentumVector.size();
	MaterialTables[Material]->AddProperty(Property.c_str(),&g4MomentumVector[0], &g4PropertyVector[0],NoOfElements); 
	mf::LogInfo("MaterialPropertyLoader") << "Added property "
					      <<Property
					      <<" to material table " 
					      << Material;
      }
    }
    
    //Add each const property element
    for(std::map<std::string,std::map<std::string,double > >::const_iterator i = fConstPropertyList.begin(); i!=fConstPropertyList.end(); i++){
      std::string Material=i->first;
      for(std::map<std::string,double>::const_iterator j = i->second.begin(); j!=i->second.end(); j++){
	std::string Property=j->first;
	G4double PropertyValue=j->second;
	MaterialTables[Material]->AddConstProperty(Property.c_str(), PropertyValue); 
	mf::LogInfo("MaterialPropertyLoader") << "Added const property "
					      <<Property
					      <<" to material table " 
					      << Material;
      }
    }
    
    //Loop through geometry elements and apply relevant material table where materials match
    for ( G4LogicalVolumeStore::iterator i = lvs->begin(); i != lvs->end(); ++i ){
      G4LogicalVolume* volume = (*i);
      G4Material* TheMaterial = volume->GetMaterial();
      std::string Material = TheMaterial->GetName();
      for(std::map<std::string,G4MaterialPropertiesTable*>::const_iterator j=MaterialTables.begin(); j!=MaterialTables.end(); j++){
	if(Material==j->first){
	  TheMaterial->SetMaterialPropertiesTable(j->second);
	  //Birks Constant, for some reason, must be set separately
	  if(fBirksConstants[Material]!=0)
	    TheMaterial->GetIonisation()->SetBirksConstant(fBirksConstants[Material]);
	  volume->SetMaterial(TheMaterial);
	}
      }
    }
  }


  void MaterialPropertyLoader::SetReflectances(std::string /*Material*/, std::map<std::string,std::map<double, double> > Reflectances,  std::map<std::string,std::map<double, double> >  DiffuseFractions)
  {
    std::map<double, double> ReflectanceToStore;
    std::map<double, double> DiffuseToStore;
    
    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=Reflectances.begin();
	itMat!=Reflectances.end();
	++itMat)
      {
	std::string ReflectancePropName = std::string("REFLECTANCE_") + itMat->first;
	ReflectanceToStore.clear();
	for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
	    itEn!=itMat->second.end();
	    ++itEn)	  
	  {
	    ReflectanceToStore[itEn->first]=itEn->second;
	  }    
	SetMaterialProperty("LAr", ReflectancePropName, ReflectanceToStore,1);
      }

 for(std::map<std::string,std::map<double,double> >::const_iterator itMat=DiffuseFractions.begin();
	itMat!=DiffuseFractions.end();
	++itMat)
      {
	std::string DiffusePropName = std::string("DIFFUSE_REFLECTANCE_FRACTION_") + itMat->first;
	DiffuseToStore.clear();
	for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
	    itEn!=itMat->second.end();
	    ++itEn)	  
	  {
	    DiffuseToStore[itEn->first]=itEn->second;
	  }    
	SetMaterialProperty("LAr", DiffusePropName, DiffuseToStore,1);
      }
    
  } 
  

  void MaterialPropertyLoader::GetPropertiesFromServices()
  {
    art::ServiceHandle<util::LArProperties>   LarProp;
    
    // wavelength dependent quantities

    SetMaterialProperty( "LAr", "FASTCOMPONENT", LarProp->FastScintSpectrum(), 1  );
    SetMaterialProperty( "LAr", "SLOWCOMPONENT", LarProp->SlowScintSpectrum(), 1  );
    SetMaterialProperty( "LAr", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "LAr", "ABSLENGTH",     LarProp->AbsLengthSpectrum(), cm );
    SetMaterialProperty( "LAr", "RAYLEIGH",      LarProp->RayleighSpectrum(),  cm );


    // scalar properties

    SetMaterialConstProperty("LAr", "SCINTILLATIONYIELD",  LarProp->ScintYield(true),       1/MeV ); // true = scaled down by prescale in larproperties
    SetMaterialConstProperty("LAr", "RESOLUTIONSCALE",     LarProp->ScintResolutionScale(), 1);
    SetMaterialConstProperty("LAr", "FASTTIMECONSTANT",    LarProp->ScintFastTimeConst(),   ns);
    SetMaterialConstProperty("LAr", "SLOWTIMECONSTANT",    LarProp->ScintSlowTimeConst(),   ns);
    SetMaterialConstProperty("LAr", "YIELDRATIO",          LarProp->ScintYieldRatio(),      1);
    SetMaterialConstProperty("LAr", "ELECTRICFIELD",       LarProp->Efield(),               kilovolt/cm);

    SetBirksConstant("LAr",LarProp->ScintBirksConstant(), cm/MeV);
    
    SetReflectances("LAr", LarProp->SurfaceReflectances(), LarProp->SurfaceReflectanceDiffuseFractions());


    // If we are using scint by particle type, load these

    if(LarProp->ScintByParticleType())
      {
        // true = scaled down by prescale in larproperties
	SetMaterialConstProperty("LAr", "PROTONSCINTILLATIONYIELD",  LarProp->ProtonScintYield(true),    1./MeV );
	SetMaterialConstProperty("LAr", "PROTONYIELDRATIO",          LarProp->ProtonScintYieldRatio(),   1.);
	SetMaterialConstProperty("LAr", "MUONSCINTILLATIONYIELD",    LarProp->MuonScintYield(true),      1./MeV );
	SetMaterialConstProperty("LAr", "MUONYIELDRATIO",            LarProp->MuonScintYieldRatio(),     1.);
	SetMaterialConstProperty("LAr", "KAONSCINTILLATIONYIELD",    LarProp->KaonScintYield(true),      1./MeV );
	SetMaterialConstProperty("LAr", "KAONYIELDRATIO",            LarProp->KaonScintYieldRatio(),     1.);
	SetMaterialConstProperty("LAr", "PIONSCINTILLATIONYIELD",    LarProp->PionScintYield(true),      1./MeV );
	SetMaterialConstProperty("LAr", "PIONYIELDRATIO",            LarProp->PionScintYieldRatio(),     1.);
	SetMaterialConstProperty("LAr", "ELECTRONSCINTILLATIONYIELD",LarProp->ElectronScintYield(true),  1./MeV );
	SetMaterialConstProperty("LAr", "ELECTRONYIELDRATIO",        LarProp->ElectronScintYieldRatio(), 1.);
      	SetMaterialConstProperty("LAr", "ALPHASCINTILLATIONYIELD",   LarProp->AlphaScintYield(true),     1./MeV );
	SetMaterialConstProperty("LAr", "ALPHAYIELDRATIO",           LarProp->AlphaScintYieldRatio(),    1.);
      }
  }

}
