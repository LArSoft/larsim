////////////////////////////////////////////////////////////////////////
/// \file MaterialPropertyLoader.cxx
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently mainly used to set optical properties
// for LAr and other optical components
//

// TODO convert tabs into spaces

// TODO verify the inclusion list
#include "larsim/LArG4/MaterialPropertyLoader.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4LogicalSkinSurface.hh"
#include "Geant4/G4OpticalSurface.hh"

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
        PropVectorWithUnit[it->first*CLHEP::eV]=it->second*Unit;
      }
    fPropertyList[Material][Property]=PropVectorWithUnit;
    // replace with MF_LOGDEBUG()
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
    // replace with MF_LOGDEBUG()
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
    // replace with MF_LOGDEBUG()
    mf::LogInfo("MaterialPropertyLoader") << "Set Birks constant "
                                          << Material;
  }

  //----------------------------------------------
  void MaterialPropertyLoader::UpdateGeometry(G4LogicalVolumeStore * lvs)
  {
    std::map<std::string,G4MaterialPropertiesTable*> MaterialTables;
    std::map<std::string,bool> MaterialsSet;

    // TODO replace console output with messagefacility output
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
        // replace with mf::LogVerbatim()
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
        // replace with mf::LogVerbatim()
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

      //
      // create reflective surfaces corresponding to the volumes made of some
      // selected materials
      //

      //--------------------------> FIXME <-----------------(parameters from fcl files(?))
      G4MaterialPropertyVector* PropertyPointer = 0;
      if(MaterialTables[Material])
        PropertyPointer = MaterialTables[Material]->GetProperty("REFLECTIVITY");

      if(Material=="Copper"){
        std::cout<< "copper foil surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining Copper optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfc = new G4OpticalSurface("Surface copper",glisur,ground,dielectric_metal);
          refl_opsurfc->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfc->SetPolish(0.2);
          new G4LogicalSkinSurface("refl_surfacec",volume, refl_opsurfc);
        }
        else
          std::cout<< "Warning: Copper surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }

      if(Material=="G10"){
        std::cout<< "G10 surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining G10 optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfg = new G4OpticalSurface("g10 Surface",glisur,ground,dielectric_metal);
          refl_opsurfg->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfg->SetPolish(0.1);
          new G4LogicalSkinSurface("refl_surfaceg",volume, refl_opsurfg);
        }
        else
          std::cout<< "Warning: G10 surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }

      if(Material=="vm2000"){
        std::cout<< "vm2000 surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining vm2000 optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurf = new G4OpticalSurface("Reflector Surface",unified,groundfrontpainted,dielectric_dielectric);
          refl_opsurf->SetMaterialPropertiesTable(MaterialTables[Material]);
          G4double sigma_alpha = 0.8;
          refl_opsurf->SetSigmaAlpha(sigma_alpha);
          new G4LogicalSkinSurface("refl_surface",volume, refl_opsurf);
        }
        else
          std::cout<< "Warning: vm2000 surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      if(Material=="ALUMINUM_Al"){
        std::cout<< "ALUMINUM_Al surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining ALUMINUM_Al optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface Aluminum",glisur,ground,dielectric_metal);
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.5);
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: ALUMINUM_Al surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      if(Material=="STEEL_STAINLESS_Fe7Cr2Ni"){
        std::cout<< "STEEL_STAINLESS_Fe7Cr2Ni surface set "<<volume->GetName()<<std::endl;
        if(PropertyPointer) {
          std::cout<< "defining STEEL_STAINLESS_Fe7Cr2Ni optical boundary "<<std::endl;
          G4OpticalSurface* refl_opsurfs = new G4OpticalSurface("Surface Steel",glisur,ground,dielectric_metal);
          refl_opsurfs->SetMaterialPropertiesTable(MaterialTables[Material]);
          refl_opsurfs->SetPolish(0.5);
          new G4LogicalSkinSurface("refl_surfaces",volume, refl_opsurfs);
        }
        else
          std::cout<< "Warning: STEEL_STAINLESS_Fe7Cr2Ni surface in the geometry without REFLECTIVITY assigned"<<std::endl;
      }
      //-----------------------------------------------------------------------------

      //
      // apply the remaining material properties
      //
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


  void MaterialPropertyLoader::SetReflectances(std::map<std::string,std::map<double, double> > Reflectances)
  {
    std::map<double, double> ReflectanceToStore;

    for(std::map<std::string,std::map<double,double> >::const_iterator itMat=Reflectances.begin();
        itMat!=Reflectances.end();
        ++itMat)
      {
        ReflectanceToStore.clear();
        for(std::map<double,double>::const_iterator itEn=itMat->second.begin();
            itEn!=itMat->second.end();
            ++itEn)
          {
            ReflectanceToStore[itEn->first]=itEn->second;
          }
        SetMaterialProperty(itMat->first, "REFLECTIVITY", ReflectanceToStore,1);
      }
  }


  void MaterialPropertyLoader::GetPropertiesFromServices()
  {
    const detinfo::LArProperties* LarProp = lar::providerFrom<detinfo::LArPropertiesService>();
    const detinfo::DetectorProperties* DetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // wavelength dependent quantities

    SetMaterialProperty( "LAr", "FASTCOMPONENT", LarProp->FastScintSpectrum(), 1  );
    SetMaterialProperty( "LAr", "SLOWCOMPONENT", LarProp->SlowScintSpectrum(), 1  );
    SetMaterialProperty( "LAr", "RINDEX",        LarProp->RIndexSpectrum(),    1  );
    SetMaterialProperty( "LAr", "ABSLENGTH",     LarProp->AbsLengthSpectrum(), CLHEP::cm );
    SetMaterialProperty( "LAr", "RAYLEIGH",      LarProp->RayleighSpectrum(),  CLHEP::cm );


    // scalar properties

    SetMaterialConstProperty("LAr", "SCINTILLATIONYIELD",  LarProp->ScintYield(true),       1/CLHEP::MeV ); // true = scaled down by prescale in larproperties
    SetMaterialConstProperty("LAr", "RESOLUTIONSCALE",     LarProp->ScintResolutionScale(), 1);
    SetMaterialConstProperty("LAr", "FASTTIMECONSTANT",    LarProp->ScintFastTimeConst(),   CLHEP::ns);
    SetMaterialConstProperty("LAr", "SLOWTIMECONSTANT",    LarProp->ScintSlowTimeConst(),   CLHEP::ns);
    SetMaterialConstProperty("LAr", "YIELDRATIO",          LarProp->ScintYieldRatio(),      1);
    SetMaterialConstProperty("LAr", "ELECTRICFIELD",       DetProp->Efield(),               CLHEP::kilovolt/CLHEP::cm);

    SetBirksConstant("LAr",LarProp->ScintBirksConstant(), CLHEP::cm/CLHEP::MeV);
    if(DetProp->SimpleBoundary())
      SetReflectances("LAr", LarProp->SurfaceReflectances(), LarProp->SurfaceReflectanceDiffuseFractions());
    else
      SetReflectances(LarProp->SurfaceReflectances());

    // If we are using scint by particle type, load these

    if(LarProp->ScintByParticleType())
      {
        // true = scaled down by prescale in larproperties
        SetMaterialConstProperty("LAr", "PROTONSCINTILLATIONYIELD",  LarProp->ProtonScintYield(true),    1./CLHEP::MeV );
        SetMaterialConstProperty("LAr", "PROTONYIELDRATIO",          LarProp->ProtonScintYieldRatio(),   1.);
        SetMaterialConstProperty("LAr", "MUONSCINTILLATIONYIELD",    LarProp->MuonScintYield(true),      1./CLHEP::MeV );
        SetMaterialConstProperty("LAr", "MUONYIELDRATIO",            LarProp->MuonScintYieldRatio(),     1.);
        SetMaterialConstProperty("LAr", "KAONSCINTILLATIONYIELD",    LarProp->KaonScintYield(true),      1./CLHEP::MeV );
        SetMaterialConstProperty("LAr", "KAONYIELDRATIO",            LarProp->KaonScintYieldRatio(),     1.);
        SetMaterialConstProperty("LAr", "PIONSCINTILLATIONYIELD",    LarProp->PionScintYield(true),      1./CLHEP::MeV );
        SetMaterialConstProperty("LAr", "PIONYIELDRATIO",            LarProp->PionScintYieldRatio(),     1.);
        SetMaterialConstProperty("LAr", "ELECTRONSCINTILLATIONYIELD",LarProp->ElectronScintYield(true),  1./CLHEP::MeV );
        SetMaterialConstProperty("LAr", "ELECTRONYIELDRATIO",        LarProp->ElectronScintYieldRatio(), 1.);
        SetMaterialConstProperty("LAr", "ALPHASCINTILLATIONYIELD",   LarProp->AlphaScintYield(true),     1./CLHEP::MeV );
        SetMaterialConstProperty("LAr", "ALPHAYIELDRATIO",           LarProp->AlphaScintYieldRatio(),    1.);
      }

    // If we are simulating the TPB load this

    if(LarProp->ExtraMatProperties())
      {
        SetMaterialProperty("TPB", "RINDEX",       LarProp->RIndexSpectrum(), 1 );
        SetMaterialProperty("TPB", "WLSABSLENGTH", LarProp->TpbAbs(),         CLHEP::m );
        SetMaterialProperty("TPB", "WLSCOMPONENT", LarProp->TpbEm(),          1 );

        SetMaterialConstProperty("TPB", "WLSTIMECONSTANT", LarProp->TpbTimeConstant(), CLHEP::ns );
      }

  }

}
