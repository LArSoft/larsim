////////////////////////////////////////////////////////////////////////
/// \file larsim/LArG4/MaterialPropertyLoader.h
//
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////
// Class to set material properties for different materials in
// the detector. Currently used to supply optical properties to LAr
// and other optical elements of the detector.
//

// TODO uniform the header guard format to LArSoft common
#ifndef LArG4_MaterialPropertyLoader_h
#define LArG4_MaterialPropertyLoader_h

#include <map>
#include <string>

class G4LogicalVolumeStore;

namespace larg4 {

  /**
   * @brief Stores material properties and sends them to GEANT4 geometry.
   *
   * Class to set material properties for different materials in the detector.
   * Currently mainly used to set optical properties for LAr and other optical
   * components.
   *
   * Notes on the implementation of reflectivity
   * --------------------------------------------
   *
   * The reflectivity properties of the material are stored in a different way
   * depending on whether the standard GEANT4 process or the custom simple
   * boundary process are used in the simulation (determined by
   * `detinfo::DetectorProperties::SimpleBoundary()`).
   * For the former, reflectivity is saved as a property of the material itself
   * with the property name `"REFLECTIVITY"`.
   * For the latter, the properties (both a reflectivity and a diffuse
   * reflection fraction) are stored as properties of the `"LAr"` material,
   * with a name `"REFLECTANCE_<Material>"` and
   * `"DIFFUSE_REFLECTANCE_FRACTION_<Material>"`. This is the storage policy
   * expected by `larg4::OpBoundaryProcessSimple`, which implements the simple
   * boundary model.
   *
   */
  class MaterialPropertyLoader
  {
  public:

    // TODO use type aliases

    // TODO remove default constructor
    MaterialPropertyLoader() {}
    // TODO remove default destructor
    ~MaterialPropertyLoader() {}


    // TODO remove duplicate "public" label
  public:

    // TODO turn arguments into constant references
    //Accessors
    std::map<double,double> GetMaterialProperty(std::string Material,std::string Property)
      {return fPropertyList[Material][Property];}

    // TODO turn arguments into constant references
    double GetMaterialConstProperty(std::string Material, std::string Property)
      {return fConstPropertyList[Material][Property];}

    // TODO turn argument into constant reference
    std::map<std::string,double> GetMaterialConstProperties(std::string Material)
      {return fConstPropertyList[Material];}

    // TODO turn argument into constant reference
    std::map<std::string,std::map<double,double> >  GetMaterialProperties(std::string Material)
      {return fPropertyList[Material];}


    // --- BEGIN Methods to set material properties ----------------------------
    /// @name Methods to set material properties
    /// @{

    /**
     * @brief Stores the specified emergy-dependent material property.
     * @param Material name of the material to set the property of
     * @param Property name of the property
     * @param Values table of property values (see below)
     * @param Unit unit of the property values (CLHEP)
     * @see `SetMaterialConstProperty()`
     *
     * The property is stored internally and _not_ propagated to GEANT4
     * (use `UpdateGeometry()` to that purpose).
     * The previous value of the property is silently overwritten.
     *
     * The table of values is in form of (`energy`, `value`) pairs, where
     * `value` is measured in `Units` and `energy` is measured in electronvolt.
     */
    // TODO turn arguments into constant references
    void SetMaterialProperty(       std::string Material, std::string Property, std::map<double,double> Values, double Unit);

    /**
     * @brief Stores the specified material property.
     * @param Material name of the material to set the property of
     * @param Property name of the property
     * @param Value the value of the property
     * @param Unit unit of the property value (CLHEP)
     * @see `SetMaterialProperty()`
     *
     * The property is stored internally and _not_ propagated to GEANT4
     * (use `UpdateGeometry()` to that purpose).
     * The previous value of the property is silently overwritten.
     */
    // TODO turn arguments into constant references
    void SetMaterialConstProperty(  std::string Material, std::string Property, double Value,                   double Unit);

    /// @}
    // --- END Methods to set material properties ------------------------------

    // --- BEGIN Setting of specific properties --------------------------------
    /// @name Setting of specific properties
    /// @{

    // Method to set LArG4 Birks constant
    void SetBirksConstant( std::string, double, double );


    /**
     * @brief
     */
    // TODO turn the arguments into constant references (using type aliases)
    void SetReflectances( std::string, std::map<std::string, std::map<double,double> >, std::map<std::string, std::map<double, double> >);
    // TODO turn the argument into constant reference (using type aliases)
    void SetReflectances( std::map<std::string, std::map<double,double> >);

    /// @}
    // --- END Setting of specific properties ----------------------------------


    /**
     * @brief Imports properties from LArSoft services
     *
     * The properties imported include:
     *
     * * material `"LAr"`:
     *     * fast scintillation light spectrum (`"FASTCOMPONENT"`) from `detinfo::LArProperties::FastScintSpectrum()`
     *     * slow scintillation light spectrum (`"SLOWCOMPONENT"`) from `detinfo::LArProperties::SlowScintSpectrum()`
     *     * refraction index vs. photon energy (`"RINDEX"`) from `detinfo::LArProperties::RIndexSpectrum()`
     *     * absorption length vs. photon energy (`"ABSLENGTH"`) from `detinfo::LArProperties::AbsLengthSpectrum()` _[cm]_
     *     * Rayleigh scattering vs. photon energy (`"RAYLEIGH"`) from `detinfo::LArProperties::RayleighSpectrum()` _[cm]_
     *     * scintillation yield (`"SCINTILLATIONYIELD"`) from `detinfo::LArProperties::ScintYield(true)` _[1/MeV]_
     *     * "RESOLUTIONSCALE"
     *     * fast scintillation light delay (`"FASTTIMECONSTANT"`) from `detinfo::LArProperties::ScintFastTimeConst()` _[ns]_
     *     * slow scintillation light delay (`"SLOWTIMECONSTANT"`) from `detinfo::LArProperties::ScintSlowTimeConst()` _[ns]_
     *     * scintillation yield ratio (`"YIELDRATIO"`) `detinfo::LArProperties::ScintYieldRatio()`
     *     * electric field (`"ELECTRICFIELD"`) from `detinfo::DetectorProperties::Efield()` _[kV/cm]_
     *     * scintillation Birks constant vs. photon energy, from `detinfo::LArProperties::ScintBirksConstant()` _[cm/MeV]_
     *     * if using the simple reflectivity model (`detinfo::DetectorProperties::SimpleBoundary()`), for each supported material `XXX`:
     *         * reflectivity (`"REFLECTANCE_XXX"`) from `detinfo::LArProperties::SurfaceReflectances()`
     *         * diffused reflectivity fraction (`"DIFFUSE_REFLECTANCE_FRACTION_XXX"`) from `detinfo::LArProperties::SurfaceReflectances()`
     *         .
     *       The materials are the ones listed in the property from `detinfo::LArProperties`.
     *     * if using different scintillation yield for different particle types (according to `detinfo::LArProperties::ScintByParticleType()`):
     *         * scintillation yield `"<PARTICLE>SCINTILLATIONYIELD"` from `detinfo::LArProperties::<Particle>ScintYield()` (as above)
     *         * scintillation yield ratio `"<PARTICLE>YIELDRATIO"` from `detinfo::LArProperties::<Particle>ScintYieldRatio()` (as above)
     *         .
     *       with `<PARTICLE>` being proton, muon, kaon, pion, electron and &alpha; particles.
     * * if _not_ using the simple reflectivity model (`!detinfo::DetectorProperties::SimpleBoundary()`):
     *     * material `XXX` reflectivity (`"REFLECTIVITY"`) from `detinfo::LArProperties::SurfaceReflectances()`
     *     .
     *     The materials are the ones listed in the property from `detinfo::LArProperties::SurfaceReflectances()`.
     * * `TPB` material (if simulating it, according to `detinfo::LArProperties::ExtraMatProperties()`):
     *     * refraction index vs. photon energy (`"RINDEX"`) from `detinfo::LArProperties::RIndexSpectrum()`
     *     * wavelength shifter absorbtion length vs. photon energy (`"WLSABSLENGTH"`) from `detinfo::LArProperties::TpbAbs()` _[m]_
     *     * `"WLSCOMPONENT"` from `detinfo::LArProperties::TpbEm()`
     *
     */
    void GetPropertiesFromServices();

    // TODO make this method constant
    /**
     * @brief Updates the material properties with the collected values.
     * @param lvs the store of logical volumes to be updated
     *
     * Before calling this function, properties for some materials (mostly
     * liquid argon, but not only: see e.g. `GetPropertiesFromServices()`) are
     * collected and updated.
     * This method considers all volumes in the store `lvs`. For the ones made
     * of a material we have properties for, their material properties are
     * updated to reflect the values we have collected.
     */
    void UpdateGeometry( G4LogicalVolumeStore* lvs );



  private:

    //         materials                properties            values
    std::map < std::string , std::map < std::string,double> > fConstPropertyList;

    //         materials                properties               energies  values
    std::map < std::string , std::map < std::string , std::map < double ,  double > > > fPropertyList;

    std::map<std::string, double> fBirksConstants;




  }; // clas MaterialPropertyLoader


} // namespace larg4


#endif // LArG4_MaterialPropertyLoader_h
