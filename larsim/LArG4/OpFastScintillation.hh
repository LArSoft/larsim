// Class adapted for LArSoft by Ben Jones, MIT 10/10/12
//
// This class is a physics process based on the standard Geant4
// scintillation process.
//
// It has been stripped down and adapted to form the backbone of
// the LArG4 fast optical simulation.  Photons, instead of being
// produced and added to the geant4 particle stack, are logged
// and used to predict the visibility of this step to each PMT in
// the detector.
//
// The photonvisibilityservice looks up the visibility of the relevant
// xyz point, and if a photon is detected at a given PMT, one OnePhoton
// object is logged in the OpDetPhotonTable
//
// At the end of the event, the OpDetPhotonTable is read out
// by LArG4, and detected photons are stored in the event.
//
// This process can be used alongside the standard Cerenkov process,
// which does step geant4 opticalphotons.  Both the fast scintillation
// table and the geant4 sensitive detectors are read out by LArG4 to
// produce a combined SimPhoton collection.
//
//
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:        OpFastScintillation.hh  
// Description:	Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//                         of energy deposited by particle type
//                         Thanks to Zach Hartwig (Department of Nuclear
//                         Science and Engineeering - MIT)
//              2005-07-28 add G4ProcessType to constructor
//              2002-11-21 change to user G4Poisson for small MeanNumPotons
//              2002-11-07 allow for fast and slow scintillation
//              2002-11-05 make use of constant material properties
//              2002-05-16 changed to inherit from VRestDiscreteProcess
//              2002-05-09 changed IsApplicable method
//              1999-10-29 add method and class descriptors
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef OpFastScintillation_h
#define OpFastScintillation_h 1

/////////////
// Includes
/////////////

#include "Geant4/globals.hh"
#include "Geant4/templates.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleMomentum.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4VRestDiscreteProcess.hh"
#include "Geant4/G4OpticalPhoton.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4Material.hh" 
#include "Geant4/G4PhysicsTable.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4PhysicsOrderedFreeVector.hh"
#include "Geant4/G4EmSaturation.hh"

#include "fhiclcpp/ParameterSet.h"
#include "TF1.h"
#include "TVector3.h"

// Class Description:
// RestDiscrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VRestDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

namespace larg4{

class OpFastScintillation : public G4VRestDiscreteProcess
{

private:

        //////////////
        // Operators
        //////////////
  
        // OpFastScintillation& operator=(const OpFastScintillation &right);
  
public: // Without description
  
	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

        OpFastScintillation(const G4String& processName = "Scintillation", G4ProcessType type = fElectromagnetic);  
        OpFastScintillation(const OpFastScintillation &right);

	~OpFastScintillation();	

        ////////////
        // Methods
        ////////////

public: // With description

        // OpFastScintillation Process has both PostStepDoIt (for energy 
        // deposition of particles in flight) and AtRestDoIt (for energy
        // given to the medium by particles at rest)

        virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable', for any particle type except
        // for an 'opticalphoton' and for short-lived particles

	G4double GetMeanFreePath(const G4Track& aTrack,
				       G4double ,
                                       G4ForceCondition* );
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'StronglyForced' condition for the DoIt to be 
        // invoked at every step.

 
        G4double GetMeanLifeTime(const G4Track& aTrack,
                                 G4ForceCondition* );
        // Returns infinity; i. e. the process does not limit the time,
        // but sets the 'StronglyForced' condition for the DoIt to be
        // invoked at every step.

	virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
			                const G4Step&  aStep);
        virtual G4VParticleChange* AtRestDoIt (const G4Track& aTrack,
                                       const G4Step& aStep);

        // These are the methods implementing the scintillation process.

	void SetTrackSecondariesFirst(const G4bool state);
        // If set, the primary particle tracking is interrupted and any
        // produced scintillation photons are tracked next. When all 
        // have been tracked, the tracking of the primary resumes.

        void SetFiniteRiseTime(const G4bool state);
        // If set, the OpFastScintillation process expects the user to have
        // set the constant material property FAST/SLOWSCINTILLATIONRISETIME.

        G4bool GetTrackSecondariesFirst() const;
        // Returns the boolean flag for tracking secondaries first.

        G4bool GetFiniteRiseTime() const;
        // Returns the boolean flag for a finite scintillation rise time.
	
        void SetScintillationYieldFactor(const G4double yieldfactor);
        // Called to set the scintillation photon yield factor, needed when
        // the yield is different for different types of particles. This
        // scales the yield obtained from the G4MaterialPropertiesTable.

        G4double GetScintillationYieldFactor() const;
        // Returns the photon yield factor.

        void SetScintillationExcitationRatio(const G4double excitationratio);
        // Called to set the scintillation exciation ratio, needed when
        // the scintillation level excitation is different for different
        // types of particles. This overwrites the YieldRatio obtained
        // from the G4MaterialPropertiesTable.

        G4double GetScintillationExcitationRatio() const;
        // Returns the scintillation level excitation ratio.

        G4PhysicsTable* GetFastIntegralTable() const;
        // Returns the address of the fast scintillation integral table.

        G4PhysicsTable* GetSlowIntegralTable() const;
        // Returns the address of the slow scintillation integral table.

        void AddSaturation(G4EmSaturation* sat) { emSaturation = sat; }
        // Adds Birks Saturation to the process.

        void RemoveSaturation() { emSaturation = NULL; }
        // Removes the Birks Saturation from the process.

        G4EmSaturation* GetSaturation() const { return emSaturation; }
        // Returns the Birks Saturation.

        void SetScintillationByParticleType(const G4bool );
        // Called by the user to set the scintillation yield as a function
        // of energy deposited by particle type

        G4bool GetScintillationByParticleType() const
        { return scintillationByParticleType; }
        // Return the boolean that determines the method of scintillation
        // production

        void DumpPhysicsTable() const;
        // Prints the fast and slow scintillation integral tables.

        std::vector<double> GetVUVTime(double, int);
        std::vector<double> GetVisibleTimeOnlyCathode(double, int);
  	// old timings -- to be deleted       
  
        std::vector<double> getVUVTime(double, int);
	void generateparam(int index);
        // Functions for vuv component Landau + Exponential timing parameterisation, updated method
        
	std::vector<double> getVISTime(TVector3 ScintPoint, TVector3 OpDetPoint, int Nphotons);
	// Visible component timing parameterisation

        int VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type);
	// Calculates semi-analytic model number of hits for vuv component

 	int VISHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type);
	// Calculates semi-analytic model number of hits for visible component

protected:

        void BuildThePhysicsTable();
        // It builds either the fast or slow scintillation integral table; 
        // or both. 

  
        bool RecordPhotonsProduced(const G4Step& aStep, double N);
        // Note the production of N photons in at point xyz.
	//  pass on to generate detector response, etc.


        ///////////////////////
        // Class Data Members
        ///////////////////////

        G4PhysicsTable* theSlowIntegralTable;
        G4PhysicsTable* theFastIntegralTable;

        G4bool fTrackSecondariesFirst;
        G4bool fFiniteRiseTime;

        G4double YieldFactor;

        G4double ExcitationRatio;

        G4bool scintillationByParticleType;

private:

        G4double single_exp(G4double t, G4double tau2) const;
        G4double bi_exp(G4double t, G4double tau1, G4double tau2) const;

        G4double scint_time(const G4Step& aStep,
                            G4double ScintillationTime,
                            G4double ScintillationRiseTime) const;
  std::vector<double> propagation_time(G4ThreeVector x0, int OpChannel, int NPhotons, bool Reflected=false); //const;
   
        // emission time distribution when there is a finite rise time
        G4double sample_time(G4double tau1, G4double tau2) const;

        // Facility for TPB emission energies
        double reemission_energy() const;
        std::map<double,double> tpbemission;
        CLHEP::RandGeneral *rgen0;

        void average_position(G4Step const& aStep, double *xzyPos) const;
  
        G4EmSaturation* emSaturation;
        // functions and parameters for the propagation time parametrization
        TF1* ParPropTimeTF1;
        float const* ReflT0s;

        TF1 const* functions_vuv[8];
        TF1 const* functions_vis[5];                     
        double fd_break;
        double fd_max;
        double ftf1_sampling_factor;
        double ft0_max, ft0_break_point;

        //For new VUV time parametrization
        double fstep_size, fmax_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance;
        std::vector<double> fparameters[9];
        // vector containing generated VUV timing parameterisations
        std::vector<TF1> VUV_timing;
        // vector containing min and max range VUV timing parameterisations are sampled to
        std::vector<double> VUV_max;
        std::vector<double> VUV_min;
       
	// For new VIS time parameterisation
	double fvis_vmean, fn_LAr_vis, fn_LAr_vuv;
	std::vector<double> fdistances_refl;
        std::vector<std::vector<double>> fcut_off_pars;
	std::vector<std::vector<double>> ftau_pars;
 
	//For VUV semi-analytic hits
	//array of correction for the VUV Nhits estimation
	std::vector<std::vector<double> > fGHvuvpars;
	TF1* GHvuv[9];

	// For VIS semi-analytic hits
	// array of corrections for VIS Nhits estimation
	std::vector<std::vector<double>> fvispars;
	TF1* VIS_pol[9];
	double fplane_depth, fcathode_width, fcathode_height;
	std::vector<double>  fcathode_centre;
        
	// Optical detector properties for semi-analytic hits
	int foptical_detector_type; 
        double fheight, fwidth, fradius;
        int fdelta_angulo, fL_abs_vuv;
        std::vector<std::vector<double> > fOpDetCenter;

        //double fGlobalTimeOffset;  

        void ProcessStep( const G4Step& step);
        
        bool bPropagate; ///< Whether propagation of photons is enabled.
        
};

  double finter_d(double*, double*);
  double LandauPlusExpoFinal(double*, double*);
  //For new VUV time parametrization    
  double interpolate( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate );
  double* interpolate( std::vector<double> &xData, std::vector<double> &yData1, std::vector<double> &yData2,
		       std::vector<double> &yData3, double x, bool extrapolate);
  double model_close(double*, double*);
  double model_far(double*, double*);
  // gaisser-hillas function
  double GaisserHillas(double *x, double *par);
  // structure definition for solid angle of rectangle function                                                                                                                                     
  struct acc{
    // ax,ay,az = centre of rectangle; w = width; h = height                                                                                                                                
    double ax, ay, az, w, h;
  };
  // solid angle of rectangular aperture calculation functions
  double Rectangle_SolidAngle(double a, double b, double d);
  double Rectangle_SolidAngle(acc& out, TVector3 v);

  // solid angle of circular aperture calculation functions
  double Disk_SolidAngle(double *x, double *p);
  double Disk_SolidAngle(double d, double h, double b);

////////////////////
// Inline methods
////////////////////

inline 
G4bool OpFastScintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
       if (aParticleType.GetParticleName() == "opticalphoton") return false;
       if (aParticleType.IsShortLived()) return false;

       return true;
}

inline 
void OpFastScintillation::SetTrackSecondariesFirst(const G4bool state) 
{
	fTrackSecondariesFirst = state;
}

inline
void OpFastScintillation::SetFiniteRiseTime(const G4bool state)
{
        fFiniteRiseTime = state;
}

inline
G4bool OpFastScintillation::GetTrackSecondariesFirst() const
{
        return fTrackSecondariesFirst;
}

inline 
G4bool OpFastScintillation::GetFiniteRiseTime() const
{
        return fFiniteRiseTime;
}

inline
void OpFastScintillation::SetScintillationYieldFactor(const G4double yieldfactor)
{
        YieldFactor = yieldfactor;
}

inline
G4double OpFastScintillation::GetScintillationYieldFactor() const
{
        return YieldFactor;
}

inline
void OpFastScintillation::SetScintillationExcitationRatio(const G4double excitationratio)
{
        ExcitationRatio = excitationratio;
}

inline
G4double OpFastScintillation::GetScintillationExcitationRatio() const
{
        return ExcitationRatio;
}

inline
G4PhysicsTable* OpFastScintillation::GetSlowIntegralTable() const
{
        return theSlowIntegralTable;
}

inline
G4PhysicsTable* OpFastScintillation::GetFastIntegralTable() const
{
        return theFastIntegralTable;
}

inline
void OpFastScintillation::DumpPhysicsTable() const
{
        if (theFastIntegralTable) {
           G4int PhysicsTableSize = theFastIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
        	v = (G4PhysicsOrderedFreeVector*)(*theFastIntegralTable)[i];
        	v->DumpValues();
           }
         }

        if (theSlowIntegralTable) {
           G4int PhysicsTableSize = theSlowIntegralTable->entries();
           G4PhysicsOrderedFreeVector *v;

           for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
           {
                v = (G4PhysicsOrderedFreeVector*)(*theSlowIntegralTable)[i];
                v->DumpValues();
           }
         }
}

inline
G4double OpFastScintillation::single_exp(G4double t, G4double tau2) const
{
         return std::exp(-1.0*t/tau2)/tau2;
}

inline
G4double OpFastScintillation::bi_exp(G4double t, G4double tau1, G4double tau2) const
{
         return std::exp(-1.0*t/tau2)*(1-std::exp(-1.0*t/tau1))/tau2/tau2*(tau1+tau2);
}

} //namespace

#endif /* OpFastScintillation_h */
