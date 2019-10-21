////////////////////////////////////////////////////////////////////////
// Class:       PDFastSimPAR
// Plugin Type: producer
// File:        PDFastSimPAR_module.cc
// Description:
// - acts on sim::SimEnergyDeposit from LArG4Main, 
// - simulate (fast, photon visibility service) the OpDet response to optical photons
// Input: 'sim::SimEnergyDeposit'
// Output: 'sim::OpDetBacktrackerRecord'
//Fast simulation of propagating the photons created from SimEnergyDeposits.

//This module does a fast simulation of propagating the photons created from SimEnergyDeposits,
//This simulation is done using the PhotonLibrary, which stores the visibilities of each optical channel
//with respect to each optical voxel in the TPC volume, to avoid propagating single photons using Geant4.
//At the end of this module a collection of the propagated photons either as
//'sim::OpDetBacktrackerRecord' are placed into the art event.

//The steps this module takes are:
//  - to take number of photon and the vertex information from 'sim::SimEnergyDeposits',
//  - use the PhotonLibrary (visibilities) to determine the amount of visible photons at each optical channel,
//  - visible photons: the number of photons times the visibility at the middle of the Geant4 step for a given optical channel.
//  - other photon information is got from 'sim::SimEnergyDeposits'
//  - add 'sim::OpDetBacktrackerRecord' to event
// Aug. 19 by Mu Wei
////////////////////////////////////////////////////////////////////////

// Art libraries
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"
#include "nurandom/RandomUtils/NuRandomService.h"


#include "fhiclcpp/ParameterSet.h"

// LArSoft libraries
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
#include "larsim/PhotonPropagation/PhotonVisibilityTypes.h" // phot::MappedT0s_t
#include "larsim/PhotonPropagation/ScintTimeTools/ScintTime.h"

#include "larsim/Simulation/LArG4Parameters.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"


// Random number engine
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
//#include "CLHEP/Random/RandGauss.h"


#include "TLorentzVector.h"

#include "Geant4/globals.hh"
#include "Geant4/G4EmProcessSubType.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4Poisson.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/templates.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleMomentum.hh"
#include "Geant4/G4VRestDiscreteProcess.hh"
#include "Geant4/G4OpticalPhoton.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4PhysicsTable.hh"
#include "Geant4/G4MaterialPropertiesTable.hh"
#include "Geant4/G4PhysicsOrderedFreeVector.hh"
#include "Geant4/G4EmSaturation.hh"

// support libraries
#include "cetlib_except/exception.h"

#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"
#include "Math/SpecFuncMathMore.h"
#include <cmath>

#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_3.hpp"

#include <chrono>
#include <ctime>


using namespace std;

namespace
{
    //......................................................................    
    double finter_d(double *x, double *par)
    {
        double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
        double y2 = TMath::Exp(par[3]+x[0]*par[4]);
        
        return TMath::Abs(y1 - y2);
    }
//  double LandauPlusExpoFinal(double *x, double *par)
//  {
//        // par0 = joining point
//        // par1 = Landau MPV
//        // par2 = Landau widt
//        // par3 = normalization
//        // par4 = Expo cte
//        // par5 = Expo tau
//        double y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
//        double y2 = TMath::Exp(par[4]+x[0]*par[5]);
//        if(x[0] > par[0]) y1 = 0.;
//        if(x[0] < par[0]) y2 = 0.;
//        
//        return (y1 + y2);
//  }
    
    //......................................................................    
//  double finter_r(double *x, double *par) 
//  {
//        double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
//        double y2 = par[5]*TMath::Landau(x[0],par[3],par[4]);
//        
//        return TMath::Abs(y1 - y2);
//  }
    
    double model_close(double *x, double *par)
    {
        // par0 = joining point
        // par1 = Landau MPV
        // par2 = Landau width
        // par3 = normalization
        // par4 = Expo cte
        // par5 = Expo tau
        // par6 = t_min
        
        double y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
        double y2 = TMath::Exp(par[4]+x[0]*par[5]);
        if(x[0] <= par[6] || x[0] > par[0]) y1 = 0.;
        if(x[0] < par[0]) y2 = 0.;
        
        return (y1 + y2);
    }
    
    //......................................................................    
    double model_far(double *x, double *par)
    {
        // par1 = Landau MPV
        // par2 = Landau width
        // par3 = normalization
        // par0 = t_min
        
        double y = par[3]*TMath::Landau(x[0],par[1],par[2]);
        if(x[0] <= par[0]) y = 0.;
        
        return y;
    }
    
    //======================================================================
    //   Returns interpolated value at x from parallel arrays ( xData, yData )
    //   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
    //   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
    //......................................................................    
    double interpolate( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
    {
        int size = xData.size();
        int i = 0;                                        // find left end of interval for interpolation
        if ( x >= xData[size - 2] )                       // special case: beyond right end
        {
            i = size - 2;
        }
        else
        {
            while ( x > xData[i+1] ) i++;
        }
        double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
        if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
        {
            if ( x < xL ) yR = yL;
            if ( x > xR ) yL = yR;
        }
        double dydx = ( yR - yL ) / ( xR - xL );          // gradient
        return yL + dydx * ( x - xL );                    // linear interpolation
    }
    
    //......................................................................    
    double* interpolate( std::vector<double> &xData, std::vector<double> &yData1, std::vector<double> &yData2, std::vector<double> &yData3, double x, bool extrapolate)
    {
        int size = xData.size();
        int i = 0;                                        // find left end of interval for interpolation
        if ( x >= xData[size - 2] )                       // special case: beyond right end
        {
            i = size - 2;
        }
        else
        {
            while ( x > xData[i+1] ) i++;
        }
        double xL = xData[i], xR = xData[i+1];// points on either side (unless beyond ends)
        double yL1 = yData1[i], yR1 = yData1[i+1], yL2 = yData2[i], yR2 = yData2[i+1], yL3 = yData3[i], yR3 = yData3[i+1];
        
        if ( !extrapolate )                                                  // if beyond ends of array and not extrapolating
        {
            if ( x < xL ) {yR1 = yL1; yR2 = yL2; yR3 = yL3;}
            if ( x > xR ) {yL1 = yR1; yL2 = yR2; yL3 = yR3;}
        }
        double dydx1 = ( yR1 - yL1 ) / ( xR - xL );          // gradient
        double dydx2 = ( yR2 - yL2 ) / ( xR - xL );
        double dydx3 = ( yR3 - yL3 ) / ( xR - xL );
    
        double *yy = new double[3];
        yy[0] = yL1 + dydx1 * ( x - xL );// linear interpolations
        yy[1] = yL2 + dydx2 * ( x - xL );
        yy[2] = yL3 + dydx3 * ( x - xL );
        
        return yy;
    }
    
    //......................................................................    
    // solid angle of circular aperture
    double Disk_SolidAngle(double* x, double *p) 
    {
        const double d = x[0];
        const double h = x[1];
        const double b = p[0];
        if(b <= 0. || d < 0. || h <= 0.) return 0.;
        const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
        if(d == 0) 
        {
            return 2.*TMath::Pi()*(1.-aa);
        }
        const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
        const double cc = 4*b*d/((b+d)*(b+d));
        
        if(TMath::Abs(boost::math::ellint_1(bb) - bb) < 1e-10 && TMath::Abs(boost::math::ellint_3(cc,bb) - cc) <1e-10)
        {
            throw(std::runtime_error("Problem loading ELLIPTIC INTEGRALS running Disk_SolidAngle!"));
        }
        if(d < b) 
        {
            return 2.*TMath::Pi() - 2.*aa*(boost::math::ellint_1(bb) + TMath::Sqrt(1.-cc)*boost::math::ellint_3(bb,cc));
        }
        if(d == b) 
        {
            return TMath::Pi() - 2.*aa*boost::math::ellint_1(bb);
        }
        if(d > b)
        {
            return 2.*aa*(TMath::Sqrt(1.-cc)*boost::math::ellint_3(bb,cc) - boost::math::ellint_1(bb));
        }
        return 0.;
    }
    
    //......................................................................    
    double Disk_SolidAngle(double d, double h, double b) 
    {
        double x[2] = { d, h };
        double p[1] = { b };
        
        return Disk_SolidAngle(x,p);
    }

    //......................................................................    
    // structure definition for solid angle of rectangle function
    struct acc
    {
        // ax,ay,az = centre of rectangle; w = width; h = height
        double ax, ay, az, w, h;
    };
    
    //......................................................................    
    // solid angle of rectanglular aperture
    double Rectangle_SolidAngle(double a, double b, double d)
    {
        double aa = a/(2.0*d);
        double bb = b/(2.0*d);
        double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
        return 4*std::acos(std::sqrt(aux));
    }
    
    //......................................................................    
    double Rectangle_SolidAngle(acc& out, TVector3 v)
    {
        //v is the position of the track segment with respect to
        //the center position of the arapuca window
        // arapuca plane fixed in x direction
        if( v.Y()==0.0 && v.Z()==0.0)
        {
            return Rectangle_SolidAngle(out.w,out.h,v.X());
        }
        
        if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0))
        {
            double A, B, a, b, d;
            A = std::abs(v.Y())-out.w/2.0;
            B = std::abs(v.Z())-out.h/2.0;
            a = out.w;
            b = out.h;
            d = std::abs(v.X());
            double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
            return to_return;
        }
        
        if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0))
        {
            double A, B, a, b, d;
            A = -std::abs(v.Y())+out.w/2.0;
            B = -std::abs(v.Z())+out.h/2.0;
            a = out.w;
            b = out.h;
            d = std::abs(v.X());
            double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
            return to_return;
        }
        
        if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0))
        {
            double A, B, a, b, d;
            A = std::abs(v.Y())-out.w/2.0;
            B = -std::abs(v.Z())+out.h/2.0;
            a = out.w;
            b = out.h;
            d = std::abs(v.X());
            double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
            return to_return;
        }
        
        if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0))
        {
            double A, B, a, b, d;
            A = -std::abs(v.Y())+out.w/2.0;
            B = std::abs(v.Z())-out.h/2.0;
            a = out.w;
            b = out.h;
            d = std::abs(v.X());
            double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
            return to_return;
        }
        
        // error message if none of these cases, i.e. something has gone wrong!
        std::cout << "Warning: invalid solid angle call." << std::endl;
        return 0.0;
    }   
    
    //......................................................................    
    double Gaisser_Hillas(double x,double *par) 
    {
        double X_mu_0        = par[3];
        double Normalization = par[0];
        double Diff          = par[1]-X_mu_0;
        double Term          = pow((x-X_mu_0)/Diff, Diff/par[2]);
        double Exponential   = std::exp((par[1]-x)/par[2]);
        
        return (Normalization*Term*Exponential);
    }
    
    //......................................................................        
    double Pol_5(double x, double *par)
    {
        // 5th order polynomial function
        return par[0] + par[1] * x + par[2] * pow(x,2) + par[3] * pow(x,3) + par[4] * pow(x,4) + par[5] * pow(x,5);
    }
} // namespace std

namespace phot
{
    class PDFastSimPAR : public art::EDProducer
    {
    public:
        explicit PDFastSimPAR(fhicl::ParameterSet const&);
        void produce(art::Event&) override;
        
        void Initialization();
        
        int VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type);
        int VISHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type);
        
        std::vector<double> propagationtime(G4ThreeVector x0, int OpChannel, int NPhotons, bool Reflected);
        std::vector<double> getVUVTime(double distance, int number_photons);
        std::vector<double> getVISTime(TVector3 ScintPoint, TVector3 OpDetPoint, int number_photons);
        
        void generateparam(int index);
        
        void AddOpDetBTR(std::vector< sim::OpDetBacktrackerRecord > & opbtr,
                              std::map<int, int> & ChannelMap,
                              sim::OpDetBacktrackerRecord btr);
    private:
        bool                          fDoSlowComponent;
        art::InputTag                 simTag;
        std::unique_ptr<ScintTime>    fScintTime;        // Tool to retrive timinig of scintillation        
        CLHEP::HepRandomEngine&       fPhotonEngine;
        CLHEP::HepRandomEngine&       fScintTimeEngine;
        std::map<int, int>            PDChannelToSOCMap; //Where each OpChan is.
        
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
        //To account for the border effects
        std::vector<double> fborder_corr;
        double fYactive_corner, fZactive_corner, fReference_to_corner, fYcathode, fZcathode;
        double fminx, fmaxx, fminy, fmaxy, fminz, fmaxz;
        // array of corrections for VIS Nhits estimation
        std::vector<std::vector<double>> fvispars;
        std::vector<double> fvis_border_distances_x;
        std::vector<double> fvis_border_distances_r;
        std::vector<std::vector<std::vector<double>>> fvis_border_correction;
        bool fApplyVisBorderCorrection;
        std::string fVisBorderCorrectionType;
        
        double fplane_depth;
        double fcathode_zdimension;
        double fcathode_ydimension;
        TVector3  fcathode_centre;
        
        // Optical detector properties for semi-analytic hits
        double fydimension;
        double fzdimension;
        double fradius;
        int    fdelta_angulo;
        int    fL_abs_vuv;
        std::vector<std::vector<double> > fOpDetCenter;
        std::vector<int>                  fOpDetType;
        std::vector<double>               fOpDetLength;
        std::vector<double>               fOpDetHeight;
        
        MappedFunctions_t ParPropTimeTF1;
    };
    
    //......................................................................    
    PDFastSimPAR::PDFastSimPAR(fhicl::ParameterSet const& pset)
    : art::EDProducer{pset}
    , fDoSlowComponent{pset.get<bool>("DoSlowComponent")}
    , simTag{pset.get<art::InputTag>("SimulationLabel")}
    , fScintTime{art::make_tool<ScintTime>(pset.get<fhicl::ParameterSet>("ScintTimeTool"))}    
    , fPhotonEngine(art::ServiceHandle<rndm::NuRandomService>{}->createEngine(*this, "HepJamesRandom", "photon", pset, "SeedPhoton"))
    , fScintTimeEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "scinttime", pset, "SeedScintTime"))
    {
        std::cout << "PDFastSimPAR Module Construct" << std::endl;
        
        Initialization();
        produces< std::vector<sim::SimPhotonsLite> >("par");
        produces< std::vector<sim::OpDetBacktrackerRecord> >("par");     
    }
    
    //......................................................................    
    void PDFastSimPAR::produce(art::Event& event)
    {
        std::cout << "PDFastSimPAR Module Producer" << std::endl;
        
        art::ServiceHandle<PhotonVisibilityService const> pvs;
        //unused auto const* larp = lar::providerFrom<detinfo::LArPropertiesService>();
        auto const nOpChannels = pvs->NOpChannels();
        
        CLHEP::RandPoissonQ randpoisphot{fPhotonEngine};
        CLHEP::RandFlat randflatscinttime{fScintTimeEngine};
        
        std::unique_ptr< std::vector< sim::OpDetBacktrackerRecord > > opbtr  (new std::vector<sim::OpDetBacktrackerRecord>);
        std::unique_ptr< std::vector< sim::SimPhotonsLite> >          phlit  (new std::vector<sim::SimPhotonsLite>);
        
        auto& photonLiteCollection (*phlit);
        photonLiteCollection.resize(nOpChannels);
        for (unsigned int i = 0; i < nOpChannels; i ++)
        {
            photonLiteCollection[i].OpChannel = i;
        }
        
        art::Handle< std::vector<sim::SimEnergyDeposit> > edepHandle;
        if (!event.getByLabel(simTag, edepHandle))
        {
            std::cout << "PDFastSimPAR Module Cannot getByLabel: " << simTag << std::endl;
            return;
        }
        
        art::ServiceHandle<geo::Geometry> geom;
        auto const& edeps = edepHandle;
        
        int num_points    = 0;
        int num_fastph    = 0;
        int num_slowph    = 0;
        int num_fastdp    = 0;
        int num_slowdp    = 0;
        
        for (auto const& edepi: *edeps)
        {
            num_points ++;
            
            int trackID       = edepi.TrackID();
            double nphot      = edepi.NumPhotons();
            double edeposit   = edepi.Energy()/nphot;
            double pos[3]     = {edepi.MidPointX(), edepi.MidPointY(), edepi.MidPointZ()};
            
            double nphot_fast = edepi.NumFPhotons();
            double nphot_slow = edepi.NumSPhotons();
            
            num_fastph += nphot_fast;
            num_slowph += nphot_slow;
            
//          ParPropTimeTF1       = pvs->GetTimingTF1(pos);
            TVector3 ScintPoint (pos[0], pos[1], pos[2]);
            for(size_t channel = 0; channel < nOpChannels; channel ++)
            {
                sim::OpDetBacktrackerRecord tmpbtr(channel);
                
                TVector3 OpDetPoint(fOpDetCenter.at(channel)[0], fOpDetCenter.at(channel)[1], fOpDetCenter.at(channel)[2]);
                fydimension = fOpDetLength.at(channel);
                fzdimension = fOpDetHeight.at(channel);
                
                if (nphot_fast > 0)
                {                   
                    auto n = VUVHits(nphot_fast, ScintPoint, OpDetPoint, fOpDetType.at(channel));
                    num_fastdp += n;
                    for (long i = 0; i < n; ++i) 
                    {
                        //calculates the time at which the photon was produced
                        fScintTime->GenScintTime(true, fScintTimeEngine);
                        auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                        ++ photonLiteCollection[channel].DetectedPhotons[time];
                        tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);                        
                    }
                }
                
                if ((nphot_slow > 0) && fDoSlowComponent) 
                {
                    auto n = VUVHits(nphot_slow, ScintPoint, OpDetPoint, fOpDetType.at(channel));
                    num_slowdp += n;
                    for (long i = 0; i < n; ++i) 
                    {
                        fScintTime->GenScintTime(false, fScintTimeEngine);
                        auto time = static_cast<int>(edepi.StartT() + fScintTime->GetScintTime());
                        ++ photonLiteCollection[channel].DetectedPhotons[time];
                        tmpbtr.AddScintillationPhotons(trackID, time, 1, pos, edeposit);
                    }
                }
                
                AddOpDetBTR(*opbtr, PDChannelToSOCMap, tmpbtr);
            }
        }
        
        std::cout << "Total points: " << num_points << ", total fast photons: " << num_fastph << ", total slow photons: " << num_slowph << std::endl;
        std::cout << "detected fast photons: " << num_fastdp << ", detected slow photons: " << num_slowdp << std::endl;
        
        PDChannelToSOCMap.clear();
        event.put(move(phlit), "par");
        event.put(move(opbtr), "par");
        
        return;
    }
    
    //......................................................................    
    void PDFastSimPAR::AddOpDetBTR(std::vector< sim::OpDetBacktrackerRecord > & opbtr,
                                         std::map<int, int> & ChannelMap,
                                         sim::OpDetBacktrackerRecord btr) 
    {
        int iChan = btr.OpDetNum();
        std::map<int, int>::iterator channelPosition = ChannelMap.find(iChan);
        
        if (channelPosition == ChannelMap.end() )
        {
            ChannelMap[iChan] = opbtr.size();
            opbtr.emplace_back(std::move(btr));
        }
        else
        {
            unsigned int idtest = channelPosition->second;
            auto const& timePDclockSDPsMap = btr.timePDclockSDPsMap();
            
            for(auto const& timePDclockSDP : timePDclockSDPsMap)
            {
                for(auto const& sdp : timePDclockSDP.second)
                {
                    double xyz[3] = {sdp.x, sdp.y, sdp.z};
                    opbtr.at(idtest).AddScintillationPhotons(sdp.trackID,
                                                             timePDclockSDP.first,
                                                             sdp.numPhotons,
                                                             xyz,
                                                             sdp.energy);
                }
            }
        }    
    }  
    
    //......................................................................    
    void PDFastSimPAR::Initialization()
    {
        std::cout << "PDFastSimPAR Initialization" << std::endl;
        std::cout << "Initializing the geometry of the detector." << std::endl;
        std::cout << "Simulate using semi-analytic model for number of hits." << std::endl;
        
        art::ServiceHandle<phot::PhotonVisibilityService const> pvs;        
        static art::ServiceHandle<geo::Geometry const> geo;
        
        // Find boundary of active volume      
        fminx =  1e9;
        fmaxx = -1e9;
        fminy =  1e9;
        fmaxy = -1e9;
        fminz =  1e9;
        fmaxz = -1e9;       
        for (size_t i = 0; i < geo->NTPC(); ++ i)
        {
            const geo::TPCGeo &tpc = geo->TPC(i);
            if (fminx > tpc.MinX()) fminx = tpc.MinX();
            if (fmaxx < tpc.MaxX()) fmaxx = tpc.MaxX();
            if (fminy > tpc.MinY()) fminy = tpc.MinY();
            if (fmaxy < tpc.MaxY()) fmaxy = tpc.MaxY();
            if (fminz > tpc.MinZ()) fminz = tpc.MinZ();
            if (fmaxz < tpc.MaxZ()) fmaxz = tpc.MaxZ();
        }
        std::cout << "Active volume boundaries:" << std::endl;
        std::cout << "minx: " << fminx << "  maxx: " << fmaxx << std::endl;
        std::cout << "miny: " << fminy << "  maxy: " << fmaxy << std::endl;
        std::cout << "minz: " << fminz << "  maxz: " << fmaxz << std::endl;
        
        TVector3 Cathode_centre(geo->TPC(0,0).GetCathodeCenter().X(), (fminy + fmaxy)/2, (fminz + fmaxz)/2);
        std::cout << "Cathode_centre: " <<Cathode_centre.X() <<",  " << Cathode_centre.Y() << ",  " << Cathode_centre.Z() << std::endl;
        
        for(size_t i = 0; i != pvs->NOpChannels(); i++)
        {
            double OpDetCenter_i[3];
            std::vector<double> OpDetCenter_v;
            geo->OpDetGeoFromOpDet(i).GetCenter(OpDetCenter_i);
            OpDetCenter_v.assign(OpDetCenter_i, OpDetCenter_i +3);
            fOpDetCenter.push_back(OpDetCenter_v);
            
            int type_i = -1;
            if(strcmp(geo->OpDetGeoFromOpDet(i).Shape()->IsA()->GetName(), "TGeoBBox") == 0) 
            {
                type_i = 0; // Arapucas
                fOpDetLength.push_back(geo->OpDetGeoFromOpDet(i).Length());
                fOpDetHeight.push_back(geo->OpDetGeoFromOpDet(i).Height());
            }
            else 
            {
                type_i = 1; // PMTs
                //    std::cout<<"Radio: "<<geo->OpDetGeoFromOpDet(i).RMax()<<std::endl;
                fOpDetLength.push_back(-1);
                fOpDetHeight.push_back(-1);
            }
            fOpDetType.push_back(type_i);
            
            std::cout <<"OpChannel: "<<i<<"  Optical_Detector_Type: "<< type_i <<"  APERTURE_height: " <<geo->OpDetGeoFromOpDet(i).Height()<<"  APERTURE_width: "<<geo->OpDetGeoFromOpDet(i).Length()<< std::endl;
        }
        
        if(pvs->IncludePropTime()) 
        {
            std::cout << "Using parameterisation of timings." << std::endl;
            // VUV time parapetrization
            pvs->LoadTimingsForVUVPar(fparameters, fstep_size, fmax_d, fvuv_vgroup_mean, fvuv_vgroup_max, finflexion_point_distance);
            
            // create vector of empty TF1s that will be replaces with the parameterisations that are generated as they are required
            // default TF1() constructor gives function with 0 dimensions, can then check numDim to qucikly see if a parameterisation has been generated
            int num_params = (fmax_d - 25) / fstep_size;  // for d < 25cm, no parameterisaton, a delta function is used instead
            std::vector<TF1> VUV_timing_temp(num_params,TF1());
            VUV_timing = VUV_timing_temp;
            
            // initialise vectors to contain range parameterisations sampled to in each case
            // when using TF1->GetRandom(xmin,xmax), must be in same range otherwise sampling is regenerated, this is the slow part!
            std::vector<double> VUV_empty(num_params, 0);
            VUV_max = VUV_empty;
            VUV_min = VUV_empty;            
        }
        
        // LAr absorption length in cm
        std::map<double, double> abs_length_spectrum = lar::providerFrom<detinfo::LArPropertiesService>()->AbsLengthSpectrum();
        std::vector<double> x_v, y_v;
        for(auto elem : abs_length_spectrum) 
        {
            x_v.push_back(elem.first);
            y_v.push_back(elem.second);
        }
        fL_abs_vuv =  interpolate(x_v, y_v, 9.7, false);
        
        std::cout << "UseNhitsModel: " << pvs->UseNhitsModel() << std::endl;
        
        // Load Gaisser-Hillas corrections for VUV semi-analytic hits
        std::cout << "Loading the GH corrections" << std::endl;
        pvs->LoadGHForVUVCorrection(fGHvuvpars, fborder_corr, fradius);
        
        fdelta_angulo = 10.; // angle bin size
        
        //Needed for Nhits-model border corrections (in cm)
        fYactive_corner = (fmaxy - fminy)/2;
        fZactive_corner = (fmaxz - fminz)/2;
        
        fYcathode = Cathode_centre.Y();
        fZcathode = Cathode_centre.Z();
        fReference_to_corner = sqrt(pow(fYactive_corner,2) + pow(fZactive_corner,2));
        
        std::cout << "For border corrections: " << fborder_corr[0] << "  " << fborder_corr[1] << std::endl;
        std::cout << "Photocathode-plane centre (z,y) = (" << fZcathode << ", " << fYcathode << ") and corner (z, y) = (" <<fZactive_corner << ", " << fYactive_corner << ")" << std::endl;
        std::cout << "Reference_to_corner: " << fReference_to_corner << std::endl;
    }

    //......................................................................    
    int PDFastSimPAR::VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type)
    {
        // check optical channel is in same TPC as scintillation light, if not return 0 hits
        // temporary method working for SBND, uBooNE, DUNE 1x2x6; to be replaced to work in full DUNE geometry
        // check x coordinate has same sign or is close to zero, otherwise return 0 hits
        if (((ScintPoint[0] < 0) != (OpDetPoint[0] < 0)) && std::abs(OpDetPoint[0]) > 10)
        {
            return 0;
        }
        
        //semi-analytic approach only works in the active volume
        if((ScintPoint[0] < fminx) || (ScintPoint[0] > fmaxx)
         ||(ScintPoint[1] < fminy) || (ScintPoint[1] > fmaxy)
         ||(ScintPoint[2] < fminz) || (ScintPoint[2] > fmaxz))
        {
            return 0;
        }
        
        // distance and angle between ScintPoint and OpDetPoint
        double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
        double cosine   = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
        double theta    = acos(cosine)*180./CLHEP::pi;
        
        // calculate solid angle:
        double solid_angle = 0;
        double d;
        double h;
        if (optical_detector_type == 0)  // Arapucas
        {
            // set Arapuca geometry struct for solid angle function
            acc detPoint;
            detPoint.ax = OpDetPoint[0];
            detPoint.ay = OpDetPoint[1];
            detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
            detPoint.w  = fydimension;
            detPoint.h  = fzdimension;    // width and height in cm of arapuca active window
            
            // get scintillation point coordinates relative to arapuca window centre
            TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;
            
            // calculate solid angle
            solid_angle = Rectangle_SolidAngle(detPoint, ScintPoint_rel);
        }
        else if (optical_detector_type == 1) // PMTs
        {
            // offset in z-y plane
            d = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
            // drift distance (in x)
            h =  sqrt(pow(ScintPoint[0] - OpDetPoint[0],2));
            // Solid angle of a disk
            solid_angle = Disk_SolidAngle(d, h, fradius);
        }
        else 
        {
            std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" <<std:: endl;
            return 0;
        }
        
        // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
        double hits_geo = exp(-1.*distance/fL_abs_vuv) * (solid_angle / (4*CLHEP::pi)) * Nphotons_created;
        
        // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence offset angle bin, accounting for border effects
        int j                     = (theta/fdelta_angulo);
        double z_to_corner        = abs(ScintPoint[2] - fZactive_corner) - fZactive_corner;
        double y_to_corner        = abs(ScintPoint[1]) - fYactive_corner;
        double distance_to_corner = sqrt(y_to_corner*y_to_corner + z_to_corner*z_to_corner);// in the ph-cathode plane
        double pars_ini_[4]       = {fGHvuvpars[0][j] + fborder_corr[0] * (distance_to_corner - fReference_to_corner),
                                     fGHvuvpars[1][j] + fborder_corr[1] * (distance_to_corner - fReference_to_corner),
                                     fGHvuvpars[2][j],
                                     fGHvuvpars[3][j]};
        double GH_correction      = Gaisser_Hillas(distance, pars_ini_);
        
        double hits_rec           = gRandom->Poisson( GH_correction*hits_geo/cosine );
        int    hits_vuv           = std::round(hits_rec); // round to integer value, cannot have non-integer number of hits
        
        return hits_vuv;
    }
    
    //......................................................................    
    int PDFastSimPAR::VISHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint, int optical_detector_type)
    {
        // check optical channel is in same TPC as scintillation light, if not return 0 hits
        // temporary method working for SBND, DUNE 1x2x6; to be replaced to work in full DUNE geometry
        // check x coordinate has same sign or is close to zero, otherwise return 0 hits
        if (((ScintPoint[0] < 0) != (OpDetPoint[0] < 0)) && std::abs(OpDetPoint[0]) > 10)
        {
            return 0;
        }
    
        //semi-analytic approach only works in the active volume
        if((ScintPoint[0] < fminx) || (ScintPoint[0] > fmaxx)
         ||(ScintPoint[1] < fminy) || (ScintPoint[1] > fmaxy)
         ||(ScintPoint[2] < fminz) || (ScintPoint[2] > fmaxz)) 
        {
            return 0;
        }
    
        // set plane_depth for correct TPC:
        double plane_depth;
        if (ScintPoint[0] < 0) 
        {
            plane_depth = -fplane_depth;
        }
        else
        {
            plane_depth = fplane_depth;
        }
    
        // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:
        // set cathode plane struct for solid angle function
        acc cathode_plane;
        cathode_plane.ax = plane_depth;
        cathode_plane.ay = fcathode_centre[1];
        cathode_plane.az = fcathode_centre[2];      // centre coordinates of cathode plane
        cathode_plane.w  = fcathode_ydimension;
        cathode_plane.h  = fcathode_zdimension;                                         // width and height in cm
        
        // get scintpoint coords relative to centre of cathode plane
        TVector3 cathodeCentrePoint(plane_depth,fcathode_centre[1],fcathode_centre[2]);
        TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint;
        
        // calculate solid angle of cathode from the scintillation point
        double solid_angle_cathode = Rectangle_SolidAngle(cathode_plane, ScintPoint_relative);
        
        // calculate distance between ScintPoint and hotspot
        // vast majority of hits in hotspot region directly infront of scintpoint
        // therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
        double distance_cathode = std::abs(plane_depth - ScintPoint[0]);
        double cosine_cathode   = 1;
        double theta_cathode    = 0;
        
        // calculate hits on cathode plane via geometric acceptance
        double cathode_hits_geo = exp(-1.*distance_cathode/fL_abs_vuv) * (solid_angle_cathode / (4.*CLHEP::pi)) * Nphotons_created;
        
        // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence offset angle bin
        int j                  = (theta_cathode/fdelta_angulo);          
        double  pars_ini_[4]   = {fGHvuvpars[0][j],
                                  fGHvuvpars[1][j],
                                  fGHvuvpars[2][j],
                                  fGHvuvpars[3][j]};
        double GH_correction    = Gaisser_Hillas(distance_cathode,pars_ini_);
        double cathode_hits_rec = GH_correction*cathode_hits_geo/cosine_cathode;  
    
        // 2). calculate number of these hits which reach the optical detector from the hotspot via solid angle hotspot coordinates
        TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);
    
        // get hotspot coordinates relative to detpoint
        TVector3 emission_relative = hotspot - OpDetPoint;
    
        // calculate solid angle of optical channel
        double solid_angle_detector = 0;
        
        if (optical_detector_type == 0) // rectangular aperture
        {
            // set rectangular aperture geometry struct for solid angle function
            acc detPoint;
            detPoint.ax = OpDetPoint[0];
            detPoint.ay = OpDetPoint[1];
            detPoint.az = OpDetPoint[2];    // centre coordinates of optical detector
            detPoint.w  = fydimension; detPoint.h = fzdimension;                            // width and height in cm of optical detector active window [rectangular aperture]
            // calculate solid angle
            solid_angle_detector = Rectangle_SolidAngle(detPoint, emission_relative);
        }
        else if (optical_detector_type == 1) // disk aperture
        {
            // offset in z-y plane
            double d = sqrt(pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
            // drift distance (in x)
            double h =  sqrt(pow(hotspot[0] - OpDetPoint[0],2));
            // calculate solid angle
            solid_angle_detector = Disk_SolidAngle(d, h, fradius);
        }
        else 
        {
            std::cout << "Error: Invalid optical detector type. 0 = rectangular, 1 = disk" <<std::endl;
            return 0;
        }
        
        // calculate number of hits via geometeric acceptance
        double hits_geo = (solid_angle_detector / (2.*CLHEP::pi)) * cathode_hits_rec;   // 2*pi rather than 4*pi due to presence of reflective foils (vm2000)
    
        // calculate distances and angles for application of corrections
        // distance to hotspot, from hotspot to optical detector, and angle between hotspot and optical detector
        double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
        double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
        double cosine_vis   = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
        double theta_vis    = acos(cosine_vis)*180./CLHEP::pi;
        int k               = (theta_vis/fdelta_angulo);
    
        // apply geometric correction
        double pars_ini_vis[6] = { fvispars[0][k], fvispars[1][k], fvispars[2][k], fvispars[3][k], fvispars[4][k], fvispars[5][k] };
        double geo_correction  = Pol_5(distance_vuv, pars_ini_vis);
        
        double hits_rec        = gRandom->Poisson(geo_correction*hits_geo/cosine_vis);
    
        // apply border correction
        int hits_vis = 0;
        if (fApplyVisBorderCorrection)
        {
            // calculate distance for interpolation depending on model
            double r = 0;
            if (fVisBorderCorrectionType == "Radial")
            {
                r = sqrt (pow(ScintPoint[1] - fcathode_ydimension,2) + pow (ScintPoint[2] - fcathode_zdimension,2));
            }
            else if (fVisBorderCorrectionType == "Vertical") 
            {
                r = std::abs(ScintPoint[1]);
            }
            else 
            {
                std::cout << "Invalid border correction type - defaulting to using central value" << std::endl;
            }
            // interpolate in x for each r bin
            int nbins_r = fvis_border_correction[k].size();
            std::vector<double> interp_vals(nbins_r, 0.0);
            for (int i = 0; i < nbins_r; i++)
            {
                interp_vals[i] = interpolate(fvis_border_distances_x, fvis_border_correction[k][i], std::abs(ScintPoint[0]), false);
            }
            // interpolate in r
            double border_correction = interpolate(fvis_border_distances_r, interp_vals, r, false);
            // apply border correction
            double hits_rec_borders = border_correction * hits_rec / cosine_vis;
            
            // round final result
            hits_vis = std::round(hits_rec_borders);
        }
        else 
        {
            // round final result
            hits_vis = std::round(hits_rec);
        }
    
        return hits_vis;
    }
        
    //......................................................................    
    std::vector<double> PDFastSimPAR::propagationtime(G4ThreeVector x0, int OpChannel, int NPhotons, bool Reflected)
    {
    
        static art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
        
        // Initialize vector of the right length with all 0's
        std::vector<double> arrival_time_dist(NPhotons, 0);
    
        if (pvs->IncludeParPropTime() && pvs->IncludePropTime()) 
        {
            throw cet::exception("OpFastScintillation") << "Cannot have both propagation time models simultaneously.";
        }
        else if (pvs->IncludeParPropTime() && !(ParPropTimeTF1  && (ParPropTimeTF1[OpChannel].GetNdim()==1)) )
        {
            //Warning: TF1::GetNdim()==1 will tell us if the TF1 is really defined or it is the default one.
            //This will fix a segfault when using timing and interpolation.
            G4cout << "WARNING: Requested parameterized timing, but no function found. Not applying propagation time." << G4endl;
        }
        else if (pvs->IncludeParPropTime())
        {
            if (Reflected)
            {
                throw cet::exception("OpFastScintillation") << "No parameterized propagation time for reflected light";
            }
            
            for (int i = 0; i < NPhotons; i++)
            {
                arrival_time_dist[i] = ParPropTimeTF1[OpChannel].GetRandom();
            }
        }
        else if (pvs->IncludePropTime()) 
        {
            // Get VUV photons arrival time distribution from the parametrization
            G4ThreeVector OpDetPoint(fOpDetCenter.at(OpChannel)[0]*CLHEP::cm,fOpDetCenter.at(OpChannel)[1]*CLHEP::cm,fOpDetCenter.at(OpChannel)[2]*CLHEP::cm);
            
            if (!Reflected)
            {
                double distance_in_cm = (x0 - OpDetPoint).mag()/CLHEP::cm; // this must be in CENTIMETERS!
                arrival_time_dist     = getVUVTime(distance_in_cm, NPhotons); // in ns
            }
            else 
            {
                TVector3 ScintPoint( x0[0]/CLHEP::cm, x0[1]/CLHEP::cm, x0[2]/CLHEP::cm ); // in cm
                TVector3 OpDetPoint_tv3(fOpDetCenter.at(OpChannel)[0], fOpDetCenter.at(OpChannel)[1], fOpDetCenter.at(OpChannel)[2]); // in cm
                arrival_time_dist = getVISTime(ScintPoint, OpDetPoint_tv3, NPhotons); // in ns
            }
        }
    
        return arrival_time_dist;
    }
    
    //......................................................................    
    std::vector<double> PDFastSimPAR::getVUVTime(double distance, int number_photons)
    {
        // pre-allocate memory
        std::vector<double> arrival_time_distrb;
        arrival_time_distrb.clear();
        arrival_time_distrb.reserve(number_photons);
        
        if (distance < 25)  // distance < 25cm
        {
            // times are fixed shift i.e. direct path only
            double t_prop_correction = distance/fvuv_vgroup_mean;
            for (int i = 0; i < number_photons; i++)
            {
                arrival_time_distrb.push_back(t_prop_correction);
            }
        }
        else 
        {
            // determine nearest parameterisation in discretisation
            int index = std::round((distance - 25) / fstep_size);
            // check whether required parameterisation has been generated, generating if not
            if (VUV_timing[index].GetNdim() == 0) 
            {
                generateparam(index);
            }
            // randomly sample parameterisation for each photon
            for (int i = 0; i < number_photons; i++)
            {
                arrival_time_distrb.push_back(VUV_timing[index].GetRandom(VUV_min[index],VUV_max[index]));
            }
        }
        return arrival_time_distrb;
        
    }    
    
    //......................................................................    
    std::vector<double> PDFastSimPAR::getVISTime(TVector3 ScintPoint, TVector3 OpDetPoint, int number_photons)
    {
        // *************************************************************************************************
        //     Calculation of earliest arrival times and corresponding unsmeared distribution
        // *************************************************************************************************
        
        // set plane_depth for correct TPC:
        double plane_depth;
        if (ScintPoint[0] < 0) 
        {
            plane_depth = -fplane_depth;
        }
        else 
        {
            plane_depth = fplane_depth;
        }
        
        // calculate point of reflection for shortest path accounting for difference in refractive indicies vectors for storing results
        TVector3 image(0,0,0);
        TVector3 bounce_point(0,0,0);
        
        // distance to wall
        TVector3 v_to_wall(plane_depth-ScintPoint[0],0,0);
        
        // hotspot is point on wall where TPB is activated most intensely by the scintillation
        TVector3 hotspot(plane_depth,ScintPoint[1],ScintPoint[2]);
        
        // define "image" by reflecting over plane
        image = hotspot + v_to_wall*(fn_LAr_vis/fn_LAr_vuv);
        
        // find point of intersection with plane j of ray from the PMT to the image
        TVector3 tempvec = (OpDetPoint-image).Unit();
        double tempnorm  = ((image-hotspot).Mag())/std::abs(tempvec[0]);
        bounce_point     = image + tempvec*tempnorm;
        
        // calculate distance travelled by VUV light and by vis light
        double VUVdist = (bounce_point-ScintPoint).Mag();
        double Visdist = (OpDetPoint-bounce_point).Mag();
        
        // calculate times taken by each part
        std::vector<double> VUVTimes  = getVUVTime(VUVdist, number_photons);
        std::vector<double> ReflTimes(number_photons,Visdist/fvis_vmean);
        
        // sum parts to get total transport times times
        std::vector<double> transport_time_vis(number_photons,0);
        for (int i=0; i<number_photons; i++) 
        {
            transport_time_vis[i] = VUVTimes[i] + ReflTimes[i];
        }
        
        // *************************************************************************************************
        //      Smearing of arrival time distribution
        // *************************************************************************************************
        
        // calculate fastest time possible
        // vis part
        double vis_time = Visdist/fvis_vmean;
        // vuv part
        double vuv_time;
        if (VUVdist < 25)
        {
            vuv_time = VUVdist/fvuv_vgroup_mean;
        }
        else 
        {
            // find index of required parameterisation
            int index = std::round((VUVdist - 25) / fstep_size);
            // find shortest time
            vuv_time  = VUV_min[index];
        }
        // sum
        double fastest_time = vis_time + vuv_time;
        
        // calculate angle alpha between scintillation point and reflection point
        double cosine_alpha = sqrt(pow(ScintPoint[0] - bounce_point[0],2)) / VUVdist;
        double alpha        = acos(cosine_alpha)*180./CLHEP::pi;
        
        // determine smearing parameters using interpolation of generated points:
        // 1). tau = exponential smearing factor, varies with distance and angle
        // 2). cutoff = largest smeared time allowed, preventing excessively large times caused by exponential
        // distance to cathode
        double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
        // angular bin
        unsigned int alpha_bin = alpha / 10;
        if (alpha_bin >= ftau_pars.size()) 
        {
            alpha_bin = ftau_pars.size() - 1;    // default to the largest available bin if alpha larger than parameterised region; i.e. last bin effectively [last bin start value, 90] deg bin
        }
        // cut-off and tau
        double cutoff = interpolate( fdistances_refl, fcut_off_pars[alpha_bin], distance_cathode_plane, true );
        double tau    = interpolate( fdistances_refl, ftau_pars[alpha_bin], distance_cathode_plane, true );
        
        // fail-safe if tau extrapolate goes wrong, drops below zero since last distance close to zero [did not occur in testing, but possible]
        if (tau < 0)
        {
            tau = 0;
        }
        
        // apply smearing:
        for (int i=0; i < number_photons; i++)
        {
            double arrival_time = transport_time_vis[i];
            double arrival_time_smeared;
            // if time is already greater than cutoff or minimum smeared time would be greater than cutoff, do not apply smearing
            if (arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(1.0))-1) >= cutoff) 
            {
                arrival_time_smeared = arrival_time;
            }
            // otherwise smear
            else 
            {
                int counter = 0;
                // loop until time generated is within cutoff limit
                // most are within single attempt, very few take more than two
                do
                {
                    // don't attempt smearings too many times for cases near cutoff (very few cases, not smearing these makes negigible difference)
                    if (counter >= 10)
                    {
                        arrival_time_smeared = arrival_time; // don't smear
                        break;
                    }
                    else 
                    {
                        // generate random number in appropriate range
                        double x = gRandom->Uniform(0.5,1.0);
                        // apply the exponential smearing
                        arrival_time_smeared = arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(x))-1);
                    }
                    // increment counter
                    counter++;
                }  while (arrival_time_smeared > cutoff);
            }
            transport_time_vis[i] = arrival_time_smeared;
        }
        
        return transport_time_vis;
    }
    
    //......................................................................    
    void PDFastSimPAR::generateparam(int index)
    {
        // get distance
        double distance_in_cm = (index * fstep_size) + 25;
    
        // time range
        const double signal_t_range = 5000.;
    
        // parameterisation TF1
        TF1 fVUVTiming;
    
        // For very short distances the time correction is just a shift
        double t_direct_mean = distance_in_cm/fvuv_vgroup_mean;
        double t_direct_min  = distance_in_cm/fvuv_vgroup_max;
    
        // Defining the model function(s) describing the photon transportation timing vs distance
        // Getting the landau parameters from the time parametrization
        double* pars_landau = interpolate(fparameters[0], fparameters[2], fparameters[3], fparameters[1], distance_in_cm, true);
        // Deciding which time model to use (depends on the distance)
        // defining useful times for the VUV arrival time shapes
        if(distance_in_cm >= finflexion_point_distance) 
        {
            double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
            // Set model: Landau
            fVUVTiming = TF1("fVUVTiming",model_far,0,signal_t_range,4);
            fVUVTiming.SetParameters(pars_far);
        }
        else 
        {
            // Set model: Landau + Exponential
            fVUVTiming = TF1("fVUVTiming",model_close,0,signal_t_range,7);
            // Exponential parameters
            double pars_expo[2];
            // Getting the exponential parameters from the time parametrization
            pars_expo[1] = interpolate(fparameters[4], fparameters[5], distance_in_cm, true);
            //For simplicity, not considering the small dependency with the offset angle in pars_expo[0]
            //Using the value for the [30,60deg] range. fparameters[6] and fparameters[8] are the values
            //for [0,30deg] range and [60,90deg] range respectively
            pars_expo[0]  = fparameters[7].at(0) + fparameters[7].at(1)*distance_in_cm;
            pars_expo[0] *= pars_landau[2];
            pars_expo[0]  = log(pars_expo[0]);
            
            // this is to find the intersection point between the two functions:
            TF1 fint          = TF1("fint",finter_d,pars_landau[0],4*t_direct_mean,5);
            double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
            fint.SetParameters(parsInt);
            double t_int  = fint.GetMinimumX();
            double minVal = fint.Eval(t_int);
            // the functions must intersect - output warning if they don't
            if(minVal>0.015) 
            {
                std::cout<<"WARNING: Parametrization of VUV light discontinuous for distance = " << distance_in_cm << std::endl;
                std::cout<<"WARNING: This shouldn't be happening " << std::endl;
            }
            double parsfinal[7] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1], t_direct_min};
            fVUVTiming.SetParameters(parsfinal);
            delete pars_landau;
        }
    
        // set the number of points used to sample parameterisation
        // for shorter distances, peak is sharper so more sensitive sampling required
        int f_sampling;
        if (distance_in_cm < 50)
        {
            f_sampling = 10000; 
        }
        else if (distance_in_cm < 100)
        {
            f_sampling = 5000; 
        }
        else
        {
             f_sampling = 1000; 
        }
        fVUVTiming.SetNpx(f_sampling);
    
        // calculate max and min distance relevant to sample parameterisation
        // max
        const int nq_max=1;
        double xq_max[nq_max];
        double yq_max[nq_max];
        xq_max[0] = 0.99;   // include 99%
        fVUVTiming.GetQuantiles(nq_max,yq_max,xq_max);
        double max = yq_max[0];
        // min
        double min = t_direct_min;
        
        // generate the sampling
        // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
        // all subsequent calls check if it has been generated previously and are ~100+ times quicker
        // add timing to the vector of timings and range to vectors of ranges
        VUV_timing[index] = fVUVTiming;
        VUV_max[index]    = max;
        VUV_min[index]    = min;
    }
} // namespace

DEFINE_ART_MODULE(phot::PDFastSimPAR)
