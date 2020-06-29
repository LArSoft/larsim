#include "DK2NuInterface.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/NuChoice.h"
#include "dk2nu/tree/calcLocationWeights.h"

#include "TFile.h"
#include "TTree.h"

#include "TRandom.h"
#include "TMath.h"

#include <iomanip>

namespace fluxr {
  DK2NuInterface::DK2NuInterface() 
  {
  }
  
  DK2NuInterface::~DK2NuInterface() 
  {
  }
  
  void DK2NuInterface::SetRootFile(TFile* rootFile)
  {
    fDk2NuTree=dynamic_cast<TTree*>(rootFile->Get("dk2nuTree"));
    fDkMetaTree=dynamic_cast<TTree*>(rootFile->Get("dkmetaTree"));

    fDk2Nu = new bsim::Dk2Nu;
    fDkMeta= new bsim::DkMeta;
    fNuChoice=new bsim::NuChoice;
    fDk2NuTree->SetBranchAddress("dk2nu",&fDk2Nu);
    fDkMetaTree->SetBranchAddress("dkmeta",&fDkMeta);
    
    fNEntries=fDk2NuTree->GetEntries();
    fDkMetaTree->GetEntry(0);
    fRun=fDkMeta->job;
    fPOT=fDkMeta->pots;
  }

  void DK2NuInterface::Init(fhicl::ParameterSet const & ps)
  {
    // Code to handle flux window adapted from
    // GENIE_R21210/src/FluxDrivers/GNuMIFlux.cxx
    if (ps.is_key_to_table("rotmatrix")) {
      // Matrix defined with series of rotations
      fhicl::ParameterSet rotps = ps.get<fhicl::ParameterSet>("rotmatrix");
      fTempRot.RotateX(rotps.get<double>("x"));
      fTempRot.RotateY(rotps.get<double>("y"));
      fTempRot.RotateZ(rotps.get<double>("z"));
    }
    else {
      std::vector<double> rotmat = ps.get<std::vector<double> >("rotmatrix");
      TVector3 newX, newY, newZ;

      if (rotmat.size() == 9) {
        // Matrix defined with new axis values
        newX = TVector3(rotmat[0], rotmat[1], rotmat[2]);
        newY = TVector3(rotmat[3], rotmat[4], rotmat[5]);
        newZ = TVector3(rotmat[6], rotmat[7], rotmat[8]);
      }
      else if (rotmat.size() == 6) {
        // Matrix defined with theta phi rotations
        newX = AnglesToAxis(rotmat[0], rotmat[1]);
        newY = AnglesToAxis(rotmat[2], rotmat[3]);
        newZ = AnglesToAxis(rotmat[4], rotmat[5]);
        fTempRot.RotateAxes(newX,newY,newZ);
        fBeamRotXML = fTempRot;
      }

      fTempRot.RotateAxes(newX, newY, newZ);
    }

    fBeamRotXML = fTempRot.Inverse();
    fBeamRot    = TLorentzRotation(fBeamRotXML);
    fBeamRotInv = fBeamRot.Inverse();

    int w=10, p=6;
    std::cout << " fBeamRotXML: " << std::setprecision(p) << std::endl;
    std::cout << " [ " 
              << std::setw(w) << fBeamRotXML.XX() << " "
              << std::setw(w) << fBeamRotXML.XY() << " "
              << std::setw(w) << fBeamRotXML.XZ() << std::endl
              << "   " 
              << std::setw(w) << fBeamRotXML.YX() << " "
              << std::setw(w) << fBeamRotXML.YY() << " "
              << std::setw(w) << fBeamRotXML.YZ() << std::endl
              << "   " 
              << std::setw(w) << fBeamRotXML.ZX() << " "
              << std::setw(w) << fBeamRotXML.ZY() << " "
              << std::setw(w) << fBeamRotXML.ZZ() << " ] " << std::endl;
    std::cout << std::endl;

    std::vector<double> detxyz = ps.get<std::vector<double> >("userbeam");
    if (detxyz.size() == 3) {
      fBeamPosXML = TVector3(detxyz[0], detxyz[1], detxyz[2]); 
    }
    else if (detxyz.size() == 6) {
      TVector3 userpos = TVector3(detxyz[0], detxyz[1], detxyz[2]); 
      TVector3 beampos = TVector3(detxyz[3], detxyz[4], detxyz[5]); 
      fBeamPosXML = userpos - fBeamRotXML*beampos;
    }
    else {
      std::cout << "userbeam needs 3 or 6 numbers to be properly defined" << std::endl;
    }

    fBeamZero=TLorentzVector(fBeamPosXML, 0);

    w=16; p=10;
    std::cout << " fBeamPosXML: [ " << std::setprecision(p) 
              << std::setw(w) << fBeamPosXML.X() << " , "
              << std::setw(w) << fBeamPosXML.Y() << " , "
              << std::setw(w) << fBeamPosXML.Z() << " ] "
              << std::endl;

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // 
    // calculation from flux window 
    // discontinued in favor of random point in active volume
    
    /*std::vector<double> windowBase=ps.get<std::vector<double> >("windowBase");
    std::vector<double> window1   =ps.get<std::vector<double> >("window1");
    std::vector<double> window2   =ps.get<std::vector<double> >("window2");

    fFluxWindowPtUser[0]=TVector3( windowBase[0], windowBase[1], windowBase[2] );
    fFluxWindowPtUser[1]=TVector3( window1[0]   , window1[1]   , window1[2] );
    fFluxWindowPtUser[2]=TVector3( window2[0]   , window2[1]   , window2[2] );
     
    // convert from user to beam coord and from 3 points to base + 2 directions
    // apply units conversion
    TLorentzVector ptbm0, ptbm1, ptbm2;
    User2BeamPos(TLorentzVector(fFluxWindowPtUser[0],0),ptbm0);
    User2BeamPos(TLorentzVector(fFluxWindowPtUser[1],0),ptbm1);
    User2BeamPos(TLorentzVector(fFluxWindowPtUser[2],0),ptbm2);

    fFluxWindowBase = ptbm0;
    fFluxWindowDir1 = ptbm1 - ptbm0;
    fFluxWindowDir2 = ptbm2 - ptbm0;

    fFluxWindowLen1 = fFluxWindowDir1.Mag();
    fFluxWindowLen2 = fFluxWindowDir2.Mag();
    fWindowNormal = fFluxWindowDir1.Vect().Cross(fFluxWindowDir2.Vect()).Unit();
    //in genie flux driver area is divided out when calculating effective POT
    //here we will keeep the POT of dk2nu file, so boost up weights 
    //convert to m^2 (since window specified in cm)
    fWindowArea   = fFluxWindowDir1.Vect().Cross(fFluxWindowDir2.Vect()).Mag()/10000.;

    double dot = fFluxWindowDir1.Dot(fFluxWindowDir2);
    if ( TMath::Abs(dot) > 1.0e-8 ) 
      std::cout << "Dot product between window direction vectors was "
        << dot << "; please check for orthoganality"<<std::endl;
    */

    //
    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    
    // calculation from random point in detector 

    std::vector<double> x_detAV = ps.get<std::vector<double> >("x_detAV"); 
    std::vector<double> y_detAV = ps.get<std::vector<double> >("y_detAV");
    std::vector<double> z_detAV = ps.get<std::vector<double> >("z_detAV");

    // assign random point  
 
    detAV_rand_user[0] = gRandom->Uniform(x_detAV[0], x_detAV[1]);
    detAV_rand_user[1] = gRandom->Uniform(y_detAV[0], y_detAV[1]);
    detAV_rand_user[2] = gRandom->Uniform(z_detAV[0], z_detAV[1]); 
    
    // convert to beam coordinates 
    fRandUser = TLorentzVector(detAV_rand_user, 0); 
    User2BeamPos(fRandUser, fRandBeam);

    // this will serve as the input point for calcEnuWgt function 

    //////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    //overwrite dkmeta
    fDkMeta->location.clear();
    fDkMeta->location.resize(2);
    fDkMeta->location[0].x=0;
    fDkMeta->location[0].y=0;
    fDkMeta->location[0].z=0;
    TLorentzVector detVec;
    User2BeamPos(TLorentzVector(0,0,0,0),detVec);
    fDkMeta->location[1].x=detVec.Vect().X();
    fDkMeta->location[1].y=detVec.Vect().Y();
    fDkMeta->location[1].z=detVec.Vect().Z();
    
    std::cout<<"Locations "<<std::endl;
    for (unsigned int i=0;i<fDkMeta->location.size();i++) {
      std::cout<<i<<"\t"<<fDkMeta->location[i].x<<"\t"
           <<fDkMeta->location[i].y<<"\t"
           <<fDkMeta->location[i].z<<std::endl;
    }
  }
  void DK2NuInterface::User2BeamPos(const TLorentzVector& usrxyz,
                                   TLorentzVector& beamxyz) const
  {
    beamxyz = fBeamRotInv*(usrxyz-fBeamZero);
  }
  void DK2NuInterface::Beam2UserPos(const TLorentzVector& beamxyz, 
                    TLorentzVector& usrxyz) const
  {
    usrxyz = fBeamRot*beamxyz + fBeamZero;
  }
  void DK2NuInterface::Beam2UserP4 (const TLorentzVector& beamp4, 
                    TLorentzVector& usrp4 ) const
  {
    usrp4 = fBeamRot*beamp4;
  }
  
    TVector3 DK2NuInterface::AnglesToAxis(double theta, double phi)
  {
    //theta,phi in rad
    double xyz[3];
    xyz[0] = TMath::Cos(phi)*TMath::Sin(theta);
    xyz[1] = TMath::Sin(phi)*TMath::Sin(theta);
    xyz[2] = TMath::Cos(theta);
    // condition vector to eliminate most floating point errors
    for (int i=0; i<3; ++i) {
      const double eps = 1.0e-15;
      if (TMath::Abs(xyz[i])   < eps ) xyz[i] =  0;
      if (TMath::Abs(xyz[i]-1) < eps ) xyz[i] =  1;
      if (TMath::Abs(xyz[i]+1) < eps ) xyz[i] = -1;
    }
    return TVector3(xyz[0],xyz[1],xyz[2]);                    
  }

  bool DK2NuInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux)

  {
    if (!fDk2NuTree->GetEntry(ientry))
      return false;

    // error handling 
    // Veto the neutrons and unphysical muon decay
    // Returns default mcflux with nu id = -999
    if(fDk2Nu->decay.ptype == 2112) return true;
    if ( (fDk2Nu->decay.ptype == 13 || fDk2Nu->decay.ptype == -13) && (fDk2Nu->decay.ntype == 14 || fDk2Nu->decay.ntype == -14 )){
      // Calculate xnu
    	double xnu = (2 * fDk2Nu->decay.necm)/0.105658389; // 
        if (xnu > 1) return true;
    }

    // bypass flux window method in favor of random point in detector
    //TLorentzVector x4beam=fFluxWindowBase+fRnd.Uniform()*fFluxWindowDir1+fRnd.Uniform()*fFluxWindowDir2;
    TLorentzVector x4beam = fRandBeam; 

    // enu = resulting energy when the neutrino is redirected 
    // wgt = output probability weight 
    double enu,wgt;
    bsim::calcEnuWgt(fDk2Nu, x4beam.Vect(),enu,wgt);

    TLorentzVector x4usr;
    Beam2UserPos(x4beam,x4usr);

    bsim::NuRay rndnuray=fDk2Nu->nuray[0];    
    fDk2Nu->nuray.clear();

    TVector3 xyzDk(fDk2Nu->decay.vx,fDk2Nu->decay.vy,fDk2Nu->decay.vz);  // origin of decay
    TVector3 p3beam = enu * (x4beam.Vect()-xyzDk).Unit();

    /*
    //weight due to window being tilted with respect to beam direction
    double tiltwgt = p3beam.Unit().Dot( fWindowNormal );
    wgt*=tiltwgt;
    //weight for the window area and divide by pi (since wgt returned by calcEnuWgt function is flux/(pi*m^2)
    wgt*=fWindowArea/3.14159;
    */

    // divide output wgt by pi so it is in unit area (output wgt is returned as flux/(pi*m^2))
    wgt *= 1/TMath::Pi(); 

    // with the recalculated energy, compute the momentum components
    bsim::NuRay anuray(p3beam.x(), p3beam.y(), p3beam.z(), enu, wgt);    
    fDk2Nu->nuray.push_back(rndnuray);
    fDk2Nu->nuray.push_back(anuray);

    TLorentzVector p4beam(p3beam,enu);
    TLorentzVector p4usr;
    Beam2UserP4(p4beam,p4usr);

    fNuPos=TLorentzVector(x4usr);
    fNuMom=TLorentzVector(p4usr);

    x4usr.SetX(x4usr.X()/100.);  // from cm --> m 
    x4usr.SetY(x4usr.Y()/100.);
    x4usr.SetZ(x4usr.Z()/100.);

    // overwrite nuchoice    
    fNuChoice->clear();
    fNuChoice->pdgNu=fDk2Nu->decay.ntype;
    fNuChoice->xyWgt=fDk2Nu->nuray[1].wgt;
    fNuChoice->impWgt=fDk2Nu->decay.nimpwt;
    fNuChoice->x4NuBeam=x4beam;
    fNuChoice->p4NuBeam=TLorentzVector(fDk2Nu->nuray[1].px,fDk2Nu->nuray[1].py,fDk2Nu->nuray[1].pz,fDk2Nu->nuray[1].E);
    //need rotation matrix here to fill these
    fNuChoice->x4NuUser=x4usr;
    fNuChoice->p4NuUser=p4usr;
    
    flux.fntype    = fDk2Nu->decay.ntype;
    flux.fnimpwt   = fDk2Nu->decay.nimpwt;
    flux.fvx       = fDk2Nu->decay.vx;
    flux.fvy       = fDk2Nu->decay.vy;
    flux.fvz       = fDk2Nu->decay.vz;
    flux.fpdpx     = fDk2Nu->decay.pdpx;
    flux.fpdpy     = fDk2Nu->decay.pdpy;
    flux.fpdpz     = fDk2Nu->decay.pdpz;
    flux.fppdxdz   = fDk2Nu->decay.ppdxdz;
    flux.fppdydz   = fDk2Nu->decay.ppdydz;
    flux.fpppz     = fDk2Nu->decay.pppz;   
    flux.fppmedium = fDk2Nu->decay.ppmedium;
    flux.fptype    = fDk2Nu->decay.ptype;
    flux.fndecay   = fDk2Nu->decay.ndecay;
    flux.fmuparpx  = fDk2Nu->decay.muparpx;
    flux.fmuparpy  = fDk2Nu->decay.muparpy;
    flux.fmuparpz  = fDk2Nu->decay.muparpz;
    flux.fmupare   = fDk2Nu->decay.mupare;
    flux.fnecm     = fDk2Nu->decay.necm;

    flux.ftpx      = fDk2Nu->tgtexit.tpx;
    flux.ftpy      = fDk2Nu->tgtexit.tpy;
    flux.ftpz      = fDk2Nu->tgtexit.tpz;
    flux.ftptype   = fDk2Nu->tgtexit.tptype;
    flux.ftgen     = fDk2Nu->tgtexit.tgen;

    flux.frun      = fDk2Nu->job;
    flux.fevtno    = fDk2Nu->potnum;
    flux.fnenergyn = flux.fnenergyf = enu;
    flux.fnwtnear  = flux.fnwtfar = wgt; 
    flux.fdk2gen   = (x4beam.Vect()-xyzDk).Mag();
    flux.ftgptype  = fDk2Nu->ancestor[1].pdg;
    
    return true;
  }
}

