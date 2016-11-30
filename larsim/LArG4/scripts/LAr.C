//-----------------------------------------------------------------------------------------------------------------
// Formulas and numbers are based on 
// arXiv:1502.04213
// and private communication with one of the authors Emily Grace:
// emilygrace.k@gmail.com
//
// to run the script from the root command line do:
// root [0] .L LAr.C++
// root [1] init();                         // initialize
// to make some plots do:
// root [2] sellmeierLAr();                 // plot refraction index between 110 and 700 nm
// root [3] sellmeierLAr(125.,130.);        // plot refraction index between 125 and 130 nm
// root [4] rayleigh();                     // plot rayleigh length between 110 and 400 nm
// root [5] rayleigh(125,130);              // plot rayleigh length between 110 and 400 nm
//
// to print out table that can be cut and paste into the Geant4 gdml description:
// root [6] rindextable();                   // refraction index between 110 and 700 nm in nsteps=100 for T=83.81 K
// root [7] rindextable(200,300,10,3);       // refraction index between 200 and 300 nm in nsteps=10 for T=90 K
// root [8] rayleightable();                 // rayleigh length between 110 and 700 nm in nsteps=100 for T=83.81 K
// root [9] rayleightable(200,300,10,3);     // rayleigh length between 200 and 300 nm in nsteps=10 for T=90 K
//------------------------------------------------------------------------------------------------------------------
#include "math.h"
#include <iostream>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TF1.h"
using namespace std;
const double pi        = 3.14159265358979323846;
const double p3times16 = 16.*pi*pi*pi;
const double kA        = 6.02214129e23;       // Avogadro constant
const double kb        = 1.3806488e-16;       // boltzmanm constant in cm2 g s-2 K-1
const double kT        = 2.18e-10;            // isothermal compressibility in cm^2dyne^-1 where 1 dyne = 1 g·cm/s2
const double aw        = 39.948;              // g/mol.
const double c         = 299792458.;          // speed of light in m/sec
const double h         = 4.13566743E-15;      // Planck constant in eVsec
const double lambdaUV  = 106.6*1e-7;          // LAr UV Resonance lambda (cm)
const double lambdaIR  = 908.3*1e-7;          // LAr IR Resonance lambda (cm)
const double alngth    = 66.;                 // LAr attenuation length in cm at 128 nm
                                              // from N. Ishida et al :
                                              // Nuclear Instruments and Methods in Physics Research Section A: Accelerators,
                                              // Spectrometers, Detectors and Associated Equipment 384 (23) (1997) 380 - 386.
// refractive index at triple point LAr from 
// A. C. Sinnock, B. L. Smith, Refractive indices of the condensed inert gases, Phys. Rev. 181 (1969) 1297-1307.
const Int_t nr= 9;
const Double_t Lambda[nr]= {361.2,
			    365.0,
			    406.3,
			    435.8,
			    475.3,
			    508.6,
			    546.1,
			    578.0,
			    643.9
};
const Double_t R[nr]= {1.2395,
		       1.2392,
		       1.2372,
		       1.2361,
		       1.2349,
		       1.2341,
		       1.2334,
		       1.2328,
		       1.2321
};
// Liquid Argon 
// sellmeier coefficient from arXiv:1502.04213
// for different temperatures T
//

const double T[4]   = {83.81     ,86.        ,88.        ,90.};        // K
const double a0[4]  = {1.24262   ,1.23828    ,1.23358    ,1.26099};
const double aUV[4] = {0.268257  ,0.266635   ,0.266183   ,0.236486};
const double aIR[4] = {0.00047342,0.000848595,0.000846994,0.0022611};
const double rho[4] = {0.03549   ,0.03513    ,0.03481    ,0.03449}; // mol/cm3
double density[4];
double kbTrhokT[4];
void init()
{
for (int i=0;i<4;i++)
  {
    density[i]  = rho[i] *aw;  // g/cm3
    kbTrhokT[i] = kb*T[i]*kT*density[i]*density[i]; 
  }
}
double lambdatoe(double lambda)
{
  // input   photon wavelength in nm 
  // return  energy in eV
 double   E  = (h*c)/(lambda*1.e-9);
 return E;
}
double etolambda(double E)
{
  // input  photon energy in eV
  // return   wavelength in nm 
 double   lambda  = (h*c)/(E*1.e-9);
 return lambda;
}
double sellmeier_LAr(double * x,double *p)
{
  double la0  = p[0];
  double laUV = p[1];
  double laIR = p[2];
  double lambda =x[0]*1e-7;   // convert from nm to cm
  double nsquare = la0
    +(laUV*lambda*lambda)/(lambda*lambda-lambdaUV*lambdaUV)
    +(laIR*lambda*lambda)/(lambda*lambda-lambdaIR*lambdaIR);  
  return  sqrt(nsquare); 
}
double sellmeier_LAr(double  lambda,int index)
{
  lambda=lambda*1e-7;
  double nsquare = a0[index]
    +(aUV[index]*lambda*lambda)/(lambda*lambda-lambdaUV*lambdaUV)
    +(aIR[index]*lambda*lambda)/(lambda*lambda-lambdaIR*lambdaIR);  
  return sqrt(nsquare); 
}
double sellmeierpe_LAr(double * x,double *p)
{
  double la0  = p[0];
  double laUV = p[1];
  double laIR = p[2];
  double pe  = x[0];
  double lambda =  etolambda(pe);
  lambda =lambda*1e-7;   // convert from nm to cm
  double nsquare = la0
    +(laUV*lambda*lambda)/(lambda*lambda-lambdaUV*lambdaUV)
    +(laIR*lambda*lambda)/(lambda*lambda-lambdaIR*lambdaIR);  
  double nord =  sqrt(nsquare); 
  return nord;
}
double sellmeierpe_LAr(double x,int index)
{
  double lambda =  etolambda(x);
  lambda =lambda*1e-7;   // convert from nm to cm
  double nsquare = a0[index]
    +(aUV[index]*lambda*lambda)/(lambda*lambda-lambdaUV*lambdaUV)
    +(aIR[index]*lambda*lambda)/(lambda*lambda-lambdaIR*lambdaIR);  
  double nord =  sqrt(nsquare); 
  return nord;
}

// emin and emax in nm
void sellmeierLAr(double emin=110, double emax=700)
{
  const double minlambda = 110;
  const double maxlambda = 700;
  if (emin<minlambda||emax>maxlambda)
    {
      cout <<" variables out of range: "<< minlambda<<" - "<< maxlambda<<endl;
      return;
    }
  TCanvas *canvas2 = new TCanvas("canvas2", "refractive indices", 200, 10, 1000, 800);
  TF1 *sm0  = new TF1("sm0", sellmeier_LAr,emin,emax,3);
  sm0-> SetTitle("T=83.81 K");
  sm0->SetParameters(a0[0],aUV[0],aIR[0]);
  sm0->SetLineWidth(2); 
  sm0->SetLineColor(2);
  sm0->GetXaxis()->SetTitle("#lambda [nm]");
  sm0->GetYaxis()->SetTitle("Index of Refraction");
  sm0->Draw();
  TF1 *sm1  = new TF1("sm1", sellmeier_LAr,emin,emax,3);
  sm1-> SetTitle("T=86. K");
  sm1->SetParameters(a0[1],aUV[1],aIR[1]);
  sm1->SetLineWidth(2); 
  sm1->SetLineColor(3);
  sm1->Draw("SAME");
  TF1 *sm2  = new TF1("sm2", sellmeier_LAr,emin,emax,3);
  sm2-> SetTitle("T=88. K");
  sm2->SetParameters(a0[2],aUV[2],aIR[2]);
  sm2->SetLineWidth(2); 
  sm2->SetLineColor(4);
  sm2->Draw("SAME");
  TF1 *sm3  = new TF1("sm3", sellmeier_LAr,emin,emax,3);
  sm3-> SetTitle("T=90. K");
  sm3->SetParameters(a0[3],aUV[3],aIR[3]);
  sm3->SetLineWidth(2); 
  sm3->SetLineColor(8);
  sm3->Draw("SAME");
  //
  TGraph *ve = new TGraph(nr,Lambda,R);
  ve-> SetTitle("Sinnock et al");
  ve->SetLineColor(2);
  ve->SetMarkerColor(2);
  ve->SetMarkerStyle(22);
  ve->SetMarkerSize(2);
  ve->Draw("PL");
  TLegend *leg = canvas2->BuildLegend(.7, .65, 0.85, .85);
  leg->Draw();
}

void sellmeierpeLAr()
{
  TF1 *sm  = new TF1("npe", sellmeierpe_LAr,1.9,11.3,3);
  sm->SetParameters(a0[0],aUV[0],aIR[0]);
  sm->SetLineWidth(2); 
  sm->SetLineColor(2);
  sm->GetXaxis()->SetTitle("photon energy [eV]");
  sm->GetYaxis()->SetTitle("Index of Refraction");
  sm->Draw();
}

double lrayleigh(Double_t *x,Double_t *p)
{
  double lambda = x[0];
  int index=int(p[0]);
  double n  = sellmeier_LAr(lambda,index);
  lambda = lambda*1.e-7;    // convert from nm to cm
  double l  = p3times16/(6.*lambda*lambda*lambda*lambda);
  double br = kbTrhokT[index] *(((n*n-1)*(n*n+2))/(3.*density[index])*(((n*n-1)*(n*n+2))/(3.*density[index])));
  l=l*br;
  l =1./l;
  return l;
}
double lrayleigh(double  lambda,int index)
{
  double n  = sellmeier_LAr(lambda,index);
  lambda = lambda*1.e-7;    // convert from nm to cm
  double l  = p3times16/(6.*lambda*lambda*lambda*lambda);
  double br = kbTrhokT[index] *(((n*n-1)*(n*n+2))/(3.*density[index])*(((n*n-1)*(n*n+2))/(3.*density[index])));
  l=l*br;
  l =1./l;
  return l;
}
// emin and emax in nm
// nsteps number of steps
void rayleigh(double emin=110, double emax=400)
{
  const double minlambda = 110;
  const double maxlambda = 400;
  if (emin<minlambda||emax>maxlambda)
    {
      cout <<" variables out of range: "<< minlambda<<" - "<< maxlambda<<endl;
      return;
    }
  TCanvas *canvas = new TCanvas("canvas", "rayleigh scattering length", 200, 10, 1000, 800);
  TF1 *vp0  = new TF1("vp0",lrayleigh,emin,emax,1);
  vp0->SetParameters(0.,0.);
  vp0->SetLineWidth(2); 
  vp0->SetLineColor(2);
  vp0->SetTitle("T=83.81K");
  vp0->GetXaxis()->SetTitle("#lambda [nm]");
  vp0->GetYaxis()->SetTitle("Rayleigh scattering length [cm]");
  vp0->Draw();
  TF1 *vp1  = new TF1("vp1",lrayleigh,emin,emax,1);
  vp1->SetParameters(1.,0.);
  vp1->SetLineWidth(2); 
  vp1->SetLineColor(3);
  vp1->SetTitle("T=86.K");
  vp1->Draw("SAME");
  TF1 *vp2  = new TF1("vp2",lrayleigh,emin,emax,1);
  vp2->SetParameters(2.,0.);
  vp2->SetLineWidth(2); 
  vp2->SetLineColor(4);
  vp2->SetTitle("T=88.K");
  vp2->Draw("SAME");
  TF1 *vp3  = new TF1("vp3",lrayleigh,emin,emax,1);
  vp3->SetParameters(3.,0.);
  vp3->SetLineWidth(2); 
  vp3->SetLineColor(8);
  vp3->SetTitle("T=90.K");
  vp3->Draw("SAME");
  TLegend *leg = canvas->BuildLegend(.15, .75, 0.45, .85);
  leg->Draw();
}
//----------------------------------------------------------------------
// function prints out the rayleigh length in the geant 4 gdml format
// emin and emax in nm
// nsteps number of steps
// index of Temperature Array
//----------------------------------------------------------------------
void rayleightable(double emin=110, double emax=700, int nsteps=100,int index=0)
{
  const double minlambda = 110;
  const double maxlambda = 700;
  if (emin<minlambda||emax>maxlambda)
    {
      cout <<" variables out of range: "<< minlambda<<" - "<< maxlambda<<endl;
      return;
    }
  double stepsize = (emax-emin)/nsteps; 
  double pe       = emax;
  double n        = lrayleigh(pe,index);
  double photone  = lambdatoe(pe);
  cout<<"     <matrix name=\"RAYLEIGH\" coldim=\"2\" values=\"" << photone <<"*eV "<<n<<"*cm"<<endl;
  for (int i=1;i<nsteps-1;i++)
    {
      pe = emax-i*stepsize;
      n = lrayleigh(pe,index);
      photone = lambdatoe(pe);
      cout << photone <<"*eV "<<n<<"*cm"<<endl;
    }
  pe = emax-(nsteps-1)*stepsize;
  n = lrayleigh(pe,index);
  photone = lambdatoe(pe);
  cout << photone <<"*eV "<<n<<"*cm\"/>"<<endl;
}
//----------------------------------------------------------------------
// function prints out the refraction index in the geant 4 gdml format
// emin and emax in nm
// nsteps number of steps
// index of Temperature Array
// --------------------------------------------------------------------
void rindextable(double emin=110, double emax=700, int nsteps=100,int index=0)
{
  const double minlambda = 110;
  const double maxlambda = 700;
  if (emin<minlambda||emax>maxlambda)
    {
      cout <<" variables out of range: "<< minlambda<<" - "<< maxlambda<<endl;
      return;
    }
  double stepsize= (emax-emin)/nsteps;
  double pe = emax;
  double n = sellmeier_LAr(pe,index);
  double photone = lambdatoe(pe);
  cout<<"     <matrix name=\"RINDEX\" coldim=\"2\" values=\"" << photone <<"*eV "<<n<<endl;
  for (int i=1;i<nsteps-1;i++)
    {
      pe = emax-i*stepsize;
      n = sellmeier_LAr(pe,index);
      photone = lambdatoe(pe);
      cout <<photone <<"*eV "<<n<<endl;
    }
  pe = emax-(nsteps-1)*stepsize;
  n = sellmeier_LAr(pe,index);
  photone = lambdatoe(pe);
  cout << photone <<"*eV "<<n<<"\"/>"<<endl;
}
