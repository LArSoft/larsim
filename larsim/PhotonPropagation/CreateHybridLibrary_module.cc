
////////////////////////////////////////////////////////////////////////
// Chris Backhouse, UCL, Nov 2017
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include <iostream>

#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMultiGraph.h"

#define PI 3.14159265

namespace phot
{
	class CreateHybridLibrary: public art::EDAnalyzer
	{
		public:
			explicit CreateHybridLibrary(const fhicl::ParameterSet& p);

			// Plugins should not be copied or assigned.
			CreateHybridLibrary(const CreateHybridLibrary&) = delete;
			CreateHybridLibrary(CreateHybridLibrary&&) = delete;
			CreateHybridLibrary& operator=(const CreateHybridLibrary&) = delete;
			CreateHybridLibrary& operator=(CreateHybridLibrary&&) = delete;

			void analyze(const art::Event& e) override;
	};

	//--------------------------------------------------------------------
	CreateHybridLibrary::CreateHybridLibrary(const fhicl::ParameterSet& p)
		: EDAnalyzer(p)
	{ 
		// std::cout<<"BEA::Photon library analyzer constructor "<<std::endl;
		art::ServiceHandle<geo::Geometry> geom;

		art::ServiceHandle<phot::PhotonVisibilityService> pvs;
		sim::PhotonVoxelDef voxdef = pvs->GetVoxelDef();

		TFile* fout_full = new TFile("full.root", "RECREATE");
		TFile* fout_fit = new TFile("fit.root", "RECREATE");

		std::cout << voxdef.GetNVoxels() << " voxels for each of " << geom->NOpDets() << " OpDets" << std::endl;
		std::cout << std::endl;

		long totExceptions = 0;
		long totPts = 0;

		for(unsigned int opdetIdx =0; opdetIdx < geom->NOpDets(); ++opdetIdx){
		  std::cout << opdetIdx << " / " << geom->NOpDets() << std::endl;

			TDirectory* d_full = fout_full->mkdir(TString::Format("opdet_%d", opdetIdx).Data());
			TDirectory* d_fit = fout_fit->mkdir(TString::Format("opdet_%d", opdetIdx).Data());

			d_full->cd();
			TTree* tr_full = new TTree("tr", "tr");
			int vox, taken;
			float dist, vis, psi, theta;
			tr_full->Branch("vox", &vox);
			tr_full->Branch("dist", &dist);
			tr_full->Branch("vis", &vis);
			tr_full->Branch("taken", &taken);
			tr_full->Branch("psi", &psi);
			tr_full->Branch("theta", &theta);

			const geo::OpDetGeo& opdet = geom->OpDetGeoFromOpDet(opdetIdx);
			double xyzopdet[3];
			opdet.GetCenter(xyzopdet);
			//std::cout << "OpDet position " << opdetIdx << " x =  " <<xyzopdet[0] << " y =  " <<xyzopdet[1] << " z =  " <<xyzopdet[2]<< std::endl;
			const TVector3 opdetvec(xyzopdet);

			struct Visibility{
			  Visibility(int vx, int t, float d, float v, float p, float th) : vox(vx), taken(t), dist(d), vis(v), psi(p), theta(th) {}
				int vox;
				int taken;
			        float dist;
				float vis;
			        float psi;
			        float theta;
			};
			TCanvas *c1=new TCanvas("c1", "c1");
			//c1->SetCanvasSize(1500, 1500);
			c1->SetWindowSize(600, 600);
			//c1->Divide(1,2);
			TGraph g;
			TGraph g2;

			std::vector<Visibility> viss;
			viss.reserve(voxdef.GetNVoxels());

			for(int voxIdx = 0; voxIdx < voxdef.GetNVoxels(); ++voxIdx){
			     
			  //std::cout << voxIdx << " / " << voxdef.GetNVoxels() << std::endl;
				const TVector3 voxvec = voxdef.GetPhotonVoxel(voxIdx).GetCenter();
				const double xyzvox[] = {voxvec.X(), voxvec.Y(), voxvec.Z()};
				//SP const double fc_y = 600; //624cm is below the center of the first voxel outside the field cage 
				//SP const double fc_z = 1394;
				//SP const double fc_x = 350;
				//SP taken = 0;
				//DP does not need taken because all voxels are inside the Field Cage
				taken = 1;
				dist = opdet.DistanceToPoint(xyzvox);
				vis = pvs->GetVisibility(xyzvox, opdetIdx);
				//SP psi = 100;
				//SP theta = 200;
				psi = 0;
				theta = 0;

				/*SP if((xyzvox[0] - xyzopdet[0])<0){
				  psi = atan((xyzvox[1] - xyzopdet[1])/(-xyzvox[0] +xyzopdet[0]));
				  }*/
				//Test angles std::cout << "                        "<<std::endl;
				//std::cout << " Theta = " << theta * 180.0/PI<< " Psi = " << psi * 180.0/PI << " z_vox = "<< xyzvox[2] <<" z_opdet = "<< xyzopdet[2]<<" dist = "<< dist << " x_vox = "<< xyzvox[0] <<" x_opdet = "<< xyzopdet[0]<< "y_vox = "<< xyzvox[1] <<" y_opdet = "<< xyzopdet[1]<<std::endl;
				//if (dist>125){
				//SP if(xyzvox[0]<fc_x && xyzvox[0]>-fc_x){
				//if (xyzvox[1]<fc_y && xyzvox[1]>-fc_y){
				//if (xyzvox[2]> -9 && xyzvox[2]<fc_z){
				  g.SetPoint(g.GetN(), dist, vis*dist*dist); //pinta el grafico
				  //SP taken = 1;
				  //SP psi = atan((xyzvox[1] - xyzopdet[1])/(xyzvox[0] - xyzopdet[0]))* 180.0/PI ; // (-PI/2, PI/2)
				  //SP theta = acos((xyzvox[2] - xyzopdet[2])/dist) * 180.0/PI;  // (0 (beam direction := delta z = distance), PI (beam direction := delta z = -distance))          
				  /*SP
				  if((xyzvox[0] - xyzopdet[0])<0){
				  psi = atan((xyzvox[1] - xyzopdet[1])/(-xyzvox[0] +xyzopdet[0]));
				  }*/
			
				  //}
				  //}
				  //}
      				tr_full->Fill();
				viss.emplace_back(voxIdx, taken, dist, vis, psi, theta);
			} // end for voxIdx

			d_full->cd();
			tr_full->Write();
			delete tr_full;

			g.SetMarkerStyle(7);
			c1->cd();
			g.Draw("ap");
			g.GetXaxis()->SetTitle("Distance (cm)");
			g.GetYaxis()->SetTitle("Visibility #times r^{2}");
			g.Fit("expo");
			TF1* fit = g.GetFunction("expo");

			d_fit->cd();
			TVectorD fitres(2);
			fitres[0] = fit->GetParameter(0);
			fitres[1] = fit->GetParameter(1);
			fitres.Write("fit");

			gPad->SetLogy();

			TH1F h("", "", 200, 0, 20); //20 the last digit

			d_fit->cd();
			TTree* tr_fit = new TTree("tr", "tr");
			tr_fit->Branch("vox", &vox);
			tr_fit->Branch("dist", &dist);
			tr_fit->Branch("vis", &vis);
			tr_fit->Branch("taken", &taken);
			tr_fit->Branch("psi", &psi);
			tr_fit->Branch("theta", &theta);
			

			for(const Visibility& v: viss){
				// 2e-5 is the magic scaling factor to get back to integer photon
				// counts. TODO this will differ for new libraries, should work out a
				// way to communicate it or derive it.
				//const double obs = v.vis / 2e-5;
				//const double pred = fit->Eval(v.dist) / (v.dist*v.dist) / 2e-5; //calculated with parametrisation
				
				const double obs = v.vis / 1e-7;
				const double pred = fit->Eval(v.dist) / (v.dist*v.dist) / 1e-7; //calculated with parametrisation
				//std::cout << "observed = "<<obs<<" predicted (by parametrization) = "<<pred <<std::endl;
				//Minimal amount of detected photons is set to 50, bc of Landau dustribution
				//Those voxels with detected photons < 50 where set to 0

				// Log-likelihood ratio for poisson statistics
				double chisq = 2*(pred-obs);
				if(obs) chisq += 2*obs*log(obs/pred); //BEA
				/*if(obs==0 || obs==1){
				  chisq = 2; 
				  }*/
				vox = v.vox;
				dist = v.dist;
				//SP vis = pred *2e-5;
				psi = v.psi;
				theta = v.theta;
				vis = pred *1e-7;  
				

				if (v.taken==1){
				  h.Fill(chisq);	
				}


				//if(chisq > 9 && dist < 140){ //equivalent to more than 9 chisquare = 3 sigma    //maybe play around with this cutoff
				if(chisq > 9){ //equivalent to more than 9 chisquare = 3 sigma    //maybe play around with this cutoff
					g2.SetPoint(g2.GetN(), v.dist, v.vis*v.dist*v.dist);
					//SP vis = obs *2e-5;
					vis = obs *1e-7;
    					tr_fit->Fill();
					//std::cout << "V in lib = " << vis << " at vox_y = "<< xyzvox[1] << " at vox_z = "<< xyzvox[2]<<" voxIdx " << voxIdx <<std::endl;
				}

			}

			d_fit->cd();
			tr_fit->Write();
			delete tr_fit;
			g2.SetMarkerSize(1); 
			g2.SetLineWidth(3);
			g2.SetMarkerStyle(7);
			g2.SetMarkerColor(kRed);
			gStyle->SetLabelSize(.04, "XY");
			gStyle->SetTitleSize(.04,"XY");
			c1->cd();
			g2.Draw("p same"); 
			gPad->Print(TString::Format("plots/vis_vs_dist_%d.png", opdetIdx).Data());
			c1->SaveAs(TString::Format("plots/Chris_vis_vs_dist_%d.png", opdetIdx).Data());


			h.GetXaxis()->SetTitle("#chi^{2}");
			h.GetYaxis()->SetTitle("Counts");
			gStyle->SetLabelSize(.04, "XY");
			gStyle->SetTitleSize(.04,"XY");
			h.Draw("hist");
			std::cout <<"Integral(90,-1)/Integral(all) = "<< h.Integral(90, -1) / h.Integral(0, -1) << std::endl;
			std::cout <<"*****************************" << std::endl;
			totExceptions += h.Integral(90, -1);
			totPts += h.Integral(0, -1);
			gPad->Print(TString::Format("plots/chisq_opdet_%d.eps", opdetIdx).Data());

			delete c1;
			//delete c2;
			//delete c3;
			//delete c4;
			//std::cout << "Plots for voxels with coordinates (x, y, z) " << x_vox <<"  y  " << z_vox << std::endl;
			} // end for opdetIdx

			std::cout << totExceptions << " exceptions from " << totPts << " points = " << (100.*totExceptions)/totPts << "%" << std::endl;

			delete fout_full;
			delete fout_fit;
			exit(0); // We're done :)
			}

			//--------------------------------------------------------------------
			void CreateHybridLibrary::analyze(const art::Event&)
			{
			}

			DEFINE_ART_MODULE(CreateHybridLibrary)
		} // namespace
