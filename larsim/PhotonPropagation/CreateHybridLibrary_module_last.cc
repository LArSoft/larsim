////////////////////////////////////////////////////////////////////////
// Chris Backhouse, UCL, Nov 2017
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"

#include <iostream>

#include "larsim/PhotonPropagation/PhotonVisibilityService.h"

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TStyle.h"

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
      float vis, dist, costheta;
      tr_full->Branch("dist", &dist);
      tr_full->Branch("costheta", &costheta);
      tr_full->Branch("vis", &vis);

      const geo::OpDetGeo& opdet = geom->OpDetGeoFromOpDet(opdetIdx);
      double xyzopdet[3];
      opdet.GetCenter(xyzopdet);
      const TVector3 opdetvec(xyzopdet);
      //std::cout << opdetIdx << " / " << geom->NOpDets() << "position center (x, y, z)"<< opdetvec[0]<<"  "<<opdetvec[1]<<"  "<<opdetvec[2]<< std::endl;
      struct Visibility{
        Visibility(int vx, float d, float c, float v) : vox(vx), dist(d), costheta(c), vis(v) {}
        int vox;
        float dist;
        float costheta;
        float vis;
      };
      /*TCanvas *c1=new TCanvas("c1", "c1");
      c1->SetCanvasSize(1500, 1500);
      c1->SetWindowSize(600, 600);
      c1->Divide(1,2);*/
      TGraph g;

      std::vector<Visibility> viss;
      viss.reserve(voxdef.GetNVoxels());

      for(int voxIdx = 0; voxIdx < voxdef.GetNVoxels(); ++voxIdx){
	//std::cout << voxIdx << " / " << voxdef.GetNVoxels() << std::endl;
        const TVector3 voxvec = voxdef.GetPhotonVoxel(voxIdx).GetCenter();
        const double xyzvox[] = {voxvec.X(), voxvec.Y(), voxvec.Z()};
	
        dist = opdet.DistanceToPoint(xyzvox);
        vis = pvs->GetVisibility(xyzvox, opdetIdx);
        g.SetPoint(g.GetN(), dist, vis*dist*dist); //pinta el grafico

        tr_full->Fill();

        viss.emplace_back(voxIdx, dist, costheta, vis);
      } // end for voxIdx

      d_full->cd();
      tr_full->Write();
      delete tr_full;

      g.SetMarkerStyle(7);
      //c1->cd(1);
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
      int vox;
      tr_fit->Branch("vox", &vox);
      tr_fit->Branch("vis", &vis);
      tr_fit->Branch("dist", &dist);
      
      TGraph g2;
      //TCanvas *c2=new TCanvas();
      //TGraph g3; //B
      for(const Visibility& v: viss){
        // 2e-5 is the magic scaling factor to get back to integer photon
        // counts. TODO this will differ for new libraries, should work out a
        // way to communicate it or derive it.
        const double obs = v.vis / 2e-5;
        const double pred = fit->Eval(v.dist) / (v.dist*v.dist) / 2e-5; //calculated with parametrisation

        // Log-likelihood ratio for poisson statistics
        double chisq = 2*(pred-obs);
        if(obs) chisq += 2*obs*log(obs/pred);
	//g3.SetPoint(g3.GetN(), dist, chisq);
        h.Fill(chisq);
       
        if(chisq > 9){ //equivalent to more than 9 chisquare = 3 sigma    //maybe play around with this cutoff
	  g2.SetPoint(g2.GetN(), v.dist, obs*2e-5*v.dist*v.dist);
	  //g3.SetPoint(g3.GetN(), v.vox, chisq);
	  // if(2710000<v.vox && v.vox<2719000){
	  //std::cout << "Voxel ID # "<< v.vox  << "         Chisquare = "<< chisq <<std::endl;}
	  vox = v.vox;
          vis = v.vis;
          dist = v.dist;
          costheta = v.costheta;
          tr_fit->Fill();
        }
      }

      d_fit->cd();
      tr_fit->Write();
      delete tr_fit;
      g2.SetMarkerSize(1); 
      g2.SetLineWidth(3);
      g2.SetMarkerStyle(7);
      g2.SetMarkerColor(kRed);
      //gStyle->SetLabelSize(.04, "XY");
      //gStyle->SetTitleSize(.04,"XY");
      //c1->cd(1);
      g2.Draw("p same"); // p same
      gPad->Print(TString::Format("plots/vis_vs_dist_%d.png", opdetIdx).Data());
      //c1->SaveAs("plots/Bea_vis_dist.pdf");
      

      //g3.SetLineColor(2) 
      /*g3.SetMarkerSize(1); 
      g3.SetLineWidth(3);
      g3.SetMarkerStyle(7);
      g3.GetXaxis()->SetTitle("Voxels");
      g3.GetXaxis()->SetRangeUser(2715000,2719000);
      //g3.GetYaxis()->SetRangeUser(0,1000);
      g3.GetYaxis()->SetTitle("#chi^{2}");
      gStyle->SetLabelSize(.04, "XY");
      gStyle->SetTitleSize(.04,"XY");
      gPad->SetLogy();
      c1->cd(2);
      //gPad->SetLogy();
      g3.Draw("AP*");
      //gPad->Print(TString::Format("plots/chisq_vs_dist_%d.png", opdetIdx).Data());
      c1->SaveAs(TString::Format("plots/all_%d.png", opdetIdx).Data());*/
      //c2->Print("plots/Bea_chisq_dist.pdf");
      
      
      //g3.SetMarkerStyle(7);//B
      //g3.SetMarkerColor(kRed);
      //g3.Draw("p same"); // p same
      //c2.Print("plots.pdf");
      
      h.GetXaxis()->SetTitle("#chi^{2}");
      h.GetYaxis()->SetTitle("Counts");
      h.Draw("hist");
      std::cout <<"Integral(90,-1)/Integral(all) = "<< h.Integral(90, -1) / h.Integral(0, -1) << std::endl;
      std::cout <<"*****************************" << std::endl;
      totExceptions += h.Integral(90, -1);
      totPts += h.Integral(0, -1);
      gPad->Print(TString::Format("plots/chisq_opdet_%d.eps", opdetIdx).Data());
            
      //delete c1;
      //delete c2;
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
