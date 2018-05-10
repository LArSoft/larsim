void Ratio()
{
  //for (int integer =0; integer <1 ; integer++){
    TCanvas *c1=new TCanvas("c1", "c1");
    gStyle->SetOptStat(0);
    //c1->SetCanvasSize(1700, 1200);
    //c1->SetWindowSize(1700, 1200);
    
    TFile* hyb = new TFile("fit_imp.root");
    TH1D* hist1 = (TH1D*) hyb->Get("opdet_0/dist");
    //TH1D* hist1 = (TH1D*) hyb->Get(TString::Format("opdet_%d/dist", integer).Data());
    //TFile* std = new TFile("full_imp.root");
    //TH1D* hist2 = (TH1D*) std->Get(TString::Format("opdet_%d/dist", integer).Data());
    //TH1D* hist3 = (TH1D*) hist1->Clone("hist3");
    
    //hist1->SetTitle(TString::Format("Ratio of amount of exceptions/amount of visibilities (over distance) for OptDet%d", integer).Data());
    //hist3->SetTitle("Ratio Hyb/Std Zproj_206 OptDet #50");
    //hist3->SetTitle("Ratio Hybrid/Std Library");
    //hist3->Divide(hist2);
    hist1->Draw("hist");
    c1->Modified();
    //c1->SaveAs(TString::Format("plot/OpDet%d_.png", integer).Data());
    c1->SaveAs("plot/OpDet0.png");
    delete c1;
    //}
}
    
   

 
