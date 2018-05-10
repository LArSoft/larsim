void test()
{
  for (int integer =0; integer <19 ; integer++){ 
  //Double_t blue[]    = {0, 0, 1};
  //Double_t red[]     = {1, 0, 0};
  //Double_t white[]   = {1, 1, 1};

  //Double_t Red[]    = {blue[0], white[0], red[0]};
  //Double_t Green[]  = {blue[1], white[1], red[1]};
    //Double_t Blue[]   = {blue[2], white[2], red[2]};
    //Double_t Length[] = {   0.00,  0.5,  1.00};

    //TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,101);
    //gStyle->SetNumberContours(101); 
    TCanvas *c1=new TCanvas("c1", "c1");
    gStyle->SetOptStat(0);
    c1->SetCanvasSize(900, 500);
    c1->SetWindowSize(900, 500);
    
    TFile* hyb = new TFile("fit.root");
    //hyb->cd("opdet_0");
    TTree* tree1 = (TTree*) hyb->Get(TString::Format("opdet_%d/tr", integer).Data());
    
    TFile* std = new TFile("full.root");
    //std->cd("opdet_0");
    TTree* tree2 = (TTree*) std->Get(TString::Format("opdet_%d/tr", integer).Data());
    
    TH1D* hist1 = new TH1D("hist1", "hist1", 200, 0, 2000);
    tree1->Draw("dist>>hist1");
    TH1D* hist2 = new TH1D("hist2", "hist2 title", 200, 0, 2000);
    tree2->Draw("dist>>hist2");
    TH1D* hist3 = (TH1D*) hist1->Clone("hist3");
    hist3->Divide(hist2);
    //c1->SetLogy();
    //hist1->GetYaxis()->SetRangeUser(0,1000);
    
    hist3->GetYaxis()->SetLabelSize(0.04);
    hist3->GetXaxis()->SetLabelSize(0.04);
    hist3->GetXaxis()->SetTickLength(0.02);
    hist3->GetYaxis()->SetTickLength(0.02);
    hist3->GetXaxis()->SetTitle("distance [cm]");
    hist3->Draw("hist");
    //hist1->Draw("same");
    hist3->SetTitle(TString::Format("Ratio of exception points/total amount of points as a function of the distance for opdet_%d",integer).Data());
    c1->Modified();
    c1->SaveAs(TString::Format("plot/EP(distance)_imp_opdet_%d.png",integer).Data());
    delete c1;
  }
    
}
    
    

 
