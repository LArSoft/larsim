void test()
{
  //for (int integer =2; integer <49 ; integer++){ 
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
    c1->SetCanvasSize(900, 400);
    c1->SetWindowSize(900, 400);
    
    TFile *hyb = new TFile("std_wFC_hist.root");
    TH2D* hist1 = (TH2D*) hyb->Get("analyze/XProjection");
    
    hist1->GetZaxis()->SetLabelSize(0.06);
    hist1->GetYaxis()->SetLabelSize(0.06);
    hist1->GetXaxis()->SetLabelSize(0.06);
    hist1->Draw("colz");
    hist1->SetTitle("Photon Library with field cage XProjection");
    //hist1->SetTitle("Hybrid Library");
    //hist1->SetTitleSize(.09, "");
    hist1->SetMinimum(-0.0000001);
    //hist1->SetMaximum(0.15);

    //gStyle->SetNumberContours(101);
        
    c1->SaveAs("plot/XProjection_poster.png");
    delete c1;
    //}
}
    
    // Creat the histogram
    //TH1D *hYScan0 = new TH1D("hYScan0", TString::Format("Detector %i;y position of #mu (cm);A.U.",iPD), nevents, axisarray);

 
