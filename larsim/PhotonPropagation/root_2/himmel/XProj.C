void XProj()
{
  for (int integer =50; integer <51 ; integer= integer+5){ 
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
    //XPROJECTION
    //c1->SetCanvasSize(1700, 1200);
    //c1->SetWindowSize(1700, 1200);

    c1->SetCanvasSize(1000, 1200);
    c1->SetWindowSize(1000, 1200);
    //gPad->SetRightMargin(2.5);
    
    //gStyle->SetTitleOffset(.9,"y");

    c1->Divide(1,3);
    // Load the hist file we produced
    //TFile *hyb = new TFile("std_wFC_hist.root");
    TFile *hybsmall = new TFile("hybrid_onlydsmall_49.root");
    //TFile *hyb = new TFile("hybrid_with_FC_49_imp.root");
    //TFile *hyb = new TFile("hybrid_with_FC_49.root");
    //TFile *hyb = new TFile("hybrid_without_FC_49.root");
    //TH2D* hist1 = (TH2D*) hyb->Get("analyze/projZ206");
    TH2D* hist1 = (TH2D*) hybsmall->Get(TString::Format("analyze/projY%d", integer).Data());
    //TH2D* hist1 = (TH2D*) hyb->Get("analyze/ZProjection");
    
    TFile *hyb = new TFile("hybrid_with_FC_49.root");
    //TFile *std = new TFile("std_withoutFC_49.root"); 
    //TH2D* hist2 = (TH2D*) std->Get("analyze/projZ206");
    TH2D* hist2 = (TH2D*) hyb->Get(TString::Format("analyze/projY%d", integer).Data());
    //TH2D* hist2 = (TH2D*) std->Get("analyze/ZProjection");
    
    //TFile *std = new TFile("std_withoutFC_49.root");
    //TFile *std = new TFile("std_withoutFC_49.root"); 
    //TH2D* hist2 = (TH2D*) std->Get("analyze/projZ206");
    //TH2D* hist3 = (TH2D*) std->Get(TString::Format("analyze/projX%d", integer).Data());
   
    TH2D *hist3 = (TH2D*) hist1->Clone("hist3");
    
    c1->cd(1);
    hist1->GetZaxis()->SetLabelSize(0.05);
    hist1->GetYaxis()->SetLabelSize(0.05);
    hist1->GetXaxis()->SetLabelSize(0.05);
    hist1->GetYaxis()->SetRange(130,250);
    //hist1->GetXaxis()->SetRange(140,250);
    hist1->Draw("colz");
    hist1->SetTitle(TString::Format("Hybrid Library only small distances replaced projY%d", integer).Data());
    //hist1->SetTitle("Hybrid Library");
    //hist1->SetTitleSize(.09, "");
    hist1->SetMaximum(0.01);
    hist1->SetMinimum(-0.00001);
      
    c1->cd(2);
    hist2->GetZaxis()->SetLabelSize(0.05);
    hist2->GetYaxis()->SetLabelSize(0.05);
    hist2->GetXaxis()->SetLabelSize(0.05);
    hist2->GetYaxis()->SetRange(130,250);
    //hist2->GetYaxis()->SetRange(30,80);
    //hist2->GetXaxis()->SetRange(140,250);
    hist2->SetTitle(TString::Format("Hybrid Library projY%d", integer).Data());
    //hist2->SetTitle("Standard Library");
    c1->SetRightMargin(1.4);
    hist2->Draw("colz");
    c1->SetRightMargin(1.4);
    hist2->SetMaximum(0.01);
    hist2->SetMinimum(-0.00001);
    

    c1->cd(3);
    gStyle->SetNumberContours(101);
    hist3->GetZaxis()->SetLabelSize(0.05);
    hist3->GetYaxis()->SetLabelSize(0.05);
    hist3->GetXaxis()->SetLabelSize(0.05);
    hist3->GetYaxis()->SetRange(130,250);
    //hist3->GetYaxis()->SetRange(30,80);
    //hist3->GetXaxis()->SetRange(140,250);
    hist3->SetTitle(TString::Format("Hybrid only small distances replaced/ Hybrid Library projY%d OptDet #50", integer).Data());
    //hist3->SetTitle("Ratio Hyb/Std Zproj_206 OptDet #50");
    //hist3->SetTitle("Ratio Hybrid/Std Library");
    hist3->Divide(hist2);
    c1->SetRightMargin(1.4);
    hist3->Draw("colz");
    hist3->SetMaximum(2);
    hist3->SetMinimum(0.1);
    

    //gPad->Print(TString::Format("plots/projX48.png"));
    //c1->SaveAs("plots_50/50_XProjection.pdf");
    //hist3->GetXaxis()->SetTickLength(0.02);
    //hist3->GetYaxis()->SetTickLength(0.02);
    c1->SetRightMargin(1.4);
    c1->Update();
    c1->SaveAs(TString::Format("plot/projY%d #50.png", integer).Data());
    //c1->SaveAs("plot/Zproj_206 #50.png");
    delete c1;
    }
}
    
    // Creat the histogram
    //TH1D *hYScan0 = new TH1D("hYScan0", TString::Format("Detector %i;y position of #mu (cm);A.U.",iPD), nevents, axisarray);

 
