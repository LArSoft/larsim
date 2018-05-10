void Tutorial2()
{
  //for (int integer =2; integer <49 ; integer++){ 
    Double_t blue[]    = {0, 0, 1};
    Double_t red[]     = {1, 0, 0};
    Double_t white[]   = {1, 1, 1};

    Double_t Red[]    = {blue[0], white[0], red[0]};
    Double_t Green[]  = {blue[1], white[1], red[1]};
    Double_t Blue[]   = {blue[2], white[2], red[2]};
    Double_t Length[] = {   0.00,  0.5,  1.00};

    TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,101);
    //gStyle->SetNumberContours(101); 
    TCanvas *c1=new TCanvas("c1", "c1");
    gStyle->SetOptStat(0);
    c1->SetCanvasSize(900, 500);
    c1->SetWindowSize(900, 500);
    //c1->Divide(1,3);
    // Load the hist file we produced
    //TFile *hyb = new TFile("std_wFC_hist.root");
    TFile *hyb = new TFile("hybrid_without_FC_49.root");
    //TFile *hyb = new TFile("hyb_withFC_imp.root");
    TH2D* hist1 = (TH2D*) hyb->Get("analyze/projZ203");
    //TH2D* hist1 = (TH2D*) hyb->Get(TString::Format("analyze/projX%d", integer).Data());
    //TH2D* hist1 = (TH2D*) hyb->Get("analyze/ZProjection");
    
    //TFile *std = new TFile("std_wFC_hist.root");
    //TFile *std = new TFile("std_wFC_hist.root"); 
    TFile *std = new TFile("std_withoutFC_49.root");
    //TH2D* hist2 = (TH2D*) std->Get("analyze/projZ230");
    //TH2D* hist2 = (TH2D*) std->Get(TString::Format("analyze/projX%d", integer).Data());
    //TH2D* hist2 = (TH2D*) std->Get("analyze/XProjection");
    TH2D* hist2 = (TH2D*) std->Get("analyze/projZ203");
    TH2D *hist3 = (TH2D*) hist1->Clone("hist3");
    
    /* c1->cd(1);
    hist1->GetZaxis()->SetLabelSize(0.07);
    hist1->GetYaxis()->SetLabelSize(0.07);
    hist1->GetXaxis()->SetLabelSize(0.07);
    hist1->Draw("colz");
    //hist1->SetTitle("Hybrid Library Zproj_230 with FC");
    hist1->SetTitle("Hybrid Library XProjection with FC");
    //hist1->SetTitle("Hybrid Library");
    //hist1->SetTitleSize(.09, "");
    hist1->SetMinimum(-0.0000001);
    //hist1->SetMaximum(0.15);

    // if(integer>30){
    //hist1->SetMaximum(0.16);
    //}
    //else hist1->SetMaximum(0.1);
    //hist1->SetMinimum(-0.00001);    
    c1->cd(2);
    hist2->GetZaxis()->SetLabelSize(0.07);
    hist2->GetYaxis()->SetLabelSize(0.07);
    hist2->GetXaxis()->SetLabelSize(0.07);
    //hist2->SetTitle("Standard Library Zproj_230 with FC");
    hist2->SetTitle("Hybrid Library XProjection with FC improved");
    //hist2->SetTitle("Standard Library");
    hist2->Draw("colz");
    //if(integer>30){
    //hist2->SetMaximum(40);
    //}
    //else hist2->SetMaximum(0.1);
    hist2->SetMinimum(-0.0000001);
    //hist2->SetMaximum(0.15);
    */
    //c1->cd(3);
    gStyle->SetNumberContours(101);
    hist3->GetZaxis()->SetLabelSize(0.06);
    hist3->GetYaxis()->SetLabelSize(0.06);
    hist3->GetXaxis()->SetLabelSize(0.06);
    //hist3->SetTitle("Ratio Hybrid  Library not/imporved with FC XProjection");
    hist3->SetTitle("Ratio Hybrid/Std Library without FC");
    hist3->Divide(hist2);
    hist3->Draw("colz");
    hist3->SetMaximum(1.5);
    hist3->SetMinimum(0.5);
    //gPad->Print(TString::Format("plots/projX48.png"));
    //c1->SaveAs("plots_50/50_XProjection.pdf");
    //c1->SaveAs(TString::Format("Xplots_all/all_projX%d.png", integer).Data());
    c1->SaveAs("plot/ZProj_203_Opdet50.png");
    delete c1;
    //}
}
    
    // Creat the histogram
    //TH1D *hYScan0 = new TH1D("hYScan0", TString::Format("Detector %i;y position of #mu (cm);A.U.",iPD), nevents, axisarray);

 
