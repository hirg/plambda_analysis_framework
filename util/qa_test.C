void fullep_flatness_test(TFile* f, Int_t);
void halfep_flatness_test(TFile* f, Int_t);
void resolution_test(TFile*, Double_t&, Int_t);
void protonv2_test(TFile*, Double_t ep_res, Int_t);
void proton_flatness_test(TFile* f, Int_t);
void lambda_flatness_test(TFile* f, Int_t);
void lambda_bin_flatness_test(TFile* f, Int_t);
void lambdav2_test(TFile* f, Double_t ep_res, Int_t);
void gamma_test(TFile* f, Double_t ep_res, Int_t);
void gamma_eff_test(TFile* f, Double_t ep_res, Int_t);
Double_t chi(Double_t res);
Double_t resEventPlane(double chi);

void qa_test(){
    std::string input_file;
    Int_t centrality;
    std::cout << "Centrality: " << std::endl;
    std::cin >> centrality;

    std::cout << "Input file: " << std::endl;
    std::cin >> input_file;

    TFile* f = new TFile(input_file.c_str(), "read");

    // Event plane flattening test
    fullep_flatness_test(f, centrality);   

    // Sub-event plane flattening test
    halfep_flatness_test(f, centrality);

    // Event plane resolution
    Double_t ep_res = 1.;
    resolution_test(f, ep_res, centrality);

    // Proton v2
    protonv2_test(f, ep_res, centrality);   

    // Lambda v2
    lambdav2_test(f, ep_res, centrality);   

    // Proton phi distribution
    proton_flatness_test(f, centrality);

    // Lambda phi distribution
    lambda_flatness_test(f, centrality);

    // Lambda phi distribution differential version
    lambda_bin_flatness_test(f, centrality);

    // Gamma
    gamma_test(f, ep_res, centrality);

    // Gamma with eff. correction
    gamma_eff_test(f, ep_res, centrality);
}

void fullep_flatness_test(TFile* f, Int_t centrality){
    TCanvas* c0 = new TCanvas("c0", "c0", 700, 500);
    TH1F* h = (TH1F*)f->Get("h1f_before_Flattened_FullEPPhi"); 
    Double_t no_bins = h->GetSize() - 2.;
    Double_t sum = h->GetSum();
    h->Scale(1. / (sum / no_bins));
    h->GetYaxis()->SetRangeUser(0, 1.7);
    h->SetLineWidth(2);
    h->SetLineColor(1);
    h->Draw();

    TH1F* h1 = (TH1F*)f->Get("h1f_Flattened_FullEPPhi");
    h1->Scale(1. / (h1->GetSum() / no_bins)); 
    h1->GetYaxis()->SetRangeUser(0, 1.7); 
    h1->SetLineWidth(2);
    h1->SetLineColor(2);
    h1->Draw("same");

    char plotname[50];
    sprintf(plotname, "plots/fullep_flatness_test_cent%d.png", centrality);
    gPad->SaveAs(plotname);
}

void halfep_flatness_test(TFile* f, Int_t centrality){
    TCanvas* c1 = new TCanvas("c1", "c1", 700, 500);
    TH1F* h = (TH1F*)f->Get("h1f_before_Flattened_EastEPPhi"); 
    Double_t no_bins = h->GetSize() - 2.;
    Double_t sum = h->GetSum();
    h->Scale(1. / (sum / no_bins));
    h->GetYaxis()->SetRangeUser(0, 1.7);
    h->SetLineWidth(2);
    h->SetLineColor(1);
    h->Draw();

    TH1F* h1 = (TH1F*)f->Get("h1f_Flattened_EastEPPhi");
    h1->Scale(1. / (h1->GetSum() / no_bins)); 
    h1->GetYaxis()->SetRangeUser(0, 1.7); 
    h1->SetLineWidth(2);
    h1->SetLineColor(2);
    h1->Draw("same");
    char plotname0[50];
    sprintf(plotname0, "plots/east_flatness_test_cent%d.png", centrality);
    gPad->SaveAs(plotname0);

    TCanvas* c2 = new TCanvas("c2", "c2", 700, 500);
    TH1F* h2 = (TH1F*)f->Get("h1f_before_Flattened_WestEPPhi"); 
    no_bins = h2->GetSize() - 2.;
    sum = h2->GetSum();
    h2->Scale(1. / (sum / no_bins));
    h2->GetYaxis()->SetRangeUser(0, 1.7);
    h2->SetLineWidth(2);
    h2->SetLineColor(1);
    h2->Draw();

    TH1F* h3 = (TH1F*)f->Get("h1f_Flattened_WestEPPhi");
    h3->Scale(1. / (h3->GetSum() / no_bins)); 
    h3->GetYaxis()->SetRangeUser(0, 1.7); 
    h3->SetLineWidth(2);
    h3->SetLineColor(2);
    h3->Draw("same");

    char plotname[50];
    sprintf(plotname, "plots/west_flatness_test_cent%d.png", centrality);
    gPad->SaveAs(plotname);
}

void resolution_test(TFile* f, Double_t& ep_res, Int_t centrality){
    TProfile* profile_ep_resolution =  (TProfile*)f->Get("profile_eventplane_resolution");
    Float_t ep_resolution = profile_ep_resolution->GetBinContent(1);
    Double_t x = chi(sqrt(ep_resolution / 100.));
    Double_t real_resolution = resEventPlane(x * sqrt(2.));
    std::cout << "Full event plane resolution for centrality " << centrality << " is: " << real_resolution << std::endl;
    ep_res = real_resolution;
}

void protonv2_test(TFile* f, Double_t ep_res, Int_t centrality){
    TCanvas* c3 = new TCanvas("c3", "c3", 700, 500); 
    TProfile* profile_protonv2 = (TProfile*)f->Get("profile_XPt_Yv2_ProtonWeighted");
    profile_protonv2->Scale(1. / ep_res / 100.);
    profile_protonv2->GetXaxis()->SetRangeUser(0, 2.4);
    profile_protonv2->GetXaxis()->SetTitle("Pt(GeV/c)");
    profile_protonv2->GetYaxis()->SetTitle("v2");
    profile_protonv2->Draw();
    gPad->SetTicks(1, 1);
    char plotname[50];
    sprintf(plotname, "plots/protonv2_alone_cent%d.png", centrality);
    gPad->SaveAs(plotname);

    // Nasim publication 30-80%
    Double_t x[10] = {0.345, 0.512, 0.703, 0.898, 1.096, 1.294, 1.493, 1.692, 1.892, 2.092};
    Double_t v2[10] = {0.0127, 0.0252, 0.049, 0.079, 0.113, 0.146, 0.176, 0.202, 0.225, 0.242};
    TGraph* gr = new TGraph(10, x, v2);
    gr->SetMarkerStyle(33);
    gr->Draw("same P");
    char plotname1[50];
    sprintf(plotname1, "plots/protonv2_comparison_cent%d.png", centrality);
    gPad->SaveAs(plotname1);
}

void lambdav2_test(TFile* f, Double_t ep_res, Int_t centrality){
    TCanvas* c3_1 = new TCanvas("c3_1", "c3_1", 700, 500); 
    TProfile* profile_lambdav2 = (TProfile*)f->Get("profile_XPt_Yv2_LambdaWeighted");
    profile_lambdav2->Scale(1. / ep_res / 100.);
    profile_lambdav2->GetXaxis()->SetRangeUser(0, 2.4);
    profile_lambdav2->GetXaxis()->SetTitle("Pt(GeV/c)");
    profile_lambdav2->GetYaxis()->SetTitle("v2");
    profile_lambdav2->Draw();
    gPad->SetTicks(1, 1);
    char plotname[50];
    sprintf(plotname, "plots/lambdav2_alone_cent%d.png", centrality);
    gPad->SaveAs(plotname);

    // Nasim publication
    const int size = 12;
    Double_t x[size] = {0.455, 0.559, 0.655, 0.759, 0.860, 0.955, 1.103, 1.297, 1.488, 1.680,  1.880, 2.095};
    Double_t v2[size] = {0.017, 0.026, 0.036, 0.050, 0.062, 0.078, 0.1, 0.131, 0.159, 0.185, 0.208, 0.227};
    TGraph* gr = new TGraph(size, x, v2);
    gr->SetMarkerStyle(33);
    gr->Draw("same P");
    char plotname1[50];
    sprintf(plotname1, "plots/lambdav2_comparison_cent%d.png", centrality);
    gPad->SaveAs(plotname1);
}

void proton_flatness_test(TFile* f, Int_t centrality){
    TCanvas* c4 = new TCanvas("c4", "c4", 700, 500);
    TH1F* h = (TH1F*)f->Get("h1f_before_Corrections_AlphaPhi");
    Double_t no_bins = h->GetSize() - 2.;
    Double_t sum = h->GetSum();
    h->Scale(1. / (sum / no_bins));
    h->GetYaxis()->SetRangeUser(0, 2.0);
    h->SetLineWidth(2);
    h->SetLineColor(1);
    h->Draw();
    h->Fit("pol0");

    TH1F* h1 = (TH1F*)f->Get("h1f_after_Corrections_AlphaPhi");
    no_bins = h1->GetSize() - 2.;
    sum = h1->GetSum();
    h1->Scale(1. / (sum / no_bins));
    h1->SetLineWidth(2);
    h1->SetLineColor(2);
    h1->Draw("same");
    h1->Fit("pol0");
    cout << "before mean:" << h->GetMean() << " after mean: " << h1->GetMean() << endl;
    char plotname[50];
    sprintf(plotname, "plots/proton_flatness_test_cent%d.png", centrality);
    gPad->SaveAs(plotname);
}

void lambda_flatness_test(TFile* f, Int_t centrality){
    TCanvas* c5 = new TCanvas("c5", "c5", 700, 500);
    TH1F* h = (TH1F*)f->Get("h1f_before_Corrections_BetaPhi");
    Double_t no_bins = h->GetSize() - 2.;
    Double_t sum = h->GetSum();
    h->Scale(1. / (sum / no_bins));
    h->GetYaxis()->SetRangeUser(0, 2.0);
    h->SetLineWidth(2);
    h->SetLineColor(1);
    h->Draw();
    h->Fit("pol0");

    TH1F* h1 = (TH1F*)f->Get("h1f_after_Corrections_BetaPhi");
    no_bins = h1->GetSize() - 2.;
    sum = h1->GetSum();
    h1->Scale(1. / (sum / no_bins));
    h1->SetLineWidth(2);
    h1->SetLineColor(2);
    h1->Draw("same");
    h1->Fit("pol0");
    cout << "before mean:" << h->GetMean() << " after mean: " << h1->GetMean() << endl;
    char plotname[50];
    sprintf(plotname, "plots/lambda_flatness_test_cent%d.png", centrality);
    gPad->SaveAs(plotname);
}

void lambda_bin_flatness_test(TFile* f, Int_t centrality){
    TCanvas* cc = new TCanvas("cc", "cc", 700, 500);
    cc->Divide(2, 2);

    cc->cd(1); 
    TH1D* h1_b = (TH1D*)f->Get("h1d_before_Corrections_FF_PVZNeg_ChPos_BetaPhi");
    TH1D* h1_a = (TH1D*)f->Get("h1d_after_Corrections_FF_PVZNeg_ChPos_BetaPhi");
    h1_a->SetLineColor(2);
    h1_b->Draw();
    h1_a->Draw("same");

    cc->cd(2); 
    TH1D* h2_b = (TH1D*)f->Get("h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi");
    TH1D* h2_a = (TH1D*)f->Get("h1d_after_Corrections_FF_PVZPos_ChPos_BetaPhi");
    h2_a->SetLineColor(2);
    h2_b->Draw();
    h2_a->Draw("same");

    cc->cd(3); 
    TH1D* h3_b = (TH1D*)f->Get("h1d_before_Corrections_RF_PVZNeg_ChPos_BetaPhi");
    TH1D* h3_a = (TH1D*)f->Get("h1d_after_Corrections_RF_PVZNeg_ChPos_BetaPhi");
    h3_a->SetLineColor(2);
    h3_b->Draw();
    h3_a->Draw("same");

    cc->cd(4); 
    TH1D* h4_b = (TH1D*)f->Get("h1d_before_Corrections_RF_PVZPos_ChPos_BetaPhi");
    TH1D* h4_a = (TH1D*)f->Get("h1d_after_Corrections_RF_PVZPos_ChPos_BetaPhi");
    h4_a->SetLineColor(2);
    h4_b->Draw();
    h1_a->Draw("same");

    char plotname[50];
    sprintf(plotname, "plots/lambda_bin_flatness_test_cent%d.png", centrality);
    cc->SaveAs(plotname);
}

void gamma_test(TFile* f, Double_t ep_res, Int_t centrality){
    TCanvas* c = new TCanvas("c_gamma", "c_gamma", 700, 500);
    TProfile* prof = (TProfile*)f->Get("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma");
    prof->Sumw2();
    prof->Scale(1. / 100 / ep_res);
    //Double_t ss = prof->GetBinContent(1);
    //Double_t ss_err = prof->GetBinError(1);
    //Double_t os = prof->GetBinContent(2);
    //Double_t os_err = prof->GetBinError(2);

    //TGraphErrors
    prof->Draw();
    char plotname[50];
    sprintf(plotname, "plots/gamma_cent%d.png", centrality);
    gPad->SaveAs(plotname);
}

void gamma_eff_test(TFile* f, Double_t ep_res, Int_t centrality){
    TCanvas* c = new TCanvas("c_gamma_eff", "c_gamma_eff", 700, 500);
    TProfile* prof = (TProfile*)f->Get("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff");
    prof->Sumw2();
    prof->Scale(1. / 100 / ep_res);
    //Double_t ss = prof->GetBinContent(1);
    //Double_t ss_err = prof->GetBinError(1);
    //Double_t os = prof->GetBinContent(2);
    //Double_t os_err = prof->GetBinError(2);

    //TGraphErrors
    prof->Draw();
    char plotname[50];
    sprintf(plotname, "plots/gamma_eff_cent%d.png", centrality);
    gPad->SaveAs(plotname);
}

Double_t chi(double res) {
  // Calculates chi from the event plane resolution
  double chi = 2.0;
  double delta = 1.0;

  for (int i = 0; i < 15; i++) {
      while(resEventPlane(chi) < res) {chi += delta;}
      delta = delta / 2.;
      while(resEventPlane(chi) > res) {chi -= delta;}
      delta = delta / 2.;
  }

  return chi;
}

Double_t resEventPlane(double chi) {
  // Calculates the event plane resolution as a function of chi

  double con = 0.626657;                   // sqrt(pi/2)/2
  double arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) +
                                          TMath::BesselI1(arg));
  return res;
}
