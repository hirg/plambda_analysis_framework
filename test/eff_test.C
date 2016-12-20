{
    gROOT->LoadMacro("./StEffMaker.h++");
    gROOT->LoadMacro("./StErfEffMaker.h++");
    gROOT->LoadMacro("./StTpcTofEffMaker.h++");

    // Proton Efficiency
    StTpcTofEffMaker* eMaker_Alpha = new StTpcTofEffMaker(); 
    eMaker_Alpha->SetCentBins(9);
    eMaker_Alpha->SetEffFileName("./auau200GeV_run11_p_eff.fit.dat");
    eMaker_Alpha->Init();
    cout << "Alpha(Proton): " << endl;
    cout << "pt = 0.15, eff = " << eMaker_Alpha->GetEfficiency(1, 0.15) << endl;
    cout << "pt = 0.5, eff = " << eMaker_Alpha->GetEfficiency(1, 0.5) << endl;
    cout << "pt = 1.0, eff = " << eMaker_Alpha->GetEfficiency(1, 1.0) << endl;
    cout << "pt = 10.6, eff = " << eMaker_Alpha->GetEfficiency(1, 10.6) << endl;
    
    // Lambda Efficiency
    StErfEffMaker* eMaker_Beta = new StErfEffMaker();
    eMaker_Beta->SetCentBins(9);
    eMaker_Beta->SetPtThreshold(1.4);
    eMaker_Beta->SetEffFileName("./auau200GeV_run11_la_eff_comb.fit.dat"); //something like that
    eMaker_Beta->Init();
    cout << "Beta(Lambda): " << endl;
    cout << "pt = 0.6, eff = " << eMaker_Beta->GetEfficiency(1, 0.6) << endl;
    eMaker_Beta->SetPtThreshold(0.4);
    cout << "pt = 0.6, eff = " << eMaker_Beta->GetEfficiency(1, 0.6) << endl;
    cout << "pt = 2.6, eff = " << eMaker_Beta->GetEfficiency(1, 2.6) << endl;
    cout << "pt = 10.6, eff = " << eMaker_Beta->GetEfficiency(1, 10.6) << endl;

    // Lambda-bar Efficiency
    StErfEffMaker* eMaker_BetaBar = new StErfEffMaker();
    eMaker_BetaBar->SetCentBins(9);
    eMaker_BetaBar->SetPtThreshold(.6);
    eMaker_BetaBar->SetEffFileName("./auau200GeV_run11_antila_eff_comb.fit.dat"); //something like that
    eMaker_BetaBar->Init();
    cout << "BetaBar(AntiLambda): " << endl;
    cout << "pt = 0.6, eff = " << eMaker_BetaBar->GetEfficiency(1, 0.6) << endl;
    cout << "pt = 2.6, eff = " << eMaker_BetaBar->GetEfficiency(1, 2.6) << endl;
    cout << "pt = 10.6, eff = " << eMaker_BetaBar->GetEfficiency(1, 10.6) << endl;

}
