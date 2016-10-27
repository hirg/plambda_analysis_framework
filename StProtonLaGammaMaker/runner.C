void runner(){
    gROOT->LoadMacro("./StRefMultCorr.cxx++");
    gROOT->LoadMacro("./StCorrelationMaker.cpp++");
    gROOT->LoadMacro("./StProtonLaGammaMaker.cpp++");
    gROOT->LoadMacro("./StEffMaker.h++");
    gROOT->LoadMacro("./StErfEffMaker.h++");
    gROOT->LoadMacro("./StTpcTofEffMaker.h++");

    // Proton Efficiency
    StTpcTofEffMaker* eMaker_Alpha = new StTpcTofEffMaker();
    eMaker_Alpha->SetCentBins(9);
    eMaker_Alpha->SetEffFileName("./auau200GeV_run11_p_eff.fit.dat"); //something like that
    eMaker_Alpha->Init();

    // Proton-bar Efficiency
    StTpcTofEffMaker* eMaker_AlphaBar = new StTpcTofEffMaker();
    eMaker_AlphaBar->SetCentBins(9);
    eMaker_AlphaBar->SetEffFileName("./auau200GeV_run11_p_eff.fit.dat"); //same for pbar and p
    eMaker_AlphaBar->Init();

    // Lambda Efficiency
    StErfEffMaker* eMaker_Beta = new StErfEffMaker();
    eMaker_Beta->SetCentBins(9);
    eMaker_Beta->SetPtThreshold(1.6);
    eMaker_Beta->SetEffFileName("./auau200GeV_run11_la_eff_comb.fit.dat"); //something like that
    eMaker_Beta->Init();

    // Lambda-bar Efficiency
    StErfEffMaker* eMaker_BetaBar = new StErfEffMaker();
    eMaker_BetaBar->SetCentBins(9);
    eMaker_BetaBar->SetPtThreshold(1.6);
    eMaker_BetaBar->SetEffFileName("./auau200GeV_run11_antila_eff_comb.fit.dat"); //something like that
    eMaker_BetaBar->Init();

    StProtonLaGammaMaker* maker = new StProtonLaGammaMaker(1, 300000, eMaker_Alpha, eMaker_Beta, eMaker_AlphaBar, eMaker_BetaBar);
    maker->Init();
    maker->ComputePhiCorrectionParams();
    maker->reconstructEventPlaneWithPhiWeightCorrection();//TODO
    maker->reconstructShiftedSubEventPlane();//TODO
    maker->reconstructShiftedFullEventPlaneAndComputeCorrelator();//TODO: change the function name to Capital
    maker->Finish();
}
