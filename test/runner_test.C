void runner_test(){
    gROOT->LoadMacro("./StRefMultCorr.cxx++");
    gROOT->LoadMacro("./StCorrelationMaker.cpp++");
    gROOT->LoadMacro("./StProtonLaGammaMaker.cpp++");
    StProtonLaGammaMaker* maker = new StProtonLaGammaMaker(3, 200000);
    maker->Init();
    maker->ComputePhiCorrectionParams();
    maker->reconstructEventPlaneWithPhiWeightCorrection();
    maker->reconstructShiftedSubEventPlane();
    maker->reconstructShiftedFullEventPlaneAndComputeCorrelator();
    maker->Finish();
}
