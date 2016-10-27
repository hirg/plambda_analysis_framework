void compile(int cen){
    gROOT->LoadMacro("StRefMultCorr/StRefMultCorr.cxx++");
    gROOT->LoadMacro("StPriTrkCuts/StPriTrkCuts.cpp++");
    gROOT->LoadMacro("StPriTrkGeneralCuts/StPriTrkGeneralCuts.cpp++");
    gROOT->LoadMacro("StPriPTrkCuts/StPriPTrkCuts.cpp++");
    gROOT->LoadMacro("StEvtInfo/StEvtInfo.cpp++");
    gROOT->LoadMacro("StEvtCuts/StEvtCuts.cpp++");
    gROOT->LoadMacro("StV0TrkCuts/StV0TrkCuts.cpp++");
    gROOT->LoadMacro("StProtonLaGammaMaker/StCorrelationMaker.cpp++");
    gROOT->LoadMacro("StProtonLaGammaMaker/StProtonLaGammaMaker.cpp++");
    gROOT->LoadMacro("StProtonLaGammaMaker/StEffMaker.h++");
    gROOT->LoadMacro("StProtonLaGammaMaker/StErfEffMaker.h++");
    gROOT->LoadMacro("StProtonLaGammaMaker/StTpcTofEffMaker.h++");

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

    // Configure event cuts
    const int Nrun_MB1 = 134;
    const int Nrun_MB2 = 60;
    const int Nrun_MB5 = 107;
    const int Nrun_MB6 = 42;

    const int bad_Ref_day3_MB1_auau200GeV[Nrun_MB1] = {133010,133011,133022,133027,133028,134065,127019,128021,133020,137004,127039,132024,132025,132026,132032,127020,127048,127049,128028,128029,128030,128031,128032,132052,133041,136039,127002,132061,134017,134063,135004,135039,137020,127030,128007,132062,133005,133039,133040,133053,133054,134006,134007,134018,134026,134028,134038,134041,135019,135020,135045,135057,136007,136032,136064,137010,127009,127046,128038,132019,132022,132044,132045,133002,133019,133038,133052,134005,134040,134055,135002,135033,135048,135049,136069,137003,127021,127022,127023,127024,128024,128025,132020,132021,132023,132033,132048,132051,132057,132063,132065,133021,134008,135012,135021,135024,135030,135054,136044,136081,136086,127003,127010,127011,127017,127018,127032,132009,132034,132043,132066,132069,133018,134023,134057,136005,136006,136014,136017,136022,136023,136024,136025,136027,136028,136029,136030,136031,136034,136053,136054,136070,136071,138017}; //MB1
    const int bad_Ref_day3_MB2_auau200GeV[Nrun_MB2] = {139032,139043,139044,139045,142002,139042,140021,140029,142063,142064,142065,144004,138081,138082,138087,138088,138089,138090,138091,139002,139003,139006,139007,139008,139009,139010,139015,139016,139017,139018,139021,142016,142033,142061,142062,144051,138092,139019,139020,140016,140020,141003,141004,141026,141062,141065,142001,142013,142023,142034,142046,142068,142076,143009,143024,143058,144016,144028,144033,145003}; //MB2
    const int bad_Ref_day3_MB5_auau200GeV[Nrun_MB5] = {155050,155056,158010,165028,154043,154044,154045,155058,158069,158070,158072,158073,164067,154046,154047,155008,155009,156015,156062,156063,158074,162015,154067,155002,155012,155047,156008,156009,157023,157030,157052,158006,159023,160021,161006,161015,161060,162004,162028,162034,163024,163058,164009,164056,164066,164089,165013,154048,154066,155011,155021,155038,155051,155060,155062,155064,156004,156035,156056,157012,157014,157031,157038,157051,158015,158021,158026,158040,158041,158051,158054,158056,158057,158058,158061,159005,159021,159022,159024,160016,160025,161007,161014,161017,161020,161022,161053,162017,162030,162035,162055,162056,162057,162058,163006,163008,163015,164011,164037,164043,164086,165001,165003,165005,165007,165026,165031}; //MB5
    const int bad_Ref_day3_MB6_auau200GeV[Nrun_MB6] = {166051,170016,167014,170034,170050,170051,171009,171015,167049,168010,169028,169032,170007,165042,166052,166059,167002,167040,167048,169031,169059,170009,170012,170018,170020,170031,171004,171014,166002,166003,167015,167024,168009,168022,168060,168077,169033,169034,170044,170045,170054,170056}; //MB6

    set<int> badruns; 
    for(int i = 0; i != Nrun_MB5; i++)
        badruns.insert(bad_Ref_day3_MB5_auau200GeV[i]); //TODO: Need to be changed when you change dataset

    set<int> cent;
/*
    cent.insert(1);
    cent.insert(2);
    cent.insert(3);
    cent.insert(4);
    cent.insert(5);
*/
    cent.insert(cen);

    StEvtCuts* evt_cuts = new StEvtCuts();
    evt_cuts->LoadBadRunsDay3(badruns);
    evt_cuts->SetCentralities(cent);
    evt_cuts->SetRunLowerBound(0);
    evt_cuts->SetRunUpperBound(99999999);

    // Configure track cuts
    /* General primary tracks:
     *   PassEtaCuts(p) 
     *&& PassPtCuts(p)
     *&& PassDcaCuts(p)
     *&& PassFlagCuts(p)
     *&& PassPIDCuts(p) = 0 
     */

    StPriTrkGeneralCuts* general_cuts = new StPriTrkGeneralCuts();
    general_cuts->SetPtLowerBound(.15);
    general_cuts->SetPtUpperBound(2.0);
    general_cuts->SetEtaLowerBound(-1.0);
    general_cuts->SetEtaUpperBound(1.0);
    general_cuts->SetDcaLowerBound(-1.0);
    general_cuts->SetDcaUpperBound(2.0);

    StPriPTrkCuts* p_cuts = new StPriPTrkCuts();
    p_cuts->SetPtLowerBound(.4);
    p_cuts->SetPtUpperBound(2.0);
    p_cuts->SetEtaLowerBound(-1.0);
    p_cuts->SetEtaUpperBound(1.0);
    p_cuts->SetDcaLowerBound(-1.0);
    p_cuts->SetDcaUpperBound(1.0);
    p_cuts->SetNSigmaPLowerBound(-2.0);
    p_cuts->SetNSigmaPUpperBound(2.0);
    p_cuts->SetM2LowerBound(0.8);
    p_cuts->SetM2UpperBound(1.0);

    StV0TrkCuts* lambda_cuts = new StV0TrkCuts();
    lambda_cuts->SetPtLowerBound(.4);
    lambda_cuts->SetPtUpperBound(5.0);
    lambda_cuts->SetEtaLowerBound(-1.0);
    lambda_cuts->SetEtaUpperBound(1.0);
    lambda_cuts->SetDcaLowerBound(-1.0);
    lambda_cuts->SetDcaUpperBound(.6);
    lambda_cuts->SetDau1NSigmaLowerBound(-3.0);
    lambda_cuts->SetDau1NSigmaUpperBound(3.0);
    lambda_cuts->SetDau2NSigmaLowerBound(-3.0);
    lambda_cuts->SetDau2NSigmaUpperBound(3.0);
    lambda_cuts->SetDecLenLowerBound(6.0);
    lambda_cuts->SetDecLenUpperBound(999999.0);
    lambda_cuts->SetDau1DcaLowerBound(0.6);
    lambda_cuts->SetDau1DcaUpperBound(999999.0);
    lambda_cuts->SetDau2DcaLowerBound(1.8);
    lambda_cuts->SetDau2DcaUpperBound(999999.0);
    lambda_cuts->SetDca1to2LowerBound(-1.0);
    lambda_cuts->SetDca1to2UpperBound(.7);
    lambda_cuts->SetMassLowerBound(1.115683 - 0.004);
    lambda_cuts->SetMassUpperBound(1.115683 + 0.004);
    lambda_cuts->Dump();

    StProtonLaGammaMaker* maker = new StProtonLaGammaMaker(cen, 0, eMaker_Alpha, eMaker_Beta, eMaker_AlphaBar, eMaker_BetaBar);

    pair<int, int> pair_day2(8000, 18000); // TODO: change
    pair<int, int> pair_day(80, 180); //TODO: This is for event-plane correction
    maker->SetDayBoundaries(pair_day);
    maker->SetDay2Boundaries(pair_day2);

    set<string> dataset;
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data154/*.lambda.picodst.root"); 
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data155/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data156/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data157/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data158/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data159/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data160/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data161/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data162/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data163/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data164/*.lambda.picodst.root");
    dataset.insert("/media/Disk_Yan/liwen/Run11_200GeV/Data165/*.lambda.picodst.root");
    maker->LoadDataSet(dataset);

    maker->LoadEvtCuts(evt_cuts);
    maker->LoadPriTrkCuts(general_cuts);
    maker->LoadProtonTrkCuts(p_cuts);
    maker->LoadLambdaTrkCuts(lambda_cuts);

    maker->Init();
    maker->ComputePhiCorrectionParams();
    maker->ReconstructEventPlaneWithPhiWeightCorrection();//TODO
    maker->ReconstructShiftedSubEventPlane();//TODO
    maker->ReconstructShiftedFullEventPlaneAndComputeCorrelator();//TODO: change the function name to Capital
    maker->Finish();
}
