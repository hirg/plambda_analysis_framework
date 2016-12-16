#include "TChain.h"
#include "TVector2.h"
#include <string>
#include "TH1.h"
#include "TObject.h"
#include <cstdlib>
#include <iostream>
#include "TFile.h"
#include "../include/constants.h"
#include "TTree.h"
#include "TLeaf.h"
#include "StCorrelationMaker.h"
#include "../StRefMultCorr/StRefMultCorr.h"
//#include "StRefMultCorr.h"
#include "TMath.h" 
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include <map>
#include <utility>

ClassImp(StCorrelationMaker)
StCorrelationMaker::StCorrelationMaker(Int_t cent, Int_t eventNumberToProcess, StEffMaker* e_alpha, StEffMaker* e_beta, StEffMaker* e_alphabar, StEffMaker* e_betabar): m_NEvents(-1), m_Centrality(cent), m_EventNumberToProcess(eventNumberToProcess), m_Chain(new TChain("StrangenessDst")), m_EffMakerAlpha(e_alpha), m_EffMakerBeta(e_beta), m_EffMakerAlphaBar(e_alphabar), m_EffMakerBetaBar(e_betabar){

    if(m_Centrality <= 0){
        cout << "ERROR: Bad Centrality Number!!" << endl;
        exit(1);
    }
    char weightFilename[50];
    sprintf(weightFilename, "weights0.cen%d_%ld.root", m_Centrality, m_EventNumberToProcess);
    m_Weights0File = new TFile(weightFilename, "recreate");

    sprintf(weightFilename, "weights1.cen%d_%ld.root", m_Centrality, m_EventNumberToProcess);
    m_Weights1File = new TFile(weightFilename, "recreate");

    sprintf(weightFilename, "weights2.cen%d_%ld.root", m_Centrality, m_EventNumberToProcess);
    m_Weights2File = new TFile(weightFilename, "recreate");

    char resultFilename[50];
    sprintf(resultFilename, "Cen%d_%ld.result.root", m_Centrality, m_EventNumberToProcess);
    m_ResultFile = new TFile(resultFilename, "recreate");
}

// Step Zero: initialize histograms 
Int_t StCorrelationMaker::Init(){
    // Initialize generic histograms and customized histograms 
    initHistograms();

    addTrees();
    m_NEvents = m_Chain->GetEntriesFast();
    std::cout << ">>>>>>>>>>>>Number of events loaded: " << m_NEvents << "<<<<<<<<<<<<<<<<<<" << std::endl;

    if(m_NEvents < m_EventNumberToProcess)
        m_EventNumberToProcess = m_NEvents;

    // Used to iterate over all of the events
    if(m_EventNumberToProcess <= 0)
        m_EventNumberToProcess = m_NEvents;

    zeroWeightsHistograms();

    m_EvtCuts->Dump();
    return m_NEvents;
}

// Step One: compute phi weights 
Int_t StCorrelationMaker::ComputePhiCorrectionParams(){
    std::cout << ">>>>>>>>Phi Correction Starts...<<<<<<<<<" << std::endl;
    StRefMultCorr refmultCorrUtil = StRefMultCorr("refmult");

    for(Long_t i = 0; i < m_EventNumberToProcess; i++){
	if(i % 10000 == 0) std::cout << i << " Events Processed" <<  std::endl;
	m_Chain->GetEntry(i);
	StEvtInfo evtInfo(m_Chain, refmultCorrUtil);
	if(!m_EvtCuts->PassAllCuts(evtInfo)) continue;

	((TH1D*)(m_HistsCollection.find("h1d_CentralityUnWeighted")->second))->Fill(evtInfo.Centrality());
	((TH1D*)(m_HistsCollection.find("h1d_CentralityWeighted")->second))->Fill(evtInfo.Centrality(), evtInfo.EWeight());
	((TH1D*)(m_HistsCollection.find("h1d_PVertexZUnWeighted")->second))->Fill(evtInfo.Vz());
	((TH1D*)(m_HistsCollection.find("h1d_PVertexZWeighted")->second))->Fill(evtInfo.Vz(), evtInfo.EWeight());
	((TH1D*)(m_HistsCollection.find("h1d_EventTally")->second))->Fill("Total Event", 1);

	// Filling particle histograms
	fillPrimaryTracksPhiHists(m_HistsCollection, evtInfo);
	fillAlphaHists(m_HistsCollection, evtInfo);
	fillBetaHists(m_HistsCollection, evtInfo);
    }
    computePhiWeights(m_HistsCollection);

    // Save phi weights histograms
    saveWeights0File();
}

void StCorrelationMaker::ReconstructEventPlaneWithPhiWeightCorrection(){
    std::cout << ">>>>>>>>Event Plane Reconstruction Starts...<<<<<<<<<" << std::endl;
    StRefMultCorr refmultCorrUtil = StRefMultCorr("refmult");
    for(Long_t i = 0; i < m_EventNumberToProcess; i++){
	if( i % 10000 == 0) std::cout << i << " Events Processed" <<  std::endl;
	m_Chain->GetEntry(i);
	StEvtInfo evtInfo(m_Chain, refmultCorrUtil);
	if(!m_EvtCuts->PassAllCuts(evtInfo)) continue;
	reconstructSubEventPlaneWithPhiWeightHelper(m_HistsCollection, evtInfo);
    }
}

void StCorrelationMaker::ReconstructShiftedSubEventPlane(){
    std::cout << ">>>>>>>>Shifted Sub Event Plane Reconstruction Starts...<<<<<<<<<" << std::endl;
    StRefMultCorr refmultCorrUtil = StRefMultCorr("refmult");
    for(Long_t i = 0; i < m_EventNumberToProcess; i++){
	if(i % 10000 == 0) std::cout << i << " Events Processed" <<  std::endl;
	m_Chain->GetEntry(i);
	StEvtInfo evtInfo(m_Chain, refmultCorrUtil);
	if(!m_EvtCuts->PassAllCuts(evtInfo)) continue;
	reconstructShiftedSubEventPlaneHelper(m_HistsCollection, evtInfo);
    }
}

void StCorrelationMaker::ReconstructShiftedFullEventPlaneAndComputeCorrelator(){
    std::cout << ">>>>>>>>Shifted Full Event Plane Reconstruction Starts...<<<<<<<<<" << std::endl;
    StRefMultCorr refmultCorrUtil = StRefMultCorr("refmult");
    for(Long_t i = 0; i < m_EventNumberToProcess; i++){
	if(i % 10000 == 0) std::cout << i << " Events Processed" <<  std::endl;
	m_Chain->GetEntry(i);
	StEvtInfo evtInfo(m_Chain, refmultCorrUtil);
	if(!m_EvtCuts->PassAllCuts(evtInfo)) continue;
	TVector2 shiftedEPPhi_full;
	reconstructShiftedFullEventPlaneHelper(m_HistsCollection, evtInfo, shiftedEPPhi_full);
	computeCorrelatorsHelper(m_HistsCollection, evtInfo, shiftedEPPhi_full);
    }
}

void StCorrelationMaker::Finish(){
    m_ResultFile->cd();
    for(std::map<std::string, TH1*>::iterator it = m_HistsCollection.begin(); it != m_HistsCollection.end(); it++)
	it->second->Write();

    m_ResultFile->Write();
    m_ResultFile->Close();
}

//============================ Private member functions ==============================
Int_t StCorrelationMaker::initHistograms(){
    //m_ResultFile->cd();

    m_Weights0File->cd();
    TH1D* h1d_EventTally = new TH1D("h1d_EventTally", "Event Tally", 10, 0, 1); 
    h1d_EventTally->SetBit(TH1::kCanRebin);
    h1d_EventTally->SetStats(0);
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_EventTally", h1d_EventTally));

    TH1D* h1d_TriggerId = new TH1D("h1d_TriggerId", "Trigger Id", 200, 0.5, 200.5);//TODO:??
    TH1D* h1d_CentralityUnWeighted = new TH1D("h1d_CentralityUnWeighted", "Centrality UnWeighted", 10, 0, 10);
    TH1D* h1d_CentralityWeighted = new TH1D("h1d_CentralityWeighted", "Centrality Weighted", 10, 0, 10);
    TH1D* h1d_PVertexZUnWeighted = new TH1D("h1d_PVertexZUnWeighted", "PVertex Z UnWeighted", 100, -100, 100);
    TH1D* h1d_PVertexZWeighted = new TH1D("h1d_PVertexZWeighted", "PVertex Z Weighted", 100, -100, 100);
    TH2D* h2d_Mult_vs_PVertexZUnWeighted = new TH2D("h2d_Mult_vs_PVertexZUnWeighted", "Mult. vs PVertex UnWeighted", 1000, -.5, 999.5, 100, -100, 100);
    TH2D* h2d_Mult_vs_PVertexZWeighted = new TH2D("h2d_Mult_vs_PVertexZWeighted", "Mult. vs PVertex Weighted", 1000, -.5, 999.5, 100, -100, 100);
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_TriggerId", h1d_TriggerId));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_CentralityUnWeighted", h1d_CentralityUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_CentralityWeighted", h1d_CentralityWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_PVertexZUnWeighted", h1d_PVertexZUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_PVertexZWeighted", h1d_PVertexZWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_Mult_vs_PVertexZUnWeighted", h2d_Mult_vs_PVertexZUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_Mult_vs_PVertexZWeighted", h2d_Mult_vs_PVertexZWeighted));

    // Histograms used for quality assurance of event plane orientation correction
    TH2D* h2d_RawTPCEPEastPhi_vs_Day = new TH2D("h2d_RawTPCEPEastPhi_vs_Day", "Raw TPC Event Plane East Phi vs Day", 36, 0, PI, m_DayUpperBound - m_DayLowerBound, m_DayLowerBound, m_DayUpperBound);
    TH2D* h2d_RawTPCEPWestPhi_vs_Day = new TH2D("h2d_RawTPCEPWestPhi_vs_Day", "Raw TPC Event Plane West Phi vs Day", 36, 0, PI, m_DayUpperBound - m_DayLowerBound, m_DayLowerBound, m_DayUpperBound);
    TH2D* h2d_RawTPCEPFullPhi_vs_Day = new TH2D("h2d_RawTPCEPFullPhi_vs_Day", "Raw TPC Event Plane Full Phi vs Day", 36, 0, PI, m_DayUpperBound - m_DayLowerBound, m_DayLowerBound, m_DayUpperBound);
    TH1F* h1f_before_Flattened_EastEPPhi = new TH1F("h1f_before_Flattened_EastEPPhi", "TPC East Event Plane Phi before Flattening", 100, 0, PI);
    TH1F* h1f_before_Flattened_WestEPPhi = new TH1F("h1f_before_Flattened_WestEPPhi", "TPC West Event Plane Phi before Flattening", 100, 0, PI);
    TH1F* h1f_before_Flattened_FullEPPhi = new TH1F("h1f_before_Flattened_FullEPPhi", "TPC Full Event Plane Phi before Flattening", 100, 0, PI);
    TH1F* h1f_Flattened_EastEPPhi = new TH1F("h1f_Flattened_EastEPPhi", "TPC East Event Plane Phi after Flattening", 100, 0, PI);
    TH1F* h1f_Flattened_WestEPPhi = new TH1F("h1f_Flattened_WestEPPhi", "TPC West Event Plane Phi after Flattening", 100, 0, PI);
    TH1F* h1f_Flattened_FullEPPhi = new TH1F("h1f_Flattened_FullEPPhi", "TPC Full Event Plane Phi after Flattening", 100, 0, PI);

    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_RawTPCEPEastPhi_vs_Day", h2d_RawTPCEPEastPhi_vs_Day));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_RawTPCEPWestPhi_vs_Day", h2d_RawTPCEPWestPhi_vs_Day));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_RawTPCEPFullPhi_vs_Day", h2d_RawTPCEPFullPhi_vs_Day));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1f_Flattened_EastEPPhi", h1f_Flattened_EastEPPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1f_Flattened_WestEPPhi", h1f_Flattened_WestEPPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1f_Flattened_FullEPPhi", h1f_Flattened_FullEPPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1f_before_Flattened_EastEPPhi", h1f_before_Flattened_EastEPPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1f_before_Flattened_WestEPPhi", h1f_before_Flattened_WestEPPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1f_before_Flattened_FullEPPhi", h1f_before_Flattened_FullEPPhi));

    // Event plane resolution
    TProfile* profile_eventplane_resolution = new TProfile("profile_eventplane_resolution", "Event Plane Resolution", 1, .5, 1.5, -100, 100, "");
    TProfile* profile_eventplane_resolution_noEW = new TProfile("profile_eventplane_resolution_noEW", "Event Plane Resolution", 1, .5, 1.5, -100, 100, "");
    m_HistsCollection.insert(std::pair<std::string, TH1*>("profile_eventplane_resolution", profile_eventplane_resolution));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("profile_eventplane_resolution_noEW", profile_eventplane_resolution_noEW));

    // Histograms to store correction terms. Used to flatten event-plane orientation distribution. Odd-cos(2, 4, 6, 8*TPC_EP), Even-sin(2, 4, 6, 8*TPC_EP)
    TProfile2D* prof2_XOrder_YDay_ZCorrectionTerm_EastEP = new TProfile2D("prof2_XOrder_YDay_ZCorrectionTerm_EastEP", "XOrder YDay ZCorrectionTerm East", 2 * order, 0.5, .5 + 2 * order, m_DayUpperBound - m_DayLowerBound, m_DayLowerBound, m_DayUpperBound, -1, 1, ""); 
    TProfile2D* prof2_XOrder_YDay_ZCorrectionTerm_WestEP = new TProfile2D("prof2_XOrder_YDay_ZCorrectionTerm_WestEP", "XOrder YDay ZCorrectionTerm West", 2 * order, 0.5, .5 + 2 * order, m_DayUpperBound - m_DayLowerBound, m_DayLowerBound, m_DayUpperBound, -1, 1, ""); 
    TProfile2D* prof2_XOrder_YDay_ZCorrectionTerm_FullEP = new TProfile2D("prof2_XOrder_YDay_ZCorrectionTerm_FullEP", "XOrder YDay ZCorrectionTerm Full", 2 * order, 0.5, .5 + 2 * order, m_DayUpperBound - m_DayLowerBound, m_DayLowerBound, m_DayUpperBound, -1, 1, ""); 
    m_HistsCollection.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZCorrectionTerm_EastEP", prof2_XOrder_YDay_ZCorrectionTerm_EastEP));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZCorrectionTerm_WestEP", prof2_XOrder_YDay_ZCorrectionTerm_WestEP));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZCorrectionTerm_FullEP", prof2_XOrder_YDay_ZCorrectionTerm_FullEP));

    // Histograms used to QA
    TH1D* h1d_PtUnWeighted = new TH1D("h1d_PtUnWeighted", "Pt UnWeighted", 300, 0, 15);
    TH1D* h1d_PtWeighted = new TH1D("h1d_PtWeighted", "Pt Weighted", 300, 0, 15);
    TH2D* h2d_Eta_vs_PtUnWeighted = new TH2D("h2d_Eta_vs_PtUnWeighted", "Eta vs Pt UnWeighted", 30, -1.5, 1.5, 300, 0, 15);  //PrimaryTracks
    TH2D* h2d_Eta_vs_PtWeighted = new TH2D("h2d_Eta_vs_PtWeighted", "Eta vs Pt Weighted", 30, -1.5, 1.5, 300, 0, 15);  //PrimaryTracks
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_PtUnWeighted", h1d_PtUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_PtWeighted", h1d_PtWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_Eta_vs_PtUnWeighted", h2d_Eta_vs_PtUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_Eta_vs_PtWeighted", h2d_Eta_vs_PtWeighted));

    std::vec<std::string> features;
    features.push_back("FF_PVZPos_ChPos");
    features.push_back("FF_PVZPos_ChNeg");
    features.push_back("FF_PVZNeg_ChPos");
    features.push_back("FF_PVZNeg_ChNeg");
    features.push_back("RF_PVZPos_ChPos");
    features.push_back("RF_PVZPos_ChNeg");
    features.push_back("RF_PVZNeg_ChPos");
    features.push_back("RF_PVZNeg_ChNeg");

    for(int i = 0; i != features.size(); ++i){
        char histname[100]; 

	// phi distribution before correction
        sprintf(histname, "h1d_before_Corrections_%s_PrimaryTrkPhi", features[i].c_str());
	TH1D* h1d_before_tmp = new TH1D(histname, histname, phiBins, -PI, PI);
	m_HistsCollection.insert(std::pair<std::string, TH1*>(histname, h1d_before_tmp));

	// phi distribution after correction
        sprintf(histname, "h1d_after_Corrections_%s_PrimaryTrkPhi", features[i].c_str());
	TH1D* h1d_before_tmp = new TH1D(histname, histname, phiBins, -PI, PI);
	m_HistsCollection.insert(std::pair<std::string, TH1*>(histname, h1d_before_tmp));
    }

    addAlphaHists(m_HistsCollection); // Phi distribution + invariant mass
    addBetaHists(m_HistsCollection);

    addCorrelationHists(m_HistsCollection);

    return 1;
}

Int_t StCorrelationMaker::addTrees(){
    Int_t nFiles = 0;
    for(set<string>::iterator iter = m_DataSet.begin(); iter != m_DataSet.end(); iter++)
	nFiles += m_Chain->Add((*iter).c_str(), 0); // Store each header in memory 
    return nFiles;
}

void StCorrelationMaker::addAlphaHists(std::map<std::string, TH1*>& histMap){
    std::cout << histMap.size() << "size before addAlphaHists" << std::endl;
    // Used to check distribution before correction 
    std::vec<std::string> features;
    features.push_back("FF_PVZPos_ChPos");
    features.push_back("FF_PVZPos_ChNeg");
    features.push_back("FF_PVZNeg_ChPos");
    features.push_back("FF_PVZNeg_ChNeg");
    features.push_back("RF_PVZPos_ChPos");
    features.push_back("RF_PVZPos_ChNeg");
    features.push_back("RF_PVZNeg_ChPos");
    features.push_back("RF_PVZNeg_ChNeg");

    for(int i = 0; i != features.size(); ++i){
        char histname[100]; 

	// phi distribution before correction
        sprintf(histname, "h1d_before_Corrections_%s_AlphaPhi", features[i].c_str());
	TH1D* h1d_before_tmp = new TH1D(histname, histname, phiBins, -PI, PI);
	histMap.insert(std::pair<std::string, TH1*>(histname, h1d_before_tmp));

	// Used to fill correction terms to flatten the distribution of proton's phi
        sprintf(histname, "prof2_XOrder_YDay_ZPhi_%s_Alpha", features[i].c_str());
	TProfile2D* prof2_correction_term_tmp = new TProfile2D(histname, histname, 2 * order, .5, 2 * order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
	histMap.insert(std::pair<std::string, TH1*>(histname, prof2_correction_term_tmp));
	// TODO:store in a weight file

	// phi distribution after correction
        sprintf(histname, "h1d_after_Corrections_%s_AlphaPhi", features[i].c_str());
	TH1D* h1d_before_tmp = new TH1D(histname, histname, phiBins, -PI, PI);
	histMap.insert(std::pair<std::string, TH1*>(histname, h1d_before_tmp));
    }

    // Used to check after correction distribution
    TH1F* h1f_before_Corrections_AlphaPhi = new TH1F("h1f_before_Corrections_AlphaPhi", "h1f_before_Corrections_AlphaPhi", phiBins, -PI, PI);
    TH1F* h1f_after_Corrections_AlphaPhi = new TH1F("h1f_after_Corrections_AlphaPhi", "h1f_after_Corrections_AlphaPhi", phiBins, -PI, PI);
    histMap.insert(std::pair<std::string, TH1*>("h1f_before_Corrections_AlphaPhi", h1f_before_Corrections_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1f_after_Corrections_AlphaPhi", h1f_after_Corrections_AlphaPhi));

    addAdditionalAlphaHists(histMap);
    std::cout << histMap.size() << "size after addAdditionalAlphaHists" << std::endl;
}

void StCorrelationMaker::addBetaHists(std::map<std::string, TH1*>& histMap){
    std::vec<std::string> features;
    features.push_back("FF_PVZPos_ChPos");
    features.push_back("FF_PVZPos_ChNeg");
    features.push_back("FF_PVZNeg_ChPos");
    features.push_back("FF_PVZNeg_ChNeg");
    features.push_back("RF_PVZPos_ChPos");
    features.push_back("RF_PVZPos_ChNeg");
    features.push_back("RF_PVZNeg_ChPos");
    features.push_back("RF_PVZNeg_ChNeg");

    for(int i = 0; i != features.size(); ++i){
        char histname[100]; 

	// phi distribution before correction
        sprintf(histname, "h1d_before_Corrections_%s_BetaPhi", features[i].c_str());
	TH1D* h1d_before_tmp = new TH1D(histname, histname, phiBins, -PI, PI);
	histMap.insert(std::pair<std::string, TH1*>(histname, h1d_before_tmp));

	// Used to fill correction terms to flatten the distribution of proton's phi
        sprintf(histname, "prof2_XOrder_YDay_ZPhi_%s_Beta", features[i].c_str());
	TProfile2D* prof2_correction_term_tmp = new TProfile2D(histname, histname, 2 * order, .5, 2 * order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
	histMap.insert(std::pair<std::string, TH1*>(histname, prof2_correction_term_tmp));
	// TODO:store in a weight file

	// phi distribution after correction
        sprintf(histname, "h1d_after_Corrections_%s_BetaPhi", features[i].c_str());
	TH1D* h1d_before_tmp = new TH1D(histname, histname, phiBins, -PI, PI);
	histMap.insert(std::pair<std::string, TH1*>(histname, h1d_before_tmp));
    }
    // Used to fill correction terms to flatten the distribution of Lambda's phi
    // Used to check after correction distribution
    TH1F* h1f_before_Corrections_BetaPhi = new TH1F("h1f_before_Corrections_BetaPhi", "h1f_before_Corrections_BetaPhi", phiBins, -PI, PI);
    TH1F* h1f_after_Corrections_BetaPhi = new TH1F("h1f_after_Corrections_BetaPhi", "h1f_after_Corrections_BetaPhi", phiBins, -PI, PI);
    histMap.insert(std::pair<std::string, TH1*>("h1f_before_Corrections_BetaPhi", h1f_before_Corrections_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1f_after_Corrections_BetaPhi", h1f_after_Corrections_BetaPhi));

    addAdditionalBetaHists(histMap);
}

void StCorrelationMaker::zeroWeightsHistograms(){
    vec<string> features;
    features.push_back("FF_PVZPos_ChPos");
    features.push_back("FF_PVZPos_ChNeg");
    features.push_back("FF_PVZNeg_ChPos");
    features.push_back("FF_PVZNeg_ChNeg");
    features.push_back("RF_PVZPos_ChPos");
    features.push_back("RF_PVZPos_ChNeg");
    features.push_back("RF_PVZNeg_ChPos");
    features.push_back("RF_PVZNeg_ChNeg");

    string trk_types[3] = {"Alpha", "Beta", "PrimaryTrk"};

    for(int i = 0; i != features.size(); ++i){
	for(int j = 0; j != 3; ++j){
            char histname[100];
	    sprintf(histname, "h1d_before_Corrections_%s_%sPhi", features[i].c_str(), trk_types.c_str());
	    std::string histname_str(histname);
	    TH1D* hist = (TH1D*)m_HistsCollection[histname_str];
	    for(int ii = 0; ii < phiBins; ++ii)
		hist->SetBinContent(ii + 1, 0.0);
	}
    }
}
