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
    sprintf(weightFilename, "Cen%d_%ld.weight.root", m_Centrality, m_EventNumberToProcess);
    m_WeightFile = new TFile(weightFilename, "recreate");

    char resultFilename[50];
    sprintf(resultFilename, "Cen%d_%ld.result.root", m_Centrality, m_EventNumberToProcess);
    m_ResultFile = new TFile(resultFilename, "recreate");
}

Int_t StCorrelationMaker::Init(){
    initiateHistograms();

    //initiateComputingParams();
    addTrees();
    m_NEvents = m_Chain->GetEntries();
    std::cout << ">>>>>>>>>>>>Number of events loaded: " << m_NEvents << "<<<<<<<<<<<<<<<<<<" << std::endl;
    if(m_NEvents < m_EventNumberToProcess)
        m_EventNumberToProcess = m_NEvents;

    // Used to iterate over all of the events
    if(m_EventNumberToProcess == 0)
        m_EventNumberToProcess = m_NEvents;
    zeroWeightsHistograms();

    m_EvtCuts->Dump();
    return 1;
}

Int_t StCorrelationMaker::ComputePhiCorrectionParams(){
    std::cout << ">>>>>>>>Phi Correction Starts...<<<<<<<<<" << std::endl;
    StRefMultCorr refmultCorrUtil = StRefMultCorr("refmult");

        for(Long_t i = 0; i < m_EventNumberToProcess; i++){
        if(i % 10000 == 0) std::cout << i << " Events Processed" <<  std::endl;
        m_Chain->GetEntry(i);
        StEvtInfo evtInfo(m_Chain, refmultCorrUtil);
        if(!m_EvtCuts->PassAllCuts(evtInfo)) continue;

	//((TH1D*)(m_HistsCollection.find("h1d_Bz")->second))->Fill(evtInfo.Bz());
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
Int_t StCorrelationMaker::initiateHistograms(){
    m_ResultFile->cd();
    //Below is only for Run11 AuAu200GeV Collisions
    TH1D* h1d_EventTally = new TH1D("h1d_EventTally", "Event Tally", 10, 0, 1); 
    h1d_EventTally->SetBit(TH1::kCanRebin);
    h1d_EventTally->SetStats(0);
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_EventTally", h1d_EventTally));

    //TH1D* h1d_Bz = new TH1D("h1d_Bz", "Magnetic Field", 10, -10, 10);
    TH1D* h1d_TriggerId = new TH1D("h1d_TriggerId", "Trigger Id", 200, 0.5, 200.5);//TODO:??
    TH1D* h1d_CentralityUnWeighted = new TH1D("h1d_CentralityUnWeighted", "Centrality UnWeighted", 10, 0, 10);
    TH1D* h1d_CentralityWeighted = new TH1D("h1d_CentralityWeighted", "Centrality Weighted", 10, 0, 10);
    TH1D* h1d_PVertexZUnWeighted = new TH1D("h1d_PVertexZUnWeighted", "PVertex Z UnWeighted", 100, -100, 100);
    TH1D* h1d_PVertexZWeighted = new TH1D("h1d_PVertexZWeighted", "PVertex Z Weighted", 100, -100, 100);
    TH2D* h2d_Mult_vs_PVertexZUnWeighted = new TH2D("h2d_Mult_vs_PVertexZUnWeighted", "Mult. vs PVertex UnWeighted", 1000, -.5, 999.5, 100, -100, 100);
    TH2D* h2d_Mult_vs_PVertexZWeighted = new TH2D("h2d_Mult_vs_PVertexZWeighted", "Mult. vs PVertex Weighted", 1000, -.5, 999.5, 100, -100, 100);
    //m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_Bz", h1d_Bz));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_TriggerId", h1d_TriggerId));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_CentralityUnWeighted", h1d_CentralityUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_CentralityWeighted", h1d_CentralityWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_PVertexZUnWeighted", h1d_PVertexZUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_PVertexZWeighted", h1d_PVertexZWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_Mult_vs_PVertexZUnWeighted", h2d_Mult_vs_PVertexZUnWeighted));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h2d_Mult_vs_PVertexZWeighted", h2d_Mult_vs_PVertexZWeighted));

    // Histograms used for phi correction for all primary tracks TODO:
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

    // Primary tracks before corrections
    TH1D* h1d_before_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi", "Primary Tracks Phi (FF-PVZPos-ChPos)", phiBins, -PI, PI); 
    TH1D* h1d_before_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi", "Primary Tracks Phi (FF-PVZPos-ChNeg)", phiBins, -PI, PI); 
    TH1D* h1d_before_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi", "Primary Tracks Phi (FF-PVZNeg-ChPos)", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi", "Primary Tracks Phi (FF-PVZNeg-ChNeg)", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi", "Primary Tracks Phi (RF-PVZPos-ChPos)", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi", "Primary Tracks Phi (RF-PVZPos-ChNeg)", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi", "Primary Tracks Phi (RF-PVZNeg-ChPos)", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi = new TH1D("h1d_before_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi", "Primary Tracks Phi (RF-PVZNeg-ChNeg)", phiBins, -PI, PI);
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi", h1d_before_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi", h1d_before_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi", h1d_before_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi", h1d_before_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi", h1d_before_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi", h1d_before_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi", h1d_before_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi", h1d_before_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi));

    TH1D* h1d_after_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi", "h1d_after_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi", "h1d_after_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi", "h1d_after_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi", "h1d_after_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi", "h1d_after_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi", "h1d_after_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi", "h1d_after_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi = new TH1D("h1d_after_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi", "h1d_after_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi", phiBins, -PI, PI);
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi", h1d_after_Corrections_FF_PVZPos_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi", h1d_after_Corrections_FF_PVZPos_ChNeg_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi", h1d_after_Corrections_FF_PVZNeg_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi", h1d_after_Corrections_FF_PVZNeg_ChNeg_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi", h1d_after_Corrections_RF_PVZPos_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi", h1d_after_Corrections_RF_PVZPos_ChNeg_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi", h1d_after_Corrections_RF_PVZNeg_ChPos_PrimaryTrkPhi));
    m_HistsCollection.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi", h1d_after_Corrections_RF_PVZNeg_ChNeg_PrimaryTrkPhi));

    // Tracks used to reconstruct event plane
    TH1D* h1d_EastWestCountDiff = new TH1D("h1d_EastWestCountDiff", "East_Count - West_Count", 500, -250, 250);
    TH1D* h1d_EastWestCountSum = new TH1D("h1d_EastWestCountSum", "East_Count + West_Count", 1000, 0, 1000);

    // TODO: 
    addAlphaHists(m_HistsCollection);// Phi distribution + invariant mass
    addBetaHists(m_HistsCollection);
    
    // TODO: add observable plots
    addCorrelationHists(m_HistsCollection);
    //TProfile* prof_GammaUnWeighted = new TProfile("prof_GammaUnWeighted", "prof_GammaUnWeighted", 2, .5, 2.5, -100, 100, "");
    //TProfile* prof_GammaWeighted = new TProfile("prof_GammaWeighted", "prof_GammaWeighted", 2, .5, 2.5, -100, 100, "");

    return 1;
}

Int_t StCorrelationMaker::addTrees(){
    Int_t nFiles = 0;
    for(set<string>::iterator iter = m_DataSet.begin(); iter != m_DataSet.end(); iter++)
        nFiles += m_Chain->Add((*iter).c_str());
    return nFiles;
    return 1; //TODO
}

void StCorrelationMaker::addAlphaHists(std::map<std::string, TH1*>& histMap){
    std::cout << histMap.size() << "size before addAlphaHists" << std::endl;
    // Used to check distribution before correction 
    TH1D* h1d_before_Corrections_FF_PVZPos_ChPos_AlphaPhi = new TH1D("h1d_before_Corrections_FF_PVZPos_ChPos_AlphaPhi", "h1d_before_Corrections_FF_PVZPos_ChPos_AlphaPhi", phiBins, -PI, PI);
    //std::cout << h1d_before_Corrections_FF_PVZPos_ChPos_AlphaPhi->GetTitle() << std::endl;
    TH1D* h1d_before_Corrections_FF_PVZPos_ChNeg_AlphaPhi = new TH1D("h1d_before_Corrections_FF_PVZPos_ChNeg_AlphaPhi", "h1d_before_Corrections_FF_PVZPos_ChNeg_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_FF_PVZNeg_ChPos_AlphaPhi = new TH1D("h1d_before_Corrections_FF_PVZNeg_ChPos_AlphaPhi", "h1d_before_Corrections_FF_PVZNeg_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_FF_PVZNeg_ChNeg_AlphaPhi = new TH1D("h1d_before_Corrections_FF_PVZNeg_ChNeg_AlphaPhi", "h1d_before_Corrections_FF_PVZNeg_ChNeg_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZPos_ChPos_AlphaPhi = new TH1D("h1d_before_Corrections_RF_PVZPos_ChPos_AlphaPhi", "h1d_before_Corrections_RF_PVZPos_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZPos_ChNeg_AlphaPhi = new TH1D("h1d_before_Corrections_RF_PVZPos_ChNeg_AlphaPhi", "h1d_before_Corrections_RF_PVZPos_ChNeg_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZNeg_ChPos_AlphaPhi = new TH1D("h1d_before_Corrections_RF_PVZNeg_ChPos_AlphaPhi", "h1d_before_Corrections_RF_PVZNeg_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZNeg_ChNeg_AlphaPhi = new TH1D("h1d_before_Corrections_RF_PVZNeg_ChNeg_AlphaPhi", "h1d_before_Corrections_RF_PVZNeg_ChNeg_AlphaPhi", phiBins, -PI, PI);

    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZPos_ChPos_AlphaPhi", h1d_before_Corrections_FF_PVZPos_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZPos_ChNeg_AlphaPhi", h1d_before_Corrections_FF_PVZPos_ChNeg_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZNeg_ChPos_AlphaPhi", h1d_before_Corrections_FF_PVZNeg_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZNeg_ChNeg_AlphaPhi", h1d_before_Corrections_FF_PVZNeg_ChNeg_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZPos_ChPos_AlphaPhi", h1d_before_Corrections_RF_PVZPos_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZPos_ChNeg_AlphaPhi", h1d_before_Corrections_RF_PVZPos_ChNeg_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZNeg_ChPos_AlphaPhi", h1d_before_Corrections_RF_PVZNeg_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZNeg_ChNeg_AlphaPhi", h1d_before_Corrections_RF_PVZNeg_ChNeg_AlphaPhi));

    // Used to fill correction terms to flatten the distribution of proton's phi
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Alpha", "prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Alpha", "prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Alpha", "prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Alpha", "prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Alpha", "prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Alpha", "prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Alpha", "prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Alpha = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Alpha", "prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Alpha", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");

    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Alpha", prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Alpha", prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Alpha", prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Alpha", prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Alpha", prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Alpha", prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Alpha", prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Alpha));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Alpha", prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Alpha));

    // Used to check after correction distribution
    TH1F* h1f_before_Corrections_AlphaPhi = new TH1F("h1f_before_Corrections_AlphaPhi", "h1f_before_Corrections_AlphaPhi", phiBins, -PI, PI);
    TH1F* h1f_after_Corrections_AlphaPhi = new TH1F("h1f_after_Corrections_AlphaPhi", "h1f_after_Corrections_AlphaPhi", phiBins, -PI, PI);
    histMap.insert(std::pair<std::string, TH1*>("h1f_before_Corrections_AlphaPhi", h1f_before_Corrections_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1f_after_Corrections_AlphaPhi", h1f_after_Corrections_AlphaPhi));

    TH1D* h1d_after_Corrections_FF_PVZPos_ChPos_AlphaPhi = new TH1D("h1d_after_Corrections_FF_PVZPos_ChPos_AlphaPhi", "h1d_after_Corrections_FF_PVZPos_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZPos_ChNeg_AlphaPhi = new TH1D("h1d_after_Corrections_FF_PVZPos_ChNeg_AlphaPhi", "h1d_after_Corrections_FF_PVZPos_ChNeg_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZNeg_ChPos_AlphaPhi = new TH1D("h1d_after_Corrections_FF_PVZNeg_ChPos_AlphaPhi", "h1d_after_Corrections_FF_PVZNeg_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZNeg_ChNeg_AlphaPhi = new TH1D("h1d_after_Corrections_FF_PVZNeg_ChNeg_AlphaPhi", "h1d_after_Corrections_FF_PVZNeg_ChNeg_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZPos_ChPos_AlphaPhi = new TH1D("h1d_after_Corrections_RF_PVZPos_ChPos_AlphaPhi", "h1d_after_Corrections_RF_PVZPos_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZPos_ChNeg_AlphaPhi = new TH1D("h1d_after_Corrections_RF_PVZPos_ChNeg_AlphaPhi", "h1d_after_Corrections_RF_PVZPos_ChNeg_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZNeg_ChPos_AlphaPhi = new TH1D("h1d_after_Corrections_RF_PVZNeg_ChPos_AlphaPhi", "h1d_after_Corrections_RF_PVZNeg_ChPos_AlphaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZNeg_ChNeg_AlphaPhi = new TH1D("h1d_after_Corrections_RF_PVZNeg_ChNeg_AlphaPhi", "h1d_after_Corrections_RF_PVZNeg_ChNeg_AlphaPhi", phiBins, -PI, PI);

    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZPos_ChPos_AlphaPhi", h1d_after_Corrections_FF_PVZPos_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZPos_ChNeg_AlphaPhi", h1d_after_Corrections_FF_PVZPos_ChNeg_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZNeg_ChPos_AlphaPhi", h1d_after_Corrections_FF_PVZNeg_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZNeg_ChNeg_AlphaPhi", h1d_after_Corrections_FF_PVZNeg_ChNeg_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZPos_ChPos_AlphaPhi", h1d_after_Corrections_RF_PVZPos_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZPos_ChNeg_AlphaPhi", h1d_after_Corrections_RF_PVZPos_ChNeg_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZNeg_ChPos_AlphaPhi", h1d_after_Corrections_RF_PVZNeg_ChPos_AlphaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZNeg_ChNeg_AlphaPhi", h1d_after_Corrections_RF_PVZNeg_ChNeg_AlphaPhi));

    std::cout << histMap.size() << "size after addAlphaHists" << std::endl;

    addAdditionalAlphaHists(histMap);
    std::cout << histMap.size() << "size after addAdditionalAlphaHists" << std::endl;
}

void StCorrelationMaker::addBetaHists(std::map<std::string, TH1*>& histMap){
    TH1D* h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi = new TH1D("h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi", "h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi", phiBins, -PI, PI);
    //std::cout << h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi->GetTitle() << std::endl;
    TH1D* h1d_before_Corrections_FF_PVZPos_ChNeg_BetaPhi = new TH1D("h1d_before_Corrections_FF_PVZPos_ChNeg_BetaPhi", "h1d_before_Corrections_FF_PVZPos_ChNeg_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_FF_PVZNeg_ChPos_BetaPhi = new TH1D("h1d_before_Corrections_FF_PVZNeg_ChPos_BetaPhi", "h1d_before_Corrections_FF_PVZNeg_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_FF_PVZNeg_ChNeg_BetaPhi = new TH1D("h1d_before_Corrections_FF_PVZNeg_ChNeg_BetaPhi", "h1d_before_Corrections_FF_PVZNeg_ChNeg_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZPos_ChPos_BetaPhi = new TH1D("h1d_before_Corrections_RF_PVZPos_ChPos_BetaPhi", "h1d_before_Corrections_RF_PVZPos_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZPos_ChNeg_BetaPhi = new TH1D("h1d_before_Corrections_RF_PVZPos_ChNeg_BetaPhi", "h1d_before_Corrections_RF_PVZPos_ChNeg_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZNeg_ChPos_BetaPhi = new TH1D("h1d_before_Corrections_RF_PVZNeg_ChPos_BetaPhi", "h1d_before_Corrections_RF_PVZNeg_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_before_Corrections_RF_PVZNeg_ChNeg_BetaPhi = new TH1D("h1d_before_Corrections_RF_PVZNeg_ChNeg_BetaPhi", "h1d_before_Corrections_RF_PVZNeg_ChNeg_BetaPhi", phiBins, -PI, PI);

    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi", h1d_before_Corrections_FF_PVZPos_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZPos_ChNeg_BetaPhi", h1d_before_Corrections_FF_PVZPos_ChNeg_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZNeg_ChPos_BetaPhi", h1d_before_Corrections_FF_PVZNeg_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_FF_PVZNeg_ChNeg_BetaPhi", h1d_before_Corrections_FF_PVZNeg_ChNeg_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZPos_ChPos_BetaPhi", h1d_before_Corrections_RF_PVZPos_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZPos_ChNeg_BetaPhi", h1d_before_Corrections_RF_PVZPos_ChNeg_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZNeg_ChPos_BetaPhi", h1d_before_Corrections_RF_PVZNeg_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_before_Corrections_RF_PVZNeg_ChNeg_BetaPhi", h1d_before_Corrections_RF_PVZNeg_ChNeg_BetaPhi));

    // Used to fill correction terms to flatten the distribution of Lambda's phi
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Beta", "prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Beta", "prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Beta", "prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Beta", "prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Beta", "prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Beta", "prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Beta", "prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");
    TProfile2D* prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Beta = new TProfile2D("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Beta", "prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Beta", 2*order, .5, 2*order + .5, m_Day2UpperBound - m_Day2LowerBound, m_Day2LowerBound, m_Day2UpperBound, -1, 1, "");

    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Beta", prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Beta", prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChNeg_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Beta", prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChPos_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Beta", prof2_XOrder_YDay_ZPhi_FF_PVZNeg_ChNeg_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Beta", prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChPos_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Beta", prof2_XOrder_YDay_ZPhi_RF_PVZPos_ChNeg_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Beta", prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChPos_Beta));
    histMap.insert(std::pair<std::string, TH1*>("prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Beta", prof2_XOrder_YDay_ZPhi_RF_PVZNeg_ChNeg_Beta));
   
    // Used to check after correction distribution
    TH1F* h1f_before_Corrections_BetaPhi = new TH1F("h1f_before_Corrections_BetaPhi", "h1f_before_Corrections_BetaPhi", phiBins, -PI, PI);
    TH1F* h1f_after_Corrections_BetaPhi = new TH1F("h1f_after_Corrections_BetaPhi", "h1f_after_Corrections_BetaPhi", phiBins, -PI, PI);
    histMap.insert(std::pair<std::string, TH1*>("h1f_before_Corrections_BetaPhi", h1f_before_Corrections_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1f_after_Corrections_BetaPhi", h1f_after_Corrections_BetaPhi));

    TH1D* h1d_after_Corrections_FF_PVZPos_ChPos_BetaPhi = new TH1D("h1d_after_Corrections_FF_PVZPos_ChPos_BetaPhi", "h1d_after_Corrections_FF_PVZPos_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZPos_ChNeg_BetaPhi = new TH1D("h1d_after_Corrections_FF_PVZPos_ChNeg_BetaPhi", "h1d_after_Corrections_FF_PVZPos_ChNeg_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZNeg_ChPos_BetaPhi = new TH1D("h1d_after_Corrections_FF_PVZNeg_ChPos_BetaPhi", "h1d_after_Corrections_FF_PVZNeg_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_FF_PVZNeg_ChNeg_BetaPhi = new TH1D("h1d_after_Corrections_FF_PVZNeg_ChNeg_BetaPhi", "h1d_after_Corrections_FF_PVZNeg_ChNeg_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZPos_ChPos_BetaPhi = new TH1D("h1d_after_Corrections_RF_PVZPos_ChPos_BetaPhi", "h1d_after_Corrections_RF_PVZPos_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZPos_ChNeg_BetaPhi = new TH1D("h1d_after_Corrections_RF_PVZPos_ChNeg_BetaPhi", "h1d_after_Corrections_RF_PVZPos_ChNeg_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZNeg_ChPos_BetaPhi = new TH1D("h1d_after_Corrections_RF_PVZNeg_ChPos_BetaPhi", "h1d_after_Corrections_RF_PVZNeg_ChPos_BetaPhi", phiBins, -PI, PI);
    TH1D* h1d_after_Corrections_RF_PVZNeg_ChNeg_BetaPhi = new TH1D("h1d_after_Corrections_RF_PVZNeg_ChNeg_BetaPhi", "h1d_after_Corrections_RF_PVZNeg_ChNeg_BetaPhi", phiBins, -PI, PI);

    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZPos_ChPos_BetaPhi", h1d_after_Corrections_FF_PVZPos_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZPos_ChNeg_BetaPhi", h1d_after_Corrections_FF_PVZPos_ChNeg_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZNeg_ChPos_BetaPhi", h1d_after_Corrections_FF_PVZNeg_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_FF_PVZNeg_ChNeg_BetaPhi", h1d_after_Corrections_FF_PVZNeg_ChNeg_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZPos_ChPos_BetaPhi", h1d_after_Corrections_RF_PVZPos_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZPos_ChNeg_BetaPhi", h1d_after_Corrections_RF_PVZPos_ChNeg_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZNeg_ChPos_BetaPhi", h1d_after_Corrections_RF_PVZNeg_ChPos_BetaPhi));
    histMap.insert(std::pair<std::string, TH1*>("h1d_after_Corrections_RF_PVZNeg_ChNeg_BetaPhi", h1d_after_Corrections_RF_PVZNeg_ChNeg_BetaPhi));

    addAdditionalBetaHists(histMap);
}

void StCorrelationMaker::zeroWeightsHistograms(){
    char particles[3] = {'a', 'b', 'p'}; 
    for(Int_t i = 0; i < 2; i++){
        for(Int_t ii = 0; ii < 2; ii++){
            for(Int_t iii = 0; iii < 2; iii++){
		std::string Eta = ((i == 0)? "FF" : "RF");
		std::string PVZ = ((ii == 0)? "PVZPos" : "PVZNeg");
		std::string Charge = ((iii == 0)? "ChPos" : "ChNeg");

		char histname[100];
   
                for(Int_t j = 0; j < 3; j++){ 
                    char particle = particles[j];
		    if(particle == 'a' || particle == 'b'){
			std::string particle_name = (particle == 'a')? "Alpha" : "Beta";
			sprintf(histname, "h1d_before_Corrections_%s_%s_%s_%sPhi", Eta.c_str(), PVZ.c_str(), Charge.c_str(), particle_name.c_str());
		    }
		    else
			sprintf(histname, "h1d_before_Corrections_%s_%s_%s_PrimaryTrkPhi", Eta.c_str(), PVZ.c_str(), Charge.c_str());

		    std::string histname_string(histname);
		    //std::cout << histname_string << std::endl;
		    TH1D* hist = (TH1D*)m_HistsCollection[histname_string];

		    for(Int_t iiii = 0; iiii < phiBins; iiii++)
                        hist->SetBinContent(iiii + 1, 0.);
		}
	    }
	}
    }
}
