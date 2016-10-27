#ifndef STPROTONLAGAMMAMAKER_H
#define STPROTONLAGAMMAMAKER_H
#include <iostream>
#include <string>
#include <cstdio>
#include "TLeaf.h"
#include "StCorrelationMaker.h"
#include "StEffMaker.h"
#include "../StEvtInfo/StEvtInfo.h"
#include "../StPriTrkCuts/StPriTrkCuts.h"
#include "../StV0TrkCuts/StV0TrkCuts.h"
#include <map>
class StProtonLaGammaMaker: public StCorrelationMaker{
    public:
        StProtonLaGammaMaker(Int_t centrality, Int_t eventnumbertoprocess, StEffMaker* emaker1, StEffMaker* emaker2, StEffMaker* emaker3, StEffMaker* emaker4): StCorrelationMaker(centrality, eventnumbertoprocess, emaker1, emaker2, emaker3, emaker4){}
	void LoadPriTrkCuts(StPriTrkCuts* cuts){ m_PriTrkCuts = cuts;}
	void LoadProtonTrkCuts(StPriTrkCuts* cuts){ m_ProtonTrkCuts = cuts; }
	void LoadLambdaTrkCuts(StV0TrkCuts* cuts){ m_LambdaTrkCuts = cuts;}
        virtual ~StProtonLaGammaMaker(){}

    private:
        StPriTrkCuts* m_PriTrkCuts;
	StPriTrkCuts* m_ProtonTrkCuts;
	StV0TrkCuts* m_LambdaTrkCuts;

        Float_t m_PrimaryTracksPhiWeight[2][2][2][phiBins]; // Eta, PVZ, Charge, phiBins
        Float_t m_AlphaPhiWeight[2][2][2][phiBins];
        Float_t m_BetaPhiWeight[2][2][2][phiBins];

	virtual void fillPrimaryTracksPhiHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);
	virtual void fillAlphaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);
	virtual void fillBetaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);

	virtual std::string getCorrectionHistName(std::string particle, Float_t Eta, Float_t PVZ, Float_t Charge, char option) const;
	virtual void addAdditionalAlphaHists(std::map<std::string, TH1*>& histMap);
	virtual void addAdditionalBetaHists(std::map<std::string, TH1*>& histMap);
        virtual void addCorrelationHists(std::map<std::string, TH1*>& histMap); 

        virtual void computePhiWeightsHelper(Int_t idx1, Int_t idx2, Int_t idx3, Char_t particle, std::map<std::string, TH1*>& histMap);

        virtual void reconstructSubEventPlaneWithPhiWeightHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo);
        virtual void reconstructShiftedSubEventPlaneHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo);
        virtual void reconstructShiftedFullEventPlaneHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo, TVector2& shiftedepphi_full);
        virtual void computeCorrelatorsHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo, TVector2 shiftedepphi_full);

};

inline
std::string StProtonLaGammaMaker::getCorrectionHistName(std::string particle, Float_t Eta, Float_t PVZ, Float_t Charge, char option) const{
    char name[100];
    std::string field = ((Eta > 0)? "FF":"RF");
    std::string vz = ((PVZ > 0)? "PVZPos":"PVZNeg");
    std::string charge = ((Charge > 0)? "ChPos":"ChNeg");
    std::string a;

    switch(option){
        // Before corrections
        case 'b': 
            sprintf(name, "h1d_before_Corrections_%s_%s_%s_%sPhi", field.c_str(), vz.c_str(), charge.c_str(), particle.c_str());
            a = name;    
            break;
        // Correction terms
        case 'c':
	    sprintf(name, "prof2_XOrder_YDay_ZPhi_%s_%s_%s_%s", field.c_str(), vz.c_str(), charge.c_str(), particle.c_str());
	    //sprintf(name, "prof2_XOrder_YDay_ZPhi_FF_PVZPos_ChPos_%s", particle.c_str());
	    a = name;
            break;
        // After corrections
        case 'a': 
            sprintf(name, "h1d_after_Corrections_%s_%s_%s_%sPhi", field.c_str(), vz.c_str(), charge.c_str(), particle.c_str());
            a = name;
            break;
        default:
            std::cout << "bad option for getCorrectionHistName!!" << std::endl;
            exit(1);
            break;
    }
    return a;
}
#endif
