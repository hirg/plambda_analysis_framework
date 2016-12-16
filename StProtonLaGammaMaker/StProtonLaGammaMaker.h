#ifndef STPROTONLAGAMMAMAKER_H
#define STPROTONLAGAMMAMAKER_H
#include <iostream>
#include <string>
#include <cstdio>
#include "TLeaf.h"
#include "../StCorrelationMaker/StCorrelationMaker.h"
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

	virtual void addAdditionalAlphaHists(std::map<std::string, TH1*>& histMap);
	virtual void addAdditionalBetaHists(std::map<std::string, TH1*>& histMap);
        virtual void addCorrelationHists(std::map<std::string, TH1*>& histMap); 

	virtual void fillPrimaryTracksPhiHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);
	virtual void fillAlphaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);
	virtual void fillBetaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);

	virtual std::string getCorrectionHistName(std::string particle, Float_t Eta, Float_t PVZ, Float_t Charge, char option) const;
	
        virtual void computePhiWeightsHelper(Int_t, Int_t, Int_t, Char_t particle, std::map<std::string, TH1*>&);
        virtual void reconstructSubEventPlaneWithPhiWeightHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);
        virtual void reconstructShiftedSubEventPlaneHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo);
        virtual void reconstructShiftedFullEventPlaneHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo, TVector2& shiftedepphi_full);
        virtual void computeCorrelatorsHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo, TVector2 shiftedepphi_full);

        virtual void saveWeights0File(std::map<std::string, TH1*>&);
        virtual void saveWeights1File(std::map<std::string, TH1*>&){}
        virtual void saveWeights2File(std::map<std::string, TH1*>&){}
};

inline 
void StProtonLaGammaMaker::saveWeights0File(std::map<std::string, TH1*>& histmap){
    TFile* weights0File = GetWeights0File();
    weights0File->cd();

    std::string eta_space[2] = {"FF", "RF"};
    std::string vz_space[2] = {"PVZPos", "PVZNeg"};
    std::string charge_space[2] = {"ChPos", "ChNeg"};
    std::string particle_space[3] = {"Alpha", "Beta", "PrimaryTrk"};
    for(int i = 0; i != 2; ++i){
        for(int j = 0; j != 2; ++j){
	    for(int k = 0; k != 2; ++k){
		for(int l = 0; l != 3; ++l){
		    char name[100];
		    sprintf(name, "h1d_before_Corrections_%s_%s_%s_%sPhi", eta_space[i].c_str(), vz_space[j].c_str(), charge_space[k].c_str(), particle_space[l].c_str());
		    std::string name_str(name);
		    histmap[name_str]->Write();
		}
	    }
	}
    }
}

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
        // Correction terms, this is only for correlated particles(Alpha and Beta)
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

// This one is called in StCorrelationMaker.cpp
inline
void StProtonLaGammaMaker::computePhiWeightsHelper(Int_t i, Int_t ii, Int_t iii, Char_t particle, std::map<std::string, TH1*>& histMap){
    std::string Eta = ((i == 0)? "FF" : "RF");
    std::string PVZ = ((ii == 0)? "PVZPos" : "PVZNeg");
    std::string Charge = ((iii == 0)? "ChPos" : "ChNeg");

    char histname[100];
   
    //std::cout << "particle is " << particle << std::endl;
    std::map<char, std::string> map_par;
    map_par.insert(std::pair<char, std::string>('a', "Alpha"));
    map_par.insert(std::pair<char, std::string>('b', "Beta"));
    map_par.insert(std::pair<char, std::string>('p', "PrimaryTrk"));
    sprintf(histname, "h1d_before_Corrections_%s_%s_%s_%sPhi", Eta.c_str(), PVZ.c_str(), Charge.c_str(), map_par[particle].c_str());

    std::string histname_string(histname);
    //std::cout << histname_string << std::endl;
    TH1D* hist = (TH1D*)histMap[histname_string];
    Float_t phi_mean = ((TH1D*)histMap[histname_string])->GetSum() / (Float_t)phiBins;

    switch(particle){
	case 'a':
	    for(Int_t j = 0; j < phiBins; ++j)
		m_AlphaPhiWeight[i][ii][iii][j] = (hist->GetBinContent(j + 1) != 0.)? phi_mean / (hist->GetBinContent(j + 1)) : 1;
	    break;
	case 'b':
	    for(Int_t j = 0; j < phiBins; ++j)
		m_BetaPhiWeight[i][ii][iii][j] = (hist->GetBinContent(j + 1) != 0.)?  phi_mean / (hist->GetBinContent(j + 1)) : 1;
	    break;
	case 'p':
	    for(Int_t j = 0; j < phiBins; ++j)
		m_PrimaryTracksPhiWeight[i][ii][iii][j] = (hist->GetBinContent(j + 1) != 0.)? phi_mean / (hist->GetBinContent(j + 1)) : 1;
	    break;
	default:
	    break;
    }  
    return;
}
#endif
