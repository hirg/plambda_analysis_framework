#include <utility>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath> 
#include "StProtonLaGammaMaker.h"
#include "../include/constants.h"
#include "../StEvtInfo/StEvtInfo.h"
#include "TChain.h"
#include "TLeaf.h"
#include <map>
#include <string>
#include "TVector2.h"
#include "TH2F.h"
#include "TProfile2D.h" 
#include "TProfile.h"
using namespace std;
void StProtonLaGammaMaker::addAdditionalAlphaHists(std::map<std::string, TH1*>& histMap){
    TH2F* h2f_XEta_YPt_ProtonUnWeighted = new TH2F("h2f_XEta_YPt_ProtonUnWeighted", "h2f_XEta_YPt_ProtonUnWeighted", 26, -1.3, 1.3, 300, 0, 15);
    TH2F* h2f_XEta_YPt_ProtonWeighted = new TH2F("h2f_XEta_YPt_ProtonWeighted", "h2f_XEta_YPt_ProtonWeighted", 26, -1.3, 1.3, 300, 0, 15);
    histMap.insert(std::pair<std::string, TH1*>("h2f_XEta_YPt_ProtonUnWeighted", h2f_XEta_YPt_ProtonUnWeighted));
    histMap.insert(std::pair<std::string, TH1*>("h2f_XEta_YPt_ProtonWeighted", h2f_XEta_YPt_ProtonWeighted));

    // Alpha v2
    TProfile* profile_XPt_Yv2_ProtonUnWeighted = new TProfile("profile_XPt_Yv2_ProtonUnWeighted", "v2 of Protons", 60, 0, 15, -100, 100, ""); 
    TProfile* profile_XPt_Yv2_ProtonWeighted = new TProfile("profile_XPt_Yv2_ProtonWeighted", "v2 of Protons with EWeight", 60, 0, 15, -100, 100, ""); 
    histMap.insert(std::pair<std::string, TH1*>("profile_XPt_Yv2_ProtonUnWeighted", profile_XPt_Yv2_ProtonUnWeighted));
    histMap.insert(std::pair<std::string, TH1*>("profile_XPt_Yv2_ProtonWeighted", profile_XPt_Yv2_ProtonWeighted));
}

void StProtonLaGammaMaker::addAdditionalBetaHists(std::map<std::string, TH1*>& histMap){
    const Float_t masswidth = 0.08;

    TH2F* h2f_XEta_YPt_LambdaUnWeighted = new TH2F("h2f_XEta_YPt_LambdaUnWeighted", "h2f_XEta_YPt_LambdaUnWeighted", 26, -1.3, 1.3, 300, 0, 15);
    TH2F* h2f_XEta_YPt_LambdaWeighted = new TH2F("h2f_XEta_YPt_LambdaWeighted", "h2f_XEta_YPt_LambdaWeighted", 26, -1.3, 1.3, 300, 0, 15);
    histMap.insert(std::pair<std::string, TH1*>("h2f_XEta_YPt_LambdaUnWeighted", h2f_XEta_YPt_LambdaUnWeighted));
    histMap.insert(std::pair<std::string, TH1*>("h2f_XEta_YPt_LambdaWeighted", h2f_XEta_YPt_LambdaWeighted));

    TH1F* h1f_LambdaMassUnWeighted = new TH1F("h1f_LambdaMassUnWeighted","Invariant Mass after cuts Uweighted", 200, LambdaPDGMass - masswidth, LambdaPDGMass + masswidth);
    TH1F* h1f_LambdaMassWeighted = new TH1F("h1f_LambdaMassWeighted","Invariant Mass after cuts Weighted", 200, LambdaPDGMass - masswidth, LambdaPDGMass + masswidth);
    histMap.insert(std::pair<std::string, TH1*>("h1f_LambdaMassUnWeighted", h1f_LambdaMassUnWeighted));
    histMap.insert(std::pair<std::string, TH1*>("h1f_LambdaMassWeighted", h1f_LambdaMassWeighted));

    // Beta v2
    TProfile* profile_XPt_Yv2_LambdaUnWeighted = new TProfile("profile_XPt_Yv2_LambdaUnWeighted", "v2 of Lambdas", 60, 0, 15, -100, 100, ""); 
    TProfile* profile_XPt_Yv2_LambdaWeighted = new TProfile("profile_XPt_Yv2_LambdaWeighted", "v2 of Lambdas w/o EWeight", 60, 0, 15, -100, 100, ""); 
    histMap.insert(std::pair<std::string, TH1*>("profile_XPt_Yv2_LambdaUnWeighted", profile_XPt_Yv2_LambdaUnWeighted));
    histMap.insert(std::pair<std::string, TH1*>("profile_XPt_Yv2_LambdaWeighted", profile_XPt_Yv2_LambdaWeighted));
}

void StProtonLaGammaMaker::addCorrelationHists(std::map<std::string, TH1*>& histMap){
    // Gamma correlator vs runId
    TProfile* profile_XDay2_YSameSignGamma = new TProfile("profile_XDay2_YSameSignGamma", "#gamma vs Day2(ss)", Day2UpperBound() - Day2LowerBound(), Day2LowerBound(), Day2UpperBound(), -100, 100, "");
    TProfile* profile_XDay2_YOppoSignGamma = new TProfile("profile_XDay2_YOppoSignGamma", "#gamma vs Day2(os)", Day2UpperBound() - Day2LowerBound(), Day2LowerBound(), Day2UpperBound(), -100, 100, "");
    histMap.insert(std::pair<std::string, TH1*>("profile_XDay2_YSameSignGamma", profile_XDay2_YSameSignGamma));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDay2_YOppoSignGamma", profile_XDay2_YOppoSignGamma));

    // Gamma correlator 
    TProfile* profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma = new TProfile("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma", "#gamma Lambda-Proton", 4, 0.5, 4.5, -100, 100, "");
    TProfile* profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_noEW = new TProfile("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_noEW", "#gamma Lambad-Proton w/ event weight", 4, 0.5, 4.5, -100, 100, "");  
    TProfile* profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff = new TProfile("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff", "#gamma Lambda-Proton", 4, 0.5, 4.5, -100, 100, "");
    histMap.insert(std::pair<std::string, TH1*>("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma", profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma));
    histMap.insert(std::pair<std::string, TH1*>("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff", profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_noEW", profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_noEW));

    // MeanPt_vs_gamma correlator
    TProfile* profile_XMeanPt_YGamma_LambdaP = new TProfile("profile_XMeanPt_YGamma_LambdaP", "#gamma v.s. #bar{p_t} for Lambda-P", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YGamma_LambdaP_Eff = new TProfile("profile_XMeanPt_YGamma_LambdaP_Eff", "#gamma v.s. #bar{p_t} for Lambda-P with eff. correction", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YGamma_LambdaPBar = new TProfile("profile_XMeanPt_YGamma_LambdaPBar", "#gamma v.s. #bar{p_t} for Lambda-PBar", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YGamma_LambdaPBar_Eff = new TProfile("profile_XMeanPt_YGamma_LambdaPBar_Eff", "#gamma v.s. #bar{p_t} for Lambda-PBar with eff. correction", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YGamma_LambdaBarP = new TProfile("profile_XMeanPt_YGamma_LambdaBarP", "#gamma v.s. #bar{p_t} for LambdaBar-P", 35, 0, 3.5, -100., 100.);
    TProfile* profile_XMeanPt_YGamma_LambdaBarP_Eff = new TProfile("profile_XMeanPt_YGamma_LambdaBarP_Eff", "#gamma v.s. #bar{p_t} for LambdaBar-P with eff. correction", 35, 0, 3.5, -100., 100.);
    TProfile* profile_XMeanPt_YGamma_LambdaBarPBar = new TProfile("profile_XMeanPt_YGamma_LambdaBarPBar", "#gamma v.s. #bar{p_t} for LambdaBar-PBar", 35, 0, 3.5, -100., 100.);
    TProfile* profile_XMeanPt_YGamma_LambdaBarPBar_Eff = new TProfile("profile_XMeanPt_YGamma_LambdaBarPBar_Eff", "#gamma v.s. #bar{p_t} for LambdaBar-PBar with eff. correction", 35, 0, 3.5, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaP", profile_XMeanPt_YGamma_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaP_Eff", profile_XMeanPt_YGamma_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaPBar", profile_XMeanPt_YGamma_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaPBar_Eff", profile_XMeanPt_YGamma_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaBarP", profile_XMeanPt_YGamma_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaBarP_Eff", profile_XMeanPt_YGamma_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaBarPBar", profile_XMeanPt_YGamma_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YGamma_LambdaBarPBar_Eff", profile_XMeanPt_YGamma_LambdaBarPBar_Eff));

    // DeltaPt_vs_gamma correlator
    TProfile* profile_XDeltaPt_YGamma_LambdaP = new TProfile("profile_XDeltaPt_YGamma_LambdaP", "#gamma v.s. #delta p_t for Lambda-P", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YGamma_LambdaP_Eff = new TProfile("profile_XDeltaPt_YGamma_LambdaP_Eff", "#gamma v.s. #delta p_t for Lambda-P with eff. correction", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YGamma_LambdaPBar = new TProfile("profile_XDeltaPt_YGamma_LambdaPBar", "#gamma v.s. #delta p_t for Lambda-PBar", 50, 0, 5.0, -100., 100.); 
    TProfile* profile_XDeltaPt_YGamma_LambdaPBar_Eff = new TProfile("profile_XDeltaPt_YGamma_LambdaPBar_Eff", "#gamma v.s. #delta p_t for Lambda-PBar with eff. correction", 50, 0, 5.0, -100., 100.);   
    TProfile* profile_XDeltaPt_YGamma_LambdaBarP = new TProfile("profile_XDeltaPt_YGamma_LambdaBarP", "#gamma v.s. #delta p_t for LambdaBar-P", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YGamma_LambdaBarP_Eff = new TProfile("profile_XDeltaPt_YGamma_LambdaBarP_Eff", "#gamma v.s. #delta p_t for LambdaBar-P with eff. correction", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YGamma_LambdaBarPBar = new TProfile("profile_XDeltaPt_YGamma_LambdaBarPBar", "#gamma v.s. #delta p_t for LambdaBar-PBar", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YGamma_LambdaBarPBar_Eff = new TProfile("profile_XDeltaPt_YGamma_LambdaBarPBar_Eff", "#gamma v.s. #delta p_t for Lambda-PBar with eff. correction", 50, 0, 5.0, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaP", profile_XDeltaPt_YGamma_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaP_Eff", profile_XDeltaPt_YGamma_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaPBar", profile_XDeltaPt_YGamma_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaPBar_Eff", profile_XDeltaPt_YGamma_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaBarP", profile_XDeltaPt_YGamma_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaBarP_Eff", profile_XDeltaPt_YGamma_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaBarPBar", profile_XDeltaPt_YGamma_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YGamma_LambdaBarPBar_Eff", profile_XDeltaPt_YGamma_LambdaBarPBar_Eff));

    // MeanEta_vs_gamma correlator
    TProfile* profile_XMeanEta_YGamma_LambdaP = new TProfile("profile_XMeanEta_YGamma_LambdaP", "#gamma v.s. #bar{#eta} for Lambda-P", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaP_Eff = new TProfile("profile_XMeanEta_YGamma_LambdaP_Eff", "#gamma v.s. #bar{#eta} for Lambda-P with eff. correction", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaPBar = new TProfile("profile_XMeanEta_YGamma_LambdaPBar", "#gamma v.s. #bar{#eta} for Lambda-PBar", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaPBar_Eff = new TProfile("profile_XMeanEta_YGamma_LambdaPBar_Eff", "#gamma v.s. #bar{#eta} for Lambda-PBar with eff. correction", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaBarP = new TProfile("profile_XMeanEta_YGamma_LambdaBarP", "#gamma v.s. #bar{#eta} for LambdaBar-P", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaBarP_Eff = new TProfile("profile_XMeanEta_YGamma_LambdaBarP_Eff", "#gamma v.s. #bar{#eta} for LambdaBar-P with eff. correction", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaBarPBar = new TProfile("profile_XMeanEta_YGamma_LambdaBarPBar", "#gamma v.s. #bar{#eta} for LambdaBar-PBar", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YGamma_LambdaBarPBar_Eff = new TProfile("profile_XMeanEta_YGamma_LambdaBarPBar_Eff", "#gamma v.s. #bar{#eta} for LambdaBar-PBar with eff. correction", 20, -1.0, 1.0, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaP", profile_XMeanEta_YGamma_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaP_Eff", profile_XMeanEta_YGamma_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaPBar", profile_XMeanEta_YGamma_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaPBar_Eff", profile_XMeanEta_YGamma_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaBarP", profile_XMeanEta_YGamma_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaBarP_Eff", profile_XMeanEta_YGamma_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaBarPBar", profile_XMeanEta_YGamma_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YGamma_LambdaBarPBar_Eff", profile_XMeanEta_YGamma_LambdaBarPBar_Eff));

    // DeltaEta_vs_gamma correlator
    TProfile* profile_XDeltaEta_YGamma_LambdaP = new TProfile("profile_XDeltaEta_YGamma_LambdaP", "#gamma v.s. #delta eta for Lambda-P", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YGamma_LambdaP_Eff = new TProfile("profile_XDeltaEta_YGamma_LambdaP_Eff", "#gamma v.s. #delta eta for Lambda-P with eff. correction", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YGamma_LambdaPBar = new TProfile("profile_XDeltaEta_YGamma_LambdaPBar", "#gamma v.s. #delta eta for Lambda-PBar", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YGamma_LambdaPBar_Eff = new TProfile("profile_XDeltaEta_YGamma_LambdaPBar_Eff", "#gamma v.s. #delta eta for Lambda-PBar with eff. correction", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YGamma_LambdaBarP = new TProfile("profile_XDeltaEta_YGamma_LambdaBarP", "#gamma v.s. #delta eta for LambdaBar-P", 20, 0, 2.0, -100., 100.);
    TProfile* profile_XDeltaEta_YGamma_LambdaBarP_Eff = new TProfile("profile_XDeltaEta_YGamma_LambdaBarP_Eff", "#gamma v.s. #delta eta for LambdaBar-P with eff. correction", 20, 0, 2.0, -100., 100.);
    TProfile* profile_XDeltaEta_YGamma_LambdaBarPBar = new TProfile("profile_XDeltaEta_YGamma_LambdaBarPBar", "#gamma v.s. #delta eta for LambdaBar-PBar", 20, 0, 2.0, -100., 100.);
    TProfile* profile_XDeltaEta_YGamma_LambdaBarPBar_Eff = new TProfile("profile_XDeltaEta_YGamma_LambdaBarPBar_Eff", "#gamma v.s. #delta eta for LambdaBar-PBar with eff. correction", 20, 0, 2.0, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaP", profile_XDeltaEta_YGamma_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaP_Eff", profile_XDeltaEta_YGamma_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaPBar", profile_XDeltaEta_YGamma_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaPBar_Eff", profile_XDeltaEta_YGamma_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaBarP", profile_XDeltaEta_YGamma_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaBarP_Eff", profile_XDeltaEta_YGamma_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaBarPBar", profile_XDeltaEta_YGamma_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YGamma_LambdaBarPBar_Eff", profile_XDeltaEta_YGamma_LambdaBarPBar_Eff));

    // Delta correlator
    TProfile* profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta = new TProfile("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta", "#delta Lambda-Proton", 4, 0.5, 4.5, -100., 100., "");
    TProfile* profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta_Eff = new TProfile("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta_Eff", "#delta Lambda-Proton with eff.", 4, 0.5, 4.5, -100., 100., "");
    histMap.insert(std::pair<std::string, TH1*>("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta", profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta));
    histMap.insert(std::pair<std::string, TH1*>("profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta_Eff", profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta_Eff));

    // Delta_vs_MeanPt correlator
    TProfile* profile_XMeanPt_YDelta_LambdaP = new TProfile("profile_XMeanPt_YDelta_LambdaP", "#delta v.s. #bar{p_t} for Lambda-P", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YDelta_LambdaP_Eff = new TProfile("profile_XMeanPt_YDelta_LambdaP_Eff", "#delta v.s. #bar{p_t} for Lambda-P with eff. correction", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YDelta_LambdaPBar = new TProfile("profile_XMeanPt_YDelta_LambdaPBar", "#delta v.s. #bar{p_t} for Lambda-PBar", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YDelta_LambdaPBar_Eff = new TProfile("profile_XMeanPt_YDelta_LambdaPBar_Eff", "#delta v.s. #bar{p_t} for Lambda-PBar with eff. correction", 35, 0, 3.5, -100., 100.);   
    TProfile* profile_XMeanPt_YDelta_LambdaBarP = new TProfile("profile_XMeanPt_YDelta_LambdaBarP", "#delta v.s. #bar{p_t} for LambdaBar-P", 35, 0, 3.5, -100., 100.);
    TProfile* profile_XMeanPt_YDelta_LambdaBarP_Eff = new TProfile("profile_XMeanPt_YDelta_LambdaBarP_Eff", "#delta v.s. #bar{p_t} for LambdaBar-P with eff. correction", 35, 0, 3.5, -100., 100.);
    TProfile* profile_XMeanPt_YDelta_LambdaBarPBar = new TProfile("profile_XMeanPt_YDelta_LambdaBarPBar", "#delta v.s. #bar{p_t} for LambdaBar-PBar", 35, 0, 3.5, -100., 100.);
    TProfile* profile_XMeanPt_YDelta_LambdaBarPBar_Eff = new TProfile("profile_XMeanPt_YDelta_LambdaBarPBar_Eff", "#delta v.s. #bar{p_t} for LambdaBar-PBar with eff. correction", 35, 0, 3.5, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaP", profile_XMeanPt_YDelta_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaP_Eff", profile_XMeanPt_YDelta_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaPBar", profile_XMeanPt_YDelta_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaPBar_Eff", profile_XMeanPt_YDelta_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaBarP", profile_XMeanPt_YDelta_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaBarP_Eff", profile_XMeanPt_YDelta_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaBarPBar", profile_XMeanPt_YDelta_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanPt_YDelta_LambdaBarPBar_Eff", profile_XMeanPt_YDelta_LambdaBarPBar_Eff));

    // Delta_vs_DeltaPt correlator
    TProfile* profile_XDeltaPt_YDelta_LambdaP = new TProfile("profile_XDeltaPt_YDelta_LambdaP", "#delta v.s. #delta p_t for Lambda-P", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YDelta_LambdaP_Eff = new TProfile("profile_XDeltaPt_YDelta_LambdaP_Eff", "#delta v.s. #delta p_t for Lambda-P with eff. correction", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YDelta_LambdaPBar = new TProfile("profile_XDeltaPt_YDelta_LambdaPBar", "#delta v.s. #delta p_t for Lambda-PBar", 50, 0, 5.0, -100., 100.); 
    TProfile* profile_XDeltaPt_YDelta_LambdaPBar_Eff = new TProfile("profile_XDeltaPt_YDelta_LambdaPBar_Eff", "#delta v.s. #delta p_t for Lambda-PBar with eff. correction", 50, 0, 5.0, -100., 100.);   
    TProfile* profile_XDeltaPt_YDelta_LambdaBarP = new TProfile("profile_XDeltaPt_YDelta_LambdaBarP", "#delta v.s. #delta p_t for LambdaBar-P", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YDelta_LambdaBarP_Eff = new TProfile("profile_XDeltaPt_YDelta_LambdaBarP_Eff", "#delta v.s. #delta p_t for LambdaBar-P with eff. correction", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YDelta_LambdaBarPBar = new TProfile("profile_XDeltaPt_YDelta_LambdaBarPBar", "#delta v.s. #delta p_t for LambdaBar-PBar", 50, 0, 5.0, -100., 100.);
    TProfile* profile_XDeltaPt_YDelta_LambdaBarPBar_Eff = new TProfile("profile_XDeltaPt_YDelta_LambdaBarPBar_Eff", "#delta v.s. #delta p_t for Lambda-PBar with eff. correction", 50, 0, 5.0, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaP", profile_XDeltaPt_YDelta_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaP_Eff", profile_XDeltaPt_YDelta_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaPBar", profile_XDeltaPt_YDelta_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaPBar_Eff", profile_XDeltaPt_YDelta_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaBarP", profile_XDeltaPt_YDelta_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaBarP_Eff", profile_XDeltaPt_YDelta_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaBarPBar", profile_XDeltaPt_YDelta_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaPt_YDelta_LambdaBarPBar_Eff", profile_XDeltaPt_YDelta_LambdaBarPBar_Eff));

    // Delta_vs_MeanEta correlator
    TProfile* profile_XMeanEta_YDelta_LambdaP = new TProfile("profile_XMeanEta_YDelta_LambdaP", "#delta v.s. #bar{#eta} for Lambda-P", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaP_Eff = new TProfile("profile_XMeanEta_YDelta_LambdaP_Eff", "#delta v.s. #bar{#eta} for Lambda-P with eff. correction", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaPBar = new TProfile("profile_XMeanEta_YDelta_LambdaPBar", "#delta v.s. #bar{#eta} for Lambda-PBar", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaPBar_Eff = new TProfile("profile_XMeanEta_YDelta_LambdaPBar_Eff", "#delta v.s. #bar{#eta} for Lambda-PBar with eff. correction", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaBarP = new TProfile("profile_XMeanEta_YDelta_LambdaBarP", "#delta v.s. #bar{#eta} for LambdaBar-P", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaBarP_Eff = new TProfile("profile_XMeanEta_YDelta_LambdaBarP_Eff", "#delta v.s. #bar{#eta} for LambdaBar-P with eff. correction", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaBarPBar = new TProfile("profile_XMeanEta_YDelta_LambdaBarPBar", "#delta v.s. #bar{#eta} for LambdaBar-PBar", 20, -1.0, 1.0, -100., 100.);
    TProfile* profile_XMeanEta_YDelta_LambdaBarPBar_Eff = new TProfile("profile_XMeanEta_YDelta_LambdaBarPBar_Eff", "#delta v.s. #bar{#eta} for LambdaBar-PBar with eff. correction", 20, -1.0, 1.0, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaP", profile_XMeanEta_YDelta_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaP_Eff", profile_XMeanEta_YDelta_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaPBar", profile_XMeanEta_YDelta_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaPBar_Eff", profile_XMeanEta_YDelta_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaBarP", profile_XMeanEta_YDelta_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaBarP_Eff", profile_XMeanEta_YDelta_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaBarPBar", profile_XMeanEta_YDelta_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XMeanEta_YDelta_LambdaBarPBar_Eff", profile_XMeanEta_YDelta_LambdaBarPBar_Eff));

    // Delta_vs_DeltaEta correlator
    TProfile* profile_XDeltaEta_YDelta_LambdaP = new TProfile("profile_XDeltaEta_YDelta_LambdaP", "#delta v.s. #delta #eta  for Lambda-P", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YDelta_LambdaP_Eff = new TProfile("profile_XDeltaEta_YDelta_LambdaP_Eff", "#delta v.s. #delta #eta  for Lambda-P with eff. correction", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YDelta_LambdaPBar = new TProfile("profile_XDeltaEta_YDelta_LambdaPBar", "#delta v.s. #delta #eta  for Lambda-PBar", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YDelta_LambdaPBar_Eff = new TProfile("profile_XDeltaEta_YDelta_LambdaPBar_Eff", "#delta v.s. #delta #eta  for Lambda-PBar with eff. correction", 20, 0, 2.0, -100., 100.);   
    TProfile* profile_XDeltaEta_YDelta_LambdaBarP = new TProfile("profile_XDeltaEta_YDelta_LambdaBarP", "#delta v.s. #delta #eta  for LambdaBar-P", 20, 0, 2.0, -100., 100.);
    TProfile* profile_XDeltaEta_YDelta_LambdaBarP_Eff = new TProfile("profile_XDeltaEta_YDelta_LambdaBarP_Eff", "#delta v.s. #delta #eta  for LambdaBar-P with eff. correction", 20, 0, 2.0, -100., 100.);
    TProfile* profile_XDeltaEta_YDelta_LambdaBarPBar = new TProfile("profile_XDeltaEta_YDelta_LambdaBarPBar", "#delta v.s. #delta #eta  for LambdaBar-PBar", 20, 0, 2.0, -100., 100.);
    TProfile* profile_XDeltaEta_YDelta_LambdaBarPBar_Eff = new TProfile("profile_XDeltaEta_YDelta_LambdaBarPBar_Eff", "#delta v.s. #delta #eta  for LambdaBar-PBar with eff. correction", 20, 0, 2.0, -100., 100.);
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaP", profile_XDeltaEta_YDelta_LambdaP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaP_Eff", profile_XDeltaEta_YDelta_LambdaP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaPBar", profile_XDeltaEta_YDelta_LambdaPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaPBar_Eff", profile_XDeltaEta_YDelta_LambdaPBar_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaBarP", profile_XDeltaEta_YDelta_LambdaBarP));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaBarP_Eff", profile_XDeltaEta_YDelta_LambdaBarP_Eff));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaBarPBar", profile_XDeltaEta_YDelta_LambdaBarPBar));
    histMap.insert(std::pair<std::string, TH1*>("profile_XDeltaEta_YDelta_LambdaBarPBar_Eff", profile_XDeltaEta_YDelta_LambdaBarPBar_Eff));
}

void StProtonLaGammaMaker::reconstructSubEventPlaneWithPhiWeightHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo){
    TVector2 mQ, mQ1, mQ2; // lack mLambdaPhi;

    vector< pair<double, StPriTrkInfo> > trkVec;
    const vector<StPriTrkInfo> vecPriTrks = evtInfo.VecPriTrks(); 
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();
    int Day = evtInfo.Day();
    double PVZ = evtInfo.Vz();
    double EWeight = evtInfo.EWeight(); 

    for(vector<StPriTrkInfo>::const_iterator iter_pri = vecPriTrks.begin(); iter_pri != vecPriTrks.end(); iter_pri++){
	// Get weights
	Int_t weight_index[3] = {1,};
	StCorrelationMaker::getWeightsIndex(weight_index, (*iter_pri).eta, PVZ, (*iter_pri).charge);
        Int_t phiBin = (Int_t)(((*iter_pri).phi + PI) / 2. / PI * phiBins);
        //std::cout << "phiBin = " << phiBin << std::endl;
        Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];
        //trk.phiWeight = phiWeight;//TODO


	//TODO:
        //if(m_ProtonTrkCuts->PassAllCuts(*iter_pri)) continue;
        if(!m_PriTrkCuts->PassAllCuts(*iter_pri)) continue;

        // Check if or not already excluded lambda daughters
	/*
        Bool_t kUse = 0;
        for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){
            Int_t dau1id = (Int_t)(*iter_la).id_dau1;
            Int_t dau2id = (Int_t)(*iter_la).id_dau2;
            if(dau1id == (*iter_pri).id || dau2id == (*iter_pri).id){
		break;
                kUse = 1; 
	    }
	} 
        if(kUse == 1)
	    continue;
	*/

        trkVec.push_back(pair<double, StPriTrkInfo>(phiWeight, *iter_pri));
    }//TODO: double check this loop

    std::random_shuffle(trkVec.begin(), trkVec.end());
    Long_t size_trkvec = trkVec.size();

    Float_t mQx = 0., mQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0.; // Full, east, west event planes
    for(Long_t i = 0; i < size_trkvec; i++){
        if(i > size_trkvec / 2){
            mQx1 += trkVec[i].first * trkVec[i].second.pt * cos(2. * trkVec[i].second.phi);
            mQy1 += trkVec[i].first * trkVec[i].second.pt * sin(2. * trkVec[i].second.phi);
	}
        else{
            mQx2 += trkVec[i].first * trkVec[i].second.pt * cos(2. * trkVec[i].second.phi);
            mQy2 += trkVec[i].first * trkVec[i].second.pt * sin(2. * trkVec[i].second.phi);
	}
    }

    //std::cout << "DEBUG: " << mQx1 << " " << mQy1 << " " << mQx2 << " " << mQy2 << std::endl;
    if(mQx1 == 0 || mQy1 == 0 || mQx2 == 0 || mQy2 == 0) return;

    mQ1.Set(mQx1, mQy1);
    mQ2.Set(mQx2, mQy2);

    Float_t TPC_RawTPCEPPhi_east = .5 * mQ1.Phi(); // Not really east, just one half of total tracks
    Float_t TPC_RawTPCEPPhi_west = .5 * mQ2.Phi();

    //std::cout << "DEBUG: Day = " << Day << endl;
    for(Int_t i = 0; i < order; i++){
	((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_EastEP"])->Fill(2 * i + 1, Day, cos(2 * (i + 1) * TPC_RawTPCEPPhi_east));
	((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_WestEP"])->Fill(2 * i + 1, Day, cos(2 * (i + 1) * TPC_RawTPCEPPhi_west));
	((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_EastEP"])->Fill(2 * i + 2, Day, sin(2 * (i + 1) * TPC_RawTPCEPPhi_east));
	((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_WestEP"])->Fill(2 * i + 2, Day, sin(2 * (i + 1) * TPC_RawTPCEPPhi_west));
    }
}

// Reconstruct shifted sub event plane
void StProtonLaGammaMaker::reconstructShiftedSubEventPlaneHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo){
    TVector2 mQ, mQ1, mQ2; // lack mLambdaPhi;
    Float_t mQx = 0., mQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0.; // Full, east, west event planes

    vector< pair<double, StPriTrkInfo> > trkVec;
    const vector<StPriTrkInfo> vecPriTrks = evtInfo.VecPriTrks();
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();
    double EWeight = evtInfo.EWeight();
    double PVZ = evtInfo.Vz();
    int Day = evtInfo.Day();
    int Day2 = evtInfo.Day2();
    int Day3 = evtInfo.Day3();
    for(vector<StPriTrkInfo>::const_iterator iter_pri = vecPriTrks.begin(); iter_pri != vecPriTrks.end(); iter_pri++){
        // Get weights
	Int_t weight_index[3] = {1,};
	StCorrelationMaker::getWeightsIndex(weight_index, (*iter_pri).eta, PVZ, (*iter_pri).charge);
        Int_t phiBin = (Int_t)(((*iter_pri).phi + PI) / 2. / PI * phiBins);
        Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];

        //TODO:
	if(!m_PriTrkCuts->PassAllCuts(*iter_pri)) continue;

        // Check if or not already excluded lambda daughters
	/*
	Bool_t kUse = 0;
        for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){
            Int_t dau1id = (Int_t)(*iter_la).id_dau1;
            Int_t dau2id = (Int_t)(*iter_la).id_dau2;
            if(dau1id == (*iter_pri).id || dau2id == (*iter_pri).id){
                kUse = 1; 
		break;
	    }
	}
        if(kUse == 1)
            continue;
        */
	trkVec.push_back(pair<double, StPriTrkInfo>(phiWeight, *iter_pri));
    }

    std::random_shuffle(trkVec.begin(), trkVec.end());
    Long_t size_trkvec = trkVec.size();

    for(Int_t i = 0; i < size_trkvec; i++){
        if(i > size_trkvec / 2){
            mQx1 += trkVec[i].first * trkVec[i].second.pt * cos(2. * trkVec[i].second.phi);
            mQy1 += trkVec[i].first * trkVec[i].second.pt * sin(2. * trkVec[i].second.phi);
	}
        else{
            mQx2 += trkVec[i].first * trkVec[i].second.pt * cos(2. * trkVec[i].second.phi);
            mQy2 += trkVec[i].first * trkVec[i].second.pt * sin(2. * trkVec[i].second.phi);
	}
    }

    if(mQx1 == 0 || mQy1 == 0 || mQx2 == 0 || mQy2 == 0) return;

    // Shifting correction up to 4th order, odd-cos, even-sin
    Double_t shift_correction_east[correction_terms], shift_correction_west[correction_terms], shift_correction_full[correction_terms]; 
    Int_t DayBin =  ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->GetYaxis()->FindBin(Day);
    for(Int_t i = 0; i < order; i++){
	shift_correction_full[2 * i] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->GetBinContent(2 * i + 1, DayBin);
	shift_correction_east[2 * i] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_EastEP"])->GetBinContent(2 * i + 1, DayBin);
	shift_correction_west[2 * i] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_WestEP"])->GetBinContent(2 * i + 1, DayBin);
	shift_correction_full[2 * i + 1] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->GetBinContent(2 * i + 2, DayBin);
	shift_correction_east[2 * i + 1] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_EastEP"])->GetBinContent(2 * i + 2, DayBin);
	shift_correction_west[2 * i + 1] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_WestEP"])->GetBinContent(2 * i + 2, DayBin);
    } 
      
    mQ1.Set(mQx1, mQy1);
    mQ2.Set(mQx2, mQy2);

    Float_t UnshiftedTPCEPPhi_full = .5 * mQ.Phi();
    Float_t UnshiftedTPCEPPhi_east = .5 * mQ1.Phi(); // Not really east, just one half of total tracks
    Float_t UnshiftedTPCEPPhi_west = .5 * mQ2.Phi();

    Float_t ShiftedTPCEPPhi_east = UnshiftedTPCEPPhi_east;
    Float_t ShiftedTPCEPPhi_west = UnshiftedTPCEPPhi_west;

    for(Int_t i = 0; i < order; i++){
        ShiftedTPCEPPhi_east += 2 * (-shift_correction_east[2 * i + 1] * cos(2 * (i + 1) * UnshiftedTPCEPPhi_east) + shift_correction_east[2 * i] * sin(2 * (i + 1) * UnshiftedTPCEPPhi_east)) / (Float_t)(2 * i + 2); 
        ShiftedTPCEPPhi_west += 2 * (-shift_correction_west[2 * i + 1] * cos(2 * (i + 1) * UnshiftedTPCEPPhi_west) + shift_correction_west[2 * i] * sin(2 * (i + 1) * UnshiftedTPCEPPhi_west)) / (Float_t)(2 * i + 2); 
    }

    if(ShiftedTPCEPPhi_east > PI) ShiftedTPCEPPhi_east -= PI;
    if(ShiftedTPCEPPhi_east < 0) ShiftedTPCEPPhi_east += PI;
    if(ShiftedTPCEPPhi_west > PI) ShiftedTPCEPPhi_west -= PI;
    if(ShiftedTPCEPPhi_west < 0) ShiftedTPCEPPhi_west += PI;

    // Combine to form full ep 
    mQx = mQ1.Mod() * cos(2 * ShiftedTPCEPPhi_east) + mQ2.Mod() * cos(2 * ShiftedTPCEPPhi_west);
    mQy = mQ1.Mod() * sin(2 * ShiftedTPCEPPhi_east) + mQ2.Mod() * sin(2 * ShiftedTPCEPPhi_west);
    mQ.Set(mQx, mQy);
    Float_t SubEPShiftedTPCEPPhi_full = .5 * mQ.Phi();
    for(Int_t i = 0; i < order; i++){
	((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->Fill(2 * i + 1, Day, cos(2 * (i + 1) * SubEPShiftedTPCEPPhi_full));
	((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->Fill(2 * i + 2, Day, sin(2 * (i + 1) * SubEPShiftedTPCEPPhi_full));
    }

    // Event resolution histogram
    ((TProfile*)histMap["profile_eventplane_resolution"])->Fill(1, 100 * cos(2 * (ShiftedTPCEPPhi_east - ShiftedTPCEPPhi_west)), EWeight);
    ((TProfile*)histMap["profile_eventplane_resolution_noEW"])->Fill(1, 100 * cos(2 * (ShiftedTPCEPPhi_east - ShiftedTPCEPPhi_west)));
}

// Reconstruct full event plane
void StProtonLaGammaMaker::reconstructShiftedFullEventPlaneHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo, TVector2& ShiftedFullTPCEPPhi){
    TVector2 mQ, mQ1, mQ2; // lack mLambdaPhi;
    Float_t mQx = 0., mQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0.; // full, east, west event planes

    std::vector< pair<double, StPriTrkInfo> > trkVec;
    const vector<StPriTrkInfo> vecPriTrks = evtInfo.VecPriTrks();
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();
    double PVZ = evtInfo.Vz();
    double EWeight = evtInfo.EWeight();
    int Day = evtInfo.Day();
    int Day2 = evtInfo.Day2();
    int Day3 = evtInfo.Day3();
    for(vector<StPriTrkInfo>::const_iterator iter_pri = vecPriTrks.begin(); iter_pri != vecPriTrks.end(); iter_pri++){
        // Get weights
	Int_t weight_index[3] = {1,};
	StCorrelationMaker::getWeightsIndex(weight_index, (*iter_pri).eta, PVZ, (*iter_pri).charge);
        Int_t phiBin = (Int_t)(((*iter_pri).phi + PI) / 2. / PI * phiBins);
        Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];
        
	//TODO:
	if(!m_PriTrkCuts->PassAllCuts(*iter_pri)) continue;

	// Check if or not already excluded lambda daughters
	/*
	Bool_t kUse = 0;
        for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){
            Int_t dau1id = (Int_t)(*iter_la).id_dau1;
            Int_t dau2id = (Int_t)(*iter_la).id_dau2;
            if(dau1id == (*iter_pri).id || dau2id == (*iter_pri).id){
                kUse = 1; 
		break;
	    }
	} 
        if(kUse == 1)
            continue;
        */

	trkVec.push_back(pair<double, StPriTrkInfo>(phiWeight, *iter_pri));
        //std::cout << "<------- End reconstructing shifted full event plane helper loop! ---------->" << std::endl;
    }

    std::random_shuffle(trkVec.begin(), trkVec.end());
    Long_t size_trkvec = trkVec.size();

    //std::cout << size_trkvec << " is the size of vector!" << std::endl;
    for(Int_t i = 0; i < size_trkvec; i++){
        if(i > size_trkvec / 2){
            mQx1 += trkVec[i].first * trkVec[i].second.pt * cos(2. * trkVec[i].second.phi);
            mQy1 += trkVec[i].first * trkVec[i].second.pt * sin(2. * trkVec[i].second.phi);
	}
        else{
            mQx2 += trkVec[i].first * trkVec[i].second.pt * cos(2. * trkVec[i].second.phi);
            mQy2 += trkVec[i].first * trkVec[i].second.pt * sin(2. * trkVec[i].second.phi);
	}
    }

    if(mQx1 == 0 || mQy1 == 0 || mQx2 == 0 || mQy2 == 0) return;

    // Shifting correction up to 4th order, odd-cos, even-sin
    Float_t shift_correction_east[correction_terms], shift_correction_west[correction_terms], shift_correction_full[correction_terms]; 
    Int_t DayBin = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->GetYaxis()->FindBin(Day);

    //cout << DayBin << " DayBin" << endl;
    for(Int_t i = 0; i < order; i++){
	shift_correction_full[2 * i] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->GetBinContent(2 * i + 1, DayBin);
	shift_correction_east[2 * i] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_EastEP"])->GetBinContent(2 * i + 1, DayBin);
	shift_correction_west[2 * i] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_WestEP"])->GetBinContent(2 * i + 1, DayBin);
	shift_correction_full[2 * i + 1] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_FullEP"])->GetBinContent(2 * i + 2, DayBin);
	shift_correction_east[2 * i + 1] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_EastEP"])->GetBinContent(2 * i + 2, DayBin);
	shift_correction_west[2 * i + 1] = ((TProfile2D*)histMap["prof2_XOrder_YDay_ZCorrectionTerm_WestEP"])->GetBinContent(2 * i + 2, DayBin);
    } 
      
    mQ1.Set(mQx1, mQy1);
    mQ2.Set(mQx2, mQy2);

    //Float_t UnshiftedTPCEPPhi_full = .5 * mQ.Phi();
    Float_t UnshiftedTPCEPPhi_east = .5 * mQ1.Phi(); // Not really east, just one half of total tracks
    Float_t UnshiftedTPCEPPhi_west = .5 * mQ2.Phi();

    Float_t correction_east = 0;
    Float_t correction_west = 0;

    // Subevent before flattening
    ((TH1F*)histMap["h1f_before_Flattened_EastEPPhi"])->Fill(UnshiftedTPCEPPhi_east, EWeight); 
    ((TH1F*)histMap["h1f_before_Flattened_WestEPPhi"])->Fill(UnshiftedTPCEPPhi_west, EWeight); 
    
    for(Int_t i = 0; i < order; i++){
        //ShiftedTPCEPPhi_east += 2 * (-shift_correction_east[2 * i + 1] * cos(2 * UnshiftedTPCEPPhi_east) + shift_correction_east[2 * i] * sin(2 * UnshiftedTPCEPPhi_east)) / (Float_t)(2 * i + 2); 
        correction_east += 2 * (-shift_correction_east[2 * i + 1] * cos(2 * (i + 1) * UnshiftedTPCEPPhi_east) + shift_correction_east[2 * i] * sin(2 * (i + 1) * UnshiftedTPCEPPhi_east)) / (Float_t)(2 * i + 2); 
	//std::cout << "DEBUG: shift_correction_east = " << shift_correction_east[2 * i + 1] << ", correction_east = " << shift_correction_east[2 * i + 1] << std::endl;
        //ShiftedTPCEPPhi_west += 2 * (-shift_correction_west[2 * i + 1] * cos(2 * UnshiftedTPCEPPhi_west) + shift_correction_west[2 * i] * sin(2 * UnshiftedTPCEPPhi_west)) / (Float_t)(2 * i + 2); 
        correction_west += 2 * (-shift_correction_west[2 * i + 1] * cos(2 * (i + 1) * UnshiftedTPCEPPhi_west) + shift_correction_west[2 * i] * sin(2 * (i + 1) * UnshiftedTPCEPPhi_west)) / (Float_t)(2 * i + 2); 
	//std::cout << "DEBUG: shift_correction_west = " << shift_correction_west[2 * i + 1] << ", correction_west = " << shift_correction_west[2 * i + 1] << std::endl;
    }

    //std::cout << "DEBUG: correction_east = " << correction_east << ", correction_west = " << correction_west << std::endl;
    Float_t ShiftedTPCEPPhi_east = UnshiftedTPCEPPhi_east + correction_east;
    Float_t ShiftedTPCEPPhi_west = UnshiftedTPCEPPhi_west + correction_west;

    if(ShiftedTPCEPPhi_east > PI) ShiftedTPCEPPhi_east -= PI;
    if(ShiftedTPCEPPhi_east < 0) ShiftedTPCEPPhi_east += PI;
    if(ShiftedTPCEPPhi_west > PI) ShiftedTPCEPPhi_west -= PI;
    if(ShiftedTPCEPPhi_west < 0) ShiftedTPCEPPhi_west += PI;

    // Subevent Flat check
    ((TH1F*)histMap["h1f_Flattened_EastEPPhi"])->Fill(ShiftedTPCEPPhi_east, EWeight);
    ((TH1F*)histMap["h1f_Flattened_WestEPPhi"])->Fill(ShiftedTPCEPPhi_west, EWeight);

    // Combine to form full ep 
    mQx = mQ1.Mod() * cos(2 * ShiftedTPCEPPhi_east) + mQ2.Mod() * cos(2 * ShiftedTPCEPPhi_west);
    mQy = mQ1.Mod() * sin(2 * ShiftedTPCEPPhi_east) + mQ2.Mod() * sin(2 * ShiftedTPCEPPhi_west);
    mQ.Set(mQx, mQy);
    Float_t UnshiftedTPCEPPhi_Full = .5 * mQ.Phi(); 

    Float_t correction_full = 0;

    ((TH1F*)histMap["h1f_before_Flattened_FullEPPhi"])->Fill(UnshiftedTPCEPPhi_Full, EWeight);
    for(Int_t i = 0; i < order; i++)
        correction_full += 2 * (-shift_correction_full[2 * i + 1] * cos(2 * (i + 1) * UnshiftedTPCEPPhi_Full) + shift_correction_full[2 * i] * sin(2 * (i + 1) * UnshiftedTPCEPPhi_Full)) / (Float_t)(2 * i + 2); 

    double newFullTPCEPPhi = UnshiftedTPCEPPhi_Full + correction_full;

    if(newFullTPCEPPhi > PI) newFullTPCEPPhi -= PI;
    if(newFullTPCEPPhi < 0) newFullTPCEPPhi += PI;

    // Full event flat check 
    ((TH1F*)histMap["h1f_Flattened_FullEPPhi"])->Fill(newFullTPCEPPhi, EWeight);
    ShiftedFullTPCEPPhi.Set(mQ.Mod() * cos(2 * newFullTPCEPPhi), mQ.Mod() * sin(2 * newFullTPCEPPhi));
}

// Compute Correlators
void StProtonLaGammaMaker::computeCorrelatorsHelper(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo, TVector2 ShiftedFullTPCEPPhi){
    //std::cout << "Start Computing Correlators!!" << std::endl;
    TVector2 mQ, mQ1, mQ2; // lack mLambdaPhi;
    Float_t mQx = 0., mQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0.; // full, east, west event planes

// Fill out some qa histograms
    const vector<StPriTrkInfo> vecPriTrks = evtInfo.VecPriTrks();
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();
    double PVZ = evtInfo.Vz();
    double EWeight = evtInfo.EWeight();
    int Day = evtInfo.Day();
    int Day2 = evtInfo.Day2();
    int Day3 = evtInfo.Day3();
    for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){

        if(!m_LambdaTrkCuts->PassAllCuts(*iter_la)) continue;
	std::string name = getCorrectionHistName("Beta", (*iter_la).eta, PVZ, 1., 'c');// Note: pay attention to fillBetaHists()
	Float_t shift_Lambda[correction_terms];
	Float_t shiftedLambdaPhi = (*iter_la).phi;

	((TH1F*)histMap["h1f_before_Corrections_BetaPhi"])->Fill((*iter_la).phi, EWeight);
    
	Int_t Day2Bin =  ((TProfile2D*)histMap[name])->GetYaxis()->FindBin(Day2);
	for(Int_t kk = 0; kk < order; kk++){
	    shift_Lambda[2 * kk] = ((TProfile2D*)histMap.find(name)->second)->GetBinContent(2 * kk + 1, Day2Bin);
	    shift_Lambda[2 * kk + 1] = ((TProfile2D*)histMap.find(name)->second)->GetBinContent(2 * kk + 2, Day2Bin);
	    shiftedLambdaPhi += (-2 * shift_Lambda[2 * kk + 1] * cos((kk + 1) * (*iter_la).phi) / (kk + 1) + 2 * shift_Lambda[2 * kk] * sin((kk + 1) * (*iter_la).phi) / (kk + 1));
	}

	if(shiftedLambdaPhi > PI) shiftedLambdaPhi -= 2 * PI;
	if(shiftedLambdaPhi < -PI) shiftedLambdaPhi += 2 * PI;

        // Check if the daugthers are used to reconstruct the full event plane;
	int flag_ep_dau1 = 0, flag_ep_dau2 = 0;
	if((*iter_la).prmatch_dau1 && m_PriTrkCuts->PassAllCuts((*iter_la), 1))
	    flag_ep_dau1 = 1;
	if((*iter_la).prmatch_dau2 && m_PriTrkCuts->PassAllCuts((*iter_la), 2))
	    flag_ep_dau2 = 1;

	TVector2 fullep = ShiftedFullTPCEPPhi;
        if(flag_ep_dau1){
	    Int_t weight_index[3] = {1,};
	    StCorrelationMaker::getWeightsIndex(weight_index, (*iter_la).eta_dau1, PVZ, (*iter_la).charge_dau1);
	    Int_t phiBin = (Int_t)(((*iter_la).phi_dau1 + PI) / 2. / PI * phiBins);
	    Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];
        
	    double mqx = fullep.X();
	    double mqy = fullep.Y();
	    mqx -= phiWeight * (*iter_la).pt_dau1 * cos(2. * (*iter_la).phi_dau1);
	    mqy -= phiWeight * (*iter_la).pt_dau1 * sin(2. * (*iter_la).phi_dau1);
	    fullep.Set(mqx, mqy);
	}

        if(flag_ep_dau2){
	    Int_t weight_index[3] = {1,};
	    StCorrelationMaker::getWeightsIndex(weight_index, (*iter_la).eta_dau2, PVZ, (*iter_la).charge_dau2);
	    Int_t phiBin = (Int_t)(((*iter_la).phi_dau2 + PI) / 2. / PI * phiBins);
	    Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];
        
	    double mqx = fullep.X();
	    double mqy = fullep.Y();
	    mqx -= phiWeight * (*iter_la).pt_dau2 * cos(2. * (*iter_la).phi_dau2);
	    mqy -= phiWeight * (*iter_la).pt_dau2 * sin(2. * (*iter_la).phi_dau2);
	    fullep.Set(mqx, mqy);
	}

        ((TH1D*)histMap[getCorrectionHistName("Beta", (*iter_la).eta, PVZ, 1., 'a')])->Fill(shiftedLambdaPhi, EWeight);

	// Fill v2
	((TProfile*)histMap["profile_XPt_Yv2_LambdaUnWeighted"])->Fill((*iter_la).pt, 100 * cos(2 * (shiftedLambdaPhi - fullep.Phi() * .5))); 
	((TProfile*)histMap["profile_XPt_Yv2_LambdaWeighted"])->Fill((*iter_la).pt, 100 * cos(2 * (shiftedLambdaPhi - fullep.Phi() * .5)), EWeight);     
    }

    for(vector<StPriTrkInfo>::const_iterator iter_pri = vecPriTrks.begin(); iter_pri != vecPriTrks.end(); iter_pri++){
        StPriTrkInfo pritrk = *iter_pri;
        Bool_t kUse = 0;
        for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){
            if(pritrk.id == (*iter_la).id_dau1 || pritrk.id == (*iter_la).id_dau2){
                kUse = 1;
                break;
	    }
	}
	if(kUse == 1) continue;

        if(m_ProtonTrkCuts->PassAllCuts(*iter_pri))
	{  
	    Float_t En_proton = sqrt(ProtonPDGMass * ProtonPDGMass + pow(pritrk.pt * cosh(pritrk.eta), 2));
            Float_t Theta_proton = 2. * atan(exp(-pritrk.eta));
            Float_t Phi_Proton = pritrk.phi;

            if(Phi_Proton > PI) Phi_Proton -= 2 * PI;
            if(Phi_Proton < -PI) Phi_Proton += 2 * PI;

	    std::string name = getCorrectionHistName("Alpha", pritrk.eta, PVZ, pritrk.charge, 'c');
	    Float_t shift_proton[correction_terms];
	    Float_t shiftedProtonPhi = pritrk.phi;
            //((TH1D*)histMap[getCorrectionHistName("Alpha", Eta, PVZ, Charge, 'b')])->Fill(Phi_Proton, EWeight);

	    ((TH1F*)histMap["h1f_before_Corrections_AlphaPhi"])->Fill((*iter_pri).phi, EWeight);

	    Int_t Day2Bin =  ((TProfile2D*)histMap[name])->GetYaxis()->FindBin(Day2);
	    for(Int_t kk = 0; kk < order; kk++){
                //std::cout << kk << "th order in proton" << std::endl;
		shift_proton[2 * kk] = ((TProfile2D*)histMap.find(name)->second)->GetBinContent(2 * kk + 1, Day2Bin);
		shift_proton[2 * kk + 1] = ((TProfile2D*)histMap.find(name)->second)->GetBinContent(2 * kk + 2, Day2Bin);
		shiftedProtonPhi += -2 * shift_proton[2 * kk + 1] * cos((kk + 1) * (*iter_pri).phi) / (kk + 1) + 2 * shift_proton[2 * kk] * sin((kk + 1) * (*iter_pri).phi) / (kk + 1);
	    }

	    if(shiftedProtonPhi > PI) shiftedProtonPhi -= 2 * PI;
	    if(shiftedProtonPhi < -PI) shiftedProtonPhi += 2 * PI;

            // Proton flat check 
	    ((TH1F*)histMap["h1f_after_Corrections_AlphaPhi"])->Fill(shiftedProtonPhi, EWeight);

            // Full event plane correction
            TVector2 fullep = ShiftedFullTPCEPPhi; 
	    if(m_PriTrkCuts->PassAllCuts(*iter_pri)){
		Int_t weight_index[3] = {1,};
		StCorrelationMaker::getWeightsIndex(weight_index, (*iter_pri).eta, PVZ, (*iter_pri).charge);
		Int_t phiBin = (Int_t)(((*iter_pri).phi + PI) / 2. / PI * phiBins);
		Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];

		double mqx = fullep.X();
		double mqy = fullep.Y();
		mqx -= phiWeight * (*iter_pri).pt * cos(2. * (*iter_pri).phi);
		mqy -= phiWeight * (*iter_pri).pt * sin(2. * (*iter_pri).phi);
		fullep.Set(mqx, mqy);
	    }

            // Fill v2
            ((TProfile*)histMap["profile_XPt_Yv2_ProtonUnWeighted"])->Fill(pritrk.pt, 100 * cos(2 * (shiftedProtonPhi - fullep.Phi() * .5)));     
            ((TProfile*)histMap["profile_XPt_Yv2_ProtonWeighted"])->Fill(pritrk.pt, 100 * cos(2 * (shiftedProtonPhi - fullep.Phi() * .5)), EWeight);     
            Float_t Charge_Proton = pritrk.charge;

	    //std::cout << "before NLambda loop" << std::endl;
	    for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){
                if(!m_LambdaTrkCuts->PassAllCuts(*iter_la)) continue; 
		std::string name = getCorrectionHistName("Beta", (*iter_la).eta, PVZ, (*iter_la).charge, 'c');
		Float_t shift_Lambda[correction_terms];
		Float_t shiftedLambdaPhi = (*iter_la).phi;

                ((TH1F*)histMap["h1f_before_Corrections_BetaPhi"])->Fill((*iter_la).phi, EWeight);

		Day2Bin =  ((TProfile2D*)histMap[name])->GetYaxis()->FindBin(Day2);
		for(Int_t kk = 0; kk < order; kk++){
		    shift_Lambda[2 * kk] = ((TProfile2D*)histMap.find(name)->second)->GetBinContent(2 * kk + 1, Day2Bin);
		    shift_Lambda[2 * kk + 1] = ((TProfile2D*)histMap.find(name)->second)->GetBinContent(2 * kk + 2, Day2Bin);
		    shiftedLambdaPhi += (-2 * shift_Lambda[2 * kk + 1] * cos((kk + 1) * (*iter_la).phi) / (kk + 1) + 2 * shift_Lambda[2 * kk] * sin((kk + 1) * (*iter_la).phi) / (kk + 1));
		}

		if(shiftedLambdaPhi > PI) shiftedLambdaPhi -= 2 * PI; 
		if(shiftedLambdaPhi < -PI) shiftedLambdaPhi += 2 * PI; 
              
                // Check if it is already flattened 
                ((TH1F*)histMap["h1f_after_Corrections_BetaPhi"])->Fill(shiftedLambdaPhi, EWeight);

                // Correct the full event-plane orientation: eliminate auto-correlation
		int flag_ep_dau1 = 0, flag_ep_dau2 = 0;
		if((*iter_la).prmatch_dau1 && m_PriTrkCuts->PassAllCuts((*iter_la), 1))
		    flag_ep_dau1 = 1;
		if((*iter_la).prmatch_dau2 && m_PriTrkCuts->PassAllCuts((*iter_la), 2))
		    flag_ep_dau2 = 1;

		TVector2 final_fullep = fullep;
		if(flag_ep_dau1){
		    Int_t weight_index[3] = {1,};
		    StCorrelationMaker::getWeightsIndex(weight_index, (*iter_la).eta_dau1, PVZ, (*iter_la).charge_dau1);
		    Int_t phiBin = (Int_t)(((*iter_la).phi_dau1 + PI) / 2. / PI * phiBins);
		    Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];

		    double mqx = final_fullep.X();
		    double mqy = final_fullep.Y();
		    mqx -= phiWeight * (*iter_la).pt_dau1 * cos(2. * (*iter_la).phi_dau1);
		    mqy -= phiWeight * (*iter_la).pt_dau1 * sin(2. * (*iter_la).phi_dau1);
		    final_fullep.Set(mqx, mqy);
		}

		if(flag_ep_dau2){
		    Int_t weight_index[3] = {1,};
		    StCorrelationMaker::getWeightsIndex(weight_index, (*iter_la).eta_dau2, PVZ, (*iter_la).charge_dau2);
		    Int_t phiBin = (Int_t)(((*iter_la).phi_dau2 + PI) / 2. / PI * phiBins);
		    Float_t phiWeight = m_PrimaryTracksPhiWeight[weight_index[0]][weight_index[1]][weight_index[2]][phiBin];

		    double mqx = final_fullep.X();
		    double mqy = final_fullep.Y();
		    mqx -= phiWeight * (*iter_la).pt_dau2 * cos(2. * (*iter_la).phi_dau2);
		    mqy -= phiWeight * (*iter_la).pt_dau2 * sin(2. * (*iter_la).phi_dau2);
		    final_fullep.Set(mqx, mqy);
		}


		// Compute correlators and fill histograms
		Float_t gamma = cos(shiftedProtonPhi + shiftedLambdaPhi - 2 * final_fullep.Phi());
                Float_t delta = cos(shiftedProtonPhi - shiftedLambdaPhi);
                
		//std::cout << "DEBUG: " << "Gamma = " << gamma << std::endl;
		if(!isnan(gamma)){
		    double Pt = (*iter_pri).pt;
		    double Eta = (*iter_pri).eta;
		    double Pt_Lambda = (*iter_la).pt;
		    double Eta_Lambda = (*iter_la).eta; 
		    //cout << (*iter_pri).charge << " ... " << (*iter_la).charge << endl;
		    if((*iter_pri).charge * (*iter_la).charge > 0){ // Same-sign
			Double_t eff = (GetAlphaEffMaker()->GetEfficiency(this->GetCentrality(), Pt)) * (GetBetaEffMaker()->GetEfficiency(this->GetCentrality(), Pt_Lambda));

			((TProfile*)histMap["profile_XDay2_YSameSignGamma"])->Fill(Day2, 100 * gamma);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma"])->Fill(1, 100 * gamma, EWeight);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_noEW"])->Fill(1, 100 * gamma);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff"])->Fill(1, 100 * gamma, EWeight / eff);
                        ((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta"])->Fill(1, 100 * delta, EWeight); 
                        ((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta_Eff"])->Fill(1, 100 * delta, EWeight / eff);

                        ((TProfile*)histMap["profile_XMeanPt_YGamma_LambdaP"])->Fill(.5 * (Pt_Lambda + Pt), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XMeanPt_YGamma_LambdaP_Eff"])->Fill(.5 * (Pt_Lambda + Pt), 100 * gamma, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaPt_YGamma_LambdaP"])->Fill(fabs(Pt_Lambda - Pt), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XDeltaPt_YGamma_LambdaP_Eff"])->Fill(fabs(Pt_Lambda - Pt), 100 * gamma, EWeight / eff);

			((TProfile*)histMap["profile_XMeanEta_YGamma_LambdaP"])->Fill(.5 * (Eta_Lambda + Eta), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XMeanEta_YGamma_LambdaP_Eff"])->Fill(.5 * (Eta_Lambda + Eta), 100 * gamma, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaEta_YGamma_LambdaP"])->Fill(fabs(Eta_Lambda - Eta), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XDeltaEta_YGamma_LambdaP_Eff"])->Fill(fabs(Eta_Lambda - Eta), 100 * gamma, EWeight / eff);

			((TProfile*)histMap["profile_XMeanPt_YDelta_LambdaP"])->Fill(.5 * (Pt_Lambda + Pt), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XMeanPt_YDelta_LambdaP_Eff"])->Fill(.5 * (Pt_Lambda + Pt), 100 * delta, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaPt_YDelta_LambdaP"])->Fill(fabs(Pt_Lambda - Pt), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XDeltaPt_YDelta_LambdaP_Eff"])->Fill(fabs(Pt_Lambda - Pt), 100 * delta, EWeight / eff);

			((TProfile*)histMap["profile_XMeanEta_YDelta_LambdaP"])->Fill(.5 * (Eta_Lambda + Eta), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XMeanEta_YDelta_LambdaP_Eff"])->Fill(.5 * (Eta_Lambda + Eta), 100 * delta, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaEta_YDelta_LambdaP"])->Fill(fabs(Eta_Lambda - Eta), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XDeltaEta_YDelta_LambdaP_Eff"])->Fill(fabs(Eta_Lambda - Eta), 100 * delta, EWeight / eff);
		    }
		    if((*iter_pri).charge * (*iter_la).charge < 0){
			Double_t eff = (GetAlphaBarEffMaker()->GetEfficiency(this->GetCentrality(), Pt)) * (GetBetaEffMaker()->GetEfficiency(this->GetCentrality(), Pt_Lambda));
			((TProfile*)histMap["profile_XDay2_YOppoSignGamma"])->Fill(Day2, 100 * gamma);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma"])->Fill(2, 100 * gamma, EWeight);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_noEW"])->Fill(2, 100 * gamma);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YGamma_Eff"])->Fill(2, 100 * gamma, EWeight / eff);
			((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta"])->Fill(2, 100 * delta, EWeight); 
                        ((TProfile*)histMap["profile_1LambdaP_2LambdaPBar_3LambdaBarP_4LambdaBarPBar_YDelta_Eff"])->Fill(2, 100 * delta, EWeight / eff);

			((TProfile*)histMap["profile_XMeanPt_YGamma_LambdaPBar"])->Fill(.5 * (Pt_Lambda + Pt), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XMeanPt_YGamma_LambdaPBar_Eff"])->Fill(.5 * (Pt_Lambda + Pt), 100 * gamma, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaPt_YGamma_LambdaPBar"])->Fill(fabs(Pt_Lambda - Pt), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XDeltaPt_YGamma_LambdaPBar_Eff"])->Fill(fabs(Pt_Lambda - Pt), 100 * gamma, EWeight / eff);

			((TProfile*)histMap["profile_XMeanEta_YGamma_LambdaPBar"])->Fill(.5 * (Eta_Lambda + Eta), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XMeanEta_YGamma_LambdaPBar_Eff"])->Fill(.5 * (Eta_Lambda + Eta), 100 * gamma, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaEta_YGamma_LambdaPBar"])->Fill(fabs(Eta_Lambda - Eta), 100 * gamma, EWeight);
                        ((TProfile*)histMap["profile_XDeltaEta_YGamma_LambdaPBar_Eff"])->Fill(fabs(Eta_Lambda - Eta), 100 * gamma, EWeight / eff);

			((TProfile*)histMap["profile_XMeanPt_YDelta_LambdaPBar"])->Fill(.5 * (Pt_Lambda + Pt), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XMeanPt_YDelta_LambdaPBar_Eff"])->Fill(.5 * (Pt_Lambda + Pt), 100 * delta, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaPt_YDelta_LambdaPBar"])->Fill(fabs(Pt_Lambda - Pt), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XDeltaPt_YDelta_LambdaPBar_Eff"])->Fill(fabs(Pt_Lambda - Pt), 100 * delta, EWeight / eff);

			((TProfile*)histMap["profile_XMeanEta_YDelta_LambdaPBar"])->Fill(.5 * (Eta_Lambda + Eta), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XMeanEta_YDelta_LambdaPBar_Eff"])->Fill(.5 * (Eta_Lambda + Eta), 100 * delta, EWeight / eff);
                        ((TProfile*)histMap["profile_XDeltaEta_YDelta_LambdaPBar"])->Fill(fabs(Eta_Lambda - Eta), 100 * delta, EWeight);
                        ((TProfile*)histMap["profile_XDeltaEta_YDelta_LambdaPBar_Eff"])->Fill(fabs(Eta_Lambda - Eta), 100 * delta, EWeight / eff);

		    }
		}
	    }
	}
    }
}

void StProtonLaGammaMaker::fillPrimaryTracksPhiHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo){ // Exclude correlation protons or Lambda daughters
    Int_t NTracks = (Int_t)evtInfo.NPTracks();
    const vector<StPriTrkInfo> vecPriTrks = evtInfo.VecPriTrks(); 
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();

    double EWeight = evtInfo.EWeight();
    double PVZ = evtInfo.Vz();
    double Day = evtInfo.Day();
    double Day2 = evtInfo.Day2();
    double Day3 = evtInfo.Day3();
    for(vector<StPriTrkInfo>::const_iterator iter_pri = vecPriTrks.begin(); 
	    iter_pri != vecPriTrks.end(); iter_pri++){

        StPriTrkInfo pritrk = *iter_pri;
	Bool_t kUse = 0; // Eliminate Lambda daughters
        for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); 
		iter_la != vecLambdaTrks.end(); iter_la++){
            if((*iter_pri).id == (*iter_la).id_dau1 || (*iter_pri).id == (*iter_la).id_dau2){
                kUse = 1; 
                break;
	    }
	}
	if(kUse == 1) continue;

        //cout << "happy so far in fill primary" << endl;
        if(m_ProtonTrkCuts->PassAllCuts(*iter_pri)) continue;
	if(!m_PriTrkCuts->PassAllCuts(*iter_pri)) continue;
	std::string name = getCorrectionHistName("PrimaryTrk", pritrk.eta, PVZ, pritrk.charge, 'b');
	    //std::string name = getCorrectionHistName("Alpha", Eta, PVZ, Charge, 'b');
	((TH1D*)(histMap.find(name)->second))->Fill(pritrk.phi, EWeight);
    }
}

//void StProtonLaGammaMaker::fillAlphaHists(std::map<std::string, TH1*>& histMap, TChain* chain, Float_t EWeight, Float_t PVZ, Long_t Run) {
void StProtonLaGammaMaker::fillAlphaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo) {
    // Iteration over all of the tracks
    double EWeight = evtInfo.EWeight(); 
    double PVZ = evtInfo.Vz();
    const vector<StPriTrkInfo> vecPriTrks = evtInfo.VecPriTrks();
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();
    int Day = evtInfo.Day();
    int Day2 = evtInfo.Day2();
    int Day3 = evtInfo.Day3();
    for(vector<StPriTrkInfo>::const_iterator iter_pri = vecPriTrks.begin(); iter_pri != vecPriTrks.end(); iter_pri++){
        // Eliminate Lambda daughters
	Bool_t kUse = 0;
	for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin(); iter_la != vecLambdaTrks.end(); iter_la++){
	    if((*iter_pri).id == (*iter_la).id_dau1 || (*iter_pri).id == (*iter_la).id_dau2){
		kUse = 1; 
		break;
	    }
	}
	if(kUse == 1) continue;

	if(m_ProtonTrkCuts->PassAllCuts(*iter_pri)){
            float pt = (*iter_pri).pt;
	    float eta = (*iter_pri).eta;
            float phi = (*iter_pri).phi;
	    float charge = (*iter_pri).charge;
	    Float_t En_proton = sqrt(ProtonPDGMass * ProtonPDGMass + pow(pt * cosh(eta), 2));
	    Float_t Theta_proton = 2. * atan(exp(-eta));
	    ((TH2F*)(histMap.find("h2f_XEta_YPt_ProtonUnWeighted")->second))->Fill(eta, pt);
	    ((TH2F*)(histMap.find("h2f_XEta_YPt_ProtonWeighted")->second))->Fill(eta, pt, EWeight);
	    // Fill the before correction phi histograms

	    std::string name = getCorrectionHistName("Alpha", eta, PVZ, charge, 'b');
	    ((TH1D*)(histMap.find(name)->second))->Fill((*iter_pri).phi, EWeight);

	    // Fill the correction term histograms
	    name = getCorrectionHistName("Alpha", eta, PVZ, charge, 'c');
	    for(Int_t kk = 0; kk < order; kk++){
		((TProfile2D*)(histMap.find(name)->second))->Fill(2 * kk + 1, Day2, cos(kk * phi + phi), EWeight);
		((TProfile2D*)(histMap.find(name)->second))->Fill(2 * kk + 2, Day2, sin(kk * phi + phi), EWeight);
	    }
	}
	else
	    continue;
    }
}

void StProtonLaGammaMaker::fillBetaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo) {
    // Iteration over all of the tracks
    double EWeight = evtInfo.EWeight(); 
    double PVZ = evtInfo.Vz();
    int Day = evtInfo.Day();
    int Day2 = evtInfo.Day2();
    const vector<StV0TrkInfo> vecLambdaTrks = evtInfo.VecBetaTrks();
    //for(Int_t i = 0; i < NLambda; i++){
    m_LambdaTrkCuts->Dump();
    for(vector<StV0TrkInfo>::const_iterator iter_la = vecLambdaTrks.begin();
                             iter_la != vecLambdaTrks.end(); iter_la++){
        double eta = (*iter_la).eta;
	double mass = (*iter_la).mass;
	double pt = (*iter_la).pt;
	int charge = (*iter_la).charge;
	double phi = (*iter_la).phi;

        if(!m_LambdaTrkCuts->PassAllCuts(*iter_la)) continue;

	//TODO:hEtaPtDist->Fill(Eta,Pt,EWeight);
	((TH2F*)(histMap.find("h2f_XEta_YPt_LambdaUnWeighted")->second))->Fill(eta, pt);
	((TH2F*)(histMap.find("h2f_XEta_YPt_LambdaWeighted")->second))->Fill(eta, pt, EWeight);
	((TH1F*)(histMap.find("h1f_LambdaMassUnWeighted")->second))->Fill(mass);
	((TH1F*)(histMap.find("h1f_LambdaMassWeighted")->second))->Fill(mass, EWeight);
	//if(fabs(Mass_Lambda - LambdaPDGMass) > LambdaCuts.ParentMassWidth_UpperLimit) continue;

	// Fill the before correction phi histograms
        std::string name = getCorrectionHistName("Beta", eta, PVZ, charge, 'b'); 
        ((TH1D*)(histMap.find(name))->second)->Fill(phi, EWeight);

        // Fill the correction term histograms
        name = getCorrectionHistName("Beta", eta, PVZ, charge, 'c');
	for(Int_t kk = 0; kk < order; kk++){
	    ((TProfile2D*)(histMap.find(name)->second))->Fill(2 * kk + 1, Day2, cos(kk * phi + phi), EWeight);
	    ((TProfile2D*)(histMap.find(name)->second))->Fill(2 * kk + 2, Day2, sin(kk * phi + phi), EWeight);
	}
    }

}

void StProtonLaGammaMaker::computePhiWeightsHelper(Int_t i, Int_t ii, Int_t iii, Char_t particle, std::map<std::string, TH1*>& histMap){
    std::string Eta = ((i == 0)? "FF" : "RF");
    std::string PVZ = ((ii == 0)? "PVZPos" : "PVZNeg");
    std::string Charge = ((iii == 0)? "ChPos" : "ChNeg");

    char histname[100];
   
    //std::cout << "particle is " << particle << std::endl;
    if(particle == 'a' || particle == 'b'){
        std::string particle_name = (particle == 'a')? "Alpha":"Beta";
	sprintf(histname, "h1d_before_Corrections_%s_%s_%s_%sPhi", Eta.c_str(), PVZ.c_str(), Charge.c_str(), particle_name.c_str());
    }
    else if(particle == 'p'){
	sprintf(histname, "h1d_before_Corrections_%s_%s_%s_PrimaryTrkPhi", Eta.c_str(), PVZ.c_str(), Charge.c_str());
    }
    else
        return;


    std::string histname_string(histname);
    //std::cout << histname_string << std::endl;
    TH1D* hist = (TH1D*)histMap[histname_string];
    Float_t phi_mean = ((TH1D*)histMap[histname_string])->GetSum() / (Float_t)phiBins;

    switch(particle){
	case 'a':
	    for(Int_t j = 0; j < phiBins; j++)
		m_AlphaPhiWeight[i][ii][iii][j] = (hist->GetBinContent(j + 1) != 0.)? phi_mean / (hist->GetBinContent(j + 1)) : 1;
	    break;
	case 'b':
	    for(Int_t j = 0; j < phiBins; j++)
		m_BetaPhiWeight[i][ii][iii][j] = (hist->GetBinContent(j + 1) != 0.)?  phi_mean / (hist->GetBinContent(j + 1)) : 1;
	    break;
	case 'p':
	    for(Int_t j = 0; j < phiBins; j++)
		m_PrimaryTracksPhiWeight[i][ii][iii][j] = (hist->GetBinContent(j + 1) != 0.)? phi_mean / (hist->GetBinContent(j + 1)) : 1;
	    break;
	default:
	    break;
    }  
    return;
}
