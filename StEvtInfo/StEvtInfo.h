/*
This class is only used in local analysis auau39GeV 2010
*/
#ifndef STEVTINFO_H
#define STEVTINFO_H
#include <iostream>
#include "../include/StPriTrkInfo.h"
#include "../include/StV0TrkInfo.h"
#include "../StRefMultCorr/StRefMultCorr.h"
#include "TChain.h"
using namespace std;
class StEvtInfo{
public:
    StEvtInfo(){}
    StEvtInfo(TChain* chain, StRefMultCorr& refmultCorrUtil):m_AlphaCharge(1.0), m_BetaCharge(1.0){
       m_RunId = chain->GetLeaf("mRunId")->GetValue(0); 
       m_EventId = chain->GetLeaf("mEventId")->GetValue(0);
       m_RefMult = chain->GetLeaf("mRefMult")->GetValue(0);
       m_Bz = chain->GetLeaf("mBz")->GetValue(0);
       m_Vx = chain->GetLeaf("mPrimaryVertexX")->GetValue(0);
       m_Vy = chain->GetLeaf("mPrimaryVertexY")->GetValue(0);
       m_Vz = chain->GetLeaf("mPrimaryVertexZ")->GetValue(0);
       m_NPTracks = chain->GetLeaf("mNoTracks")->GetValue(0);
       m_NLambdas = chain->GetLeaf("mNoLambdas")->GetValue(0);
       m_ZdcCoinRate = chain->GetLeaf("mZDCCoin")->GetValue(0); // TODO
       m_Day = int((m_RunId - 12000000) / 1000); // TODO:
       m_Day2 = int((m_RunId - 12000000) / 10); // TODO:
       m_Day3 = int((m_RunId - 12000000) / 1); // TODO:

// Add Primary Tracks, used to reconstruct event planes
       for(int i = 0; i < m_NPTracks; ++i){
           StPriTrkInfo trk(chain, i); 
	   AddPrimaryTrk(trk);
       }

// Add Alpha Tracks no?

// Add Beta Tracks
       for(int i = 0; i < m_NLambdas; ++i){
           StV0TrkInfo trk(chain, m_BetaCharge, i); // TODO: Lambda particles, the beta charge is the baryonic charge
           AddBetaTrk(trk);
       }

       refmultCorrUtil.init(m_RunId);
       refmultCorrUtil.initEvent(m_RefMult, m_Vz, m_ZdcCoinRate);

       m_Centrality = refmultCorrUtil.getCentralityBin9() + 1;
       m_EWeight = refmultCorrUtil.getWeight();
       m_badrun = (refmultCorrUtil.isBadRun(m_RunId))? true:false;
    }

    // Getters...
    int Centrality() const { return m_Centrality; }
    double EWeight() const { return m_EWeight; }
    long RunId() const { return m_RunId; }
    int EventId() const { return m_EventId; }
    int RefMult() const { return m_RefMult; }
    double Vx() const { return m_Vx; }
    double Vy() const { return m_Vy; }
    double Vz() const { return m_Vz; }
    int NPTracks() const { return m_NPTracks; }
    int NLambdas() const { return m_NLambdas; }
    int ZDC() const { return m_ZdcCoinRate; }
    int Day() const { return m_Day; }
    int Day2() const { return m_Day2; }
    int Day3() const { return m_Day3; }
    bool IsBadRunRef() const { return m_badrun; }

    // Setters...
    void SetCentrality(int cent) { m_Centrality = cent; }
    void SetEWeight(double weight) { m_EWeight = weight; }
    void SetAlphaCharge(float charge) { m_AlphaCharge = charge; }
    void SetBetaCharge(float charge) { m_BetaCharge = charge; }

    virtual ~StEvtInfo(){}
    virtual const vector<StPriTrkInfo>& VecPriTrks() const { return m_VecPriTrks; }
    virtual const vector<StPriTrkInfo>& VecAlphaTrks() const { return m_VecAlphaTrks; }
    virtual const vector<StV0TrkInfo>& VecBetaTrks() const { return m_VecBetaTrks; } // TODO:
    virtual void AddPrimaryTrk(const StPriTrkInfo& trk);
    virtual void AddAlphaTrk(const StPriTrkInfo& trk);
    virtual void AddBetaTrk(const StV0TrkInfo& trk); // TODO: type
private: 
    int m_Centrality;
    double m_EWeight;
    long m_RunId;
    int m_EventId;
    int m_RefMult;
    int m_Bz;
    double m_Vx;
    double m_Vy;
    double m_Vz; 
    int m_NPTracks;
    int m_NLambdas;
    int m_ZdcCoinRate; 
    int m_Day;
    int m_Day2;
    int m_Day3;
    bool m_badrun;
    float m_AlphaCharge;
    float m_BetaCharge;
    vector<StPriTrkInfo> m_VecPriTrks;
    vector<StPriTrkInfo> m_VecAlphaTrks;
    vector<StV0TrkInfo> m_VecBetaTrks;
    ClassDef(StEvtInfo, 1)
};
#endif
