#ifndef STCORRELATIONMAKER_H
#define STCORRELATIONMAKER_H
#include "THashTable.h" 
#include "../include/constants.h"
#include "TH1.h"
#include "TChain.h"
#include "TObject.h"
#include "TLeaf.h"
#include "StEffMaker.h"
#include "../StEvtCuts/StEvtCuts.h"
#include "../StEvtInfo/StEvtInfo.h"
#include <map>
#include <set>
#include "TVector2.h"
#include <string>
//This is the base class for event-plane related azimuthal angle correlation analysis
class StCorrelationMaker{
    public:
	StCorrelationMaker(Int_t, Int_t, StEffMaker* effMaker_Alpha, StEffMaker* effMaker_Beta, StEffMaker* effMaker_AlphaBar, StEffMaker* effMaker_BetaBar);
        virtual ~StCorrelationMaker(){}

        // Step Zero: configuration
        void LoadDataSet(const set<string> dataset) { m_DataSet = dataset; } 
	void LoadEvtCuts(StEvtCuts* cuts) { m_EvtCuts = cuts; }
	void SetDay2Boundaries(pair<int, int> bds) { m_Day2LowerBound = bds.first; m_Day2UpperBound = bds.second; }
	void SetDayBoundaries(pair<int, int> bds) { m_DayLowerBound = bds.first; m_DayUpperBound = bds.second; }
	Int_t Init();
        // Step One
	Int_t ComputePhiCorrectionParams();
        virtual void ReconstructEventPlaneWithPhiWeightCorrection();
        virtual void ReconstructShiftedSubEventPlane();
        virtual void ReconstructShiftedFullEventPlaneAndComputeCorrelator();
        void Finish();
        void getWeightsIndex(Int_t*, Float_t, Float_t, Float_t);

        StEffMaker* GetAlphaEffMaker(){ return m_EffMakerAlpha; }
        StEffMaker* GetAlphaBarEffMaker(){ return m_EffMakerAlphaBar; }
        StEffMaker* GetBetaEffMaker(){ return m_EffMakerBeta; }
        StEffMaker* GetBetaBarEffMaker(){ return m_EffMakerBetaBar; }

        int DayLowerBound(){ return m_DayLowerBound; }
        int DayUpperBound(){ return m_DayUpperBound; }
        int Day2LowerBound(){ return m_Day2LowerBound; }
	int Day2UpperBound(){ return m_Day2UpperBound; }
        Int_t GetCentrality(){ return m_Centrality; }

    private:
        set<string> m_DataSet;
        StEffMaker* m_EffMakerAlpha;
        StEffMaker* m_EffMakerAlphaBar;
        StEffMaker* m_EffMakerBeta;
        StEffMaker* m_EffMakerBetaBar;
 
        StEvtCuts* m_EvtCuts;

	virtual Int_t addTrees();
	virtual Int_t initiateHistograms();
        virtual void zeroWeightsHistograms();
	//virtual Int_t initiateComputingParams();

	// No. of events
	Long_t m_NEvents;
	// Centrality
	Int_t m_Centrality;
	// TChain
	TChain* m_Chain;
	// Get azimuthal angle weight to correct possible detector effects in analysis
	virtual void addAlphaHists(std::map<std::string, TH1*>&);
        virtual void addAdditionalAlphaHists(std::map<std::string, TH1*>&){}
	virtual void addBetaHists(std::map<std::string, TH1*>&);
        virtual void addAdditionalBetaHists(std::map<std::string, TH1*>&){}
        virtual void addCorrelationHists(std::map<std::string, TH1*>&){}

        // Compute phi weights
	//virtual void reconstructRawEventPlane(std::map<std::string, TH1*>&, TChain*, LeavesPriTrk leavesPriTrk, Float_t, Float_t, Long_t);
	virtual void fillPrimaryTracksPhiHists(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo){}
	virtual void fillAlphaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo){}
	virtual void fillBetaHists(std::map<std::string, TH1*>& histMap, const StEvtInfo& evtInfo){}

        virtual void computePhiWeights(std::map<std::string, TH1*>& histMap){
            computePrimaryTracksPhiWeights(histMap);
            computeAlphaPhiWeights(histMap);
            computeBetaPhiWeights(histMap);
	}

	virtual void computePrimaryTracksPhiWeights(std::map<std::string, TH1*>& histMap);
        virtual void computeAlphaPhiWeights(std::map<std::string, TH1*>& histMap);
	virtual void computeBetaPhiWeights(std::map<std::string, TH1*>& histMap);

        virtual void computeParticlePhiWeights(std::map<std::string, TH1*>& histMap, Char_t particle);
        virtual void computePhiWeightsHelper(Int_t, Int_t, Int_t, Char_t particle, std::map<std::string, TH1*>&){}

        // Event plane correction I
        virtual void reconstructSubEventPlaneWithPhiWeightHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo){}
        // Event plane correction II
        virtual void reconstructShiftedSubEventPlaneHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo){}
        // Event plane correction III
        virtual void reconstructShiftedFullEventPlaneHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo, TVector2& shiftedepphi_full){}
        virtual void computeCorrelatorsHelper(std::map<std::string, TH1*>& hsitMap, const StEvtInfo& evtInfo, TVector2 shiftedepphi_full){}
  

        // Get event plane weight
	//Int_t computeRawSubEpCorrectedByPhiWeight(){}

        // Histograms collections
        std::map<std::string, TH1*> m_HistsCollection;
        //std::unordered_map<std::string, TObject*> m_ParticlePhiHistsCollection;	
	//std::unoredered_map<std::string, TObject*> m_EventPlaneHistsCollection;

	// Output files
	TFile* m_WeightFile;
	TFile* m_ResultFile;

	// Event number to processed
	Long_t m_EventNumberToProcess;

        int m_DayLowerBound;
	int m_DayUpperBound;
        int m_Day2LowerBound;
        int m_Day2UpperBound;

        ClassDef(StCorrelationMaker, 1)
};

inline
void StCorrelationMaker::getWeightsIndex(Int_t* index, Float_t Eta, Float_t PVZ, Float_t Charge){
    index[0] = (Eta > 0)? 0 : 1;
    index[1] = (PVZ > 0)? 0 : 1;
    index[2] = (Charge > 0)? 0 : 1;
}

inline
void StCorrelationMaker::computeParticlePhiWeights(std::map<std::string, TH1*>& histMap, Char_t particle){
    for(Int_t i = 0; i < 2; i++){
        for(Int_t ii = 0; ii < 2; ii++){
            for(Int_t iii = 0; iii < 2; iii++)
                computePhiWeightsHelper(i, ii, iii, particle, histMap);
        }
    }
}

inline
void StCorrelationMaker::computePrimaryTracksPhiWeights(std::map<std::string, TH1*>& histMap){
    computeParticlePhiWeights(histMap, 'p');
}

inline
void StCorrelationMaker::computeAlphaPhiWeights(std::map<std::string, TH1*>& histMap){
    computeParticlePhiWeights(histMap, 'a');
}

inline
void StCorrelationMaker::computeBetaPhiWeights(std::map<std::string, TH1*>& histMap){
    computeParticlePhiWeights(histMap, 'b');
}
#endif
