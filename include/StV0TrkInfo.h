#ifndef STV0TRKINFO_H
#define STV0TRKINFO_H
//#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "TChain.h"
#include "TVector2.h"
#include "TLeaf.h"
#include "constants.h"
struct StV0TrkInfo{// nHits
    StV0TrkInfo(){}
    StV0TrkInfo(TChain* chain, float charge, int idx){
	pt = chain->GetLeaf("fW0s.mPtglobalW0")->GetValue(idx);
        px = chain->GetLeaf("fW0s.mPxglobalW0")->GetValue(idx);
        py = chain->GetLeaf("fW0s.mPyglobalW0")->GetValue(idx);
        pz = chain->GetLeaf("fW0s.mPzglobalW0")->GetValue(idx);
	eta = chain->GetLeaf("fW0s.mEtaW0")->GetValue(idx);
        mass = chain->GetLeaf("fW0s.mMassW0")->GetValue(idx);    
        nsigma_dau1 = chain->GetLeaf("fW0s.mDau1nSigmaW0")->GetValue(idx);
        nsigma_dau2 = chain->GetLeaf("fW0s.mDau2nSigmaW0")->GetValue(idx);
        id_dau1 = chain->GetLeaf("fW0s.mDau1idW0")->GetValue(idx);
        id_dau2 = chain->GetLeaf("fW0s.mDau2idW0")->GetValue(idx);
        dcaGlobal = chain->GetLeaf("fW0s.mDcaW0")->GetValue(idx);
        declen = chain->GetLeaf("fW0s.mDecaylenW0")->GetValue(idx);
	charge_dau1 = charge;
	charge_dau2 = -charge;
	pt_dau1 = chain->GetLeaf("fW0s.mDau1ptW0")->GetValue(idx);
	pt_dau2 = chain->GetLeaf("fW0s.mDau2ptW0")->GetValue(idx);
	eta_dau1 = chain->GetLeaf("fW0s.mDau1etaW0")->GetValue(idx);
	eta_dau2 = chain->GetLeaf("fW0s.mDau2etaW0")->GetValue(idx);
        dca_dau1 = chain->GetLeaf("fW0s.mDau1dcaW0")->GetValue(idx); 
        dca_dau2 = chain->GetLeaf("fW0s.mDau2dcaW0")->GetValue(idx); 
	prmatch_dau1 = chain->GetLeaf("fW0s.mDau1PrMatchW0")->GetValue(idx);
	prmatch_dau2 = chain->GetLeaf("fW0s.mDau2PrMatchW0")->GetValue(idx);
        dca1to2 = chain->GetLeaf("fW0s.mDca1to2W0")->GetValue(idx);

        TVector2 phi2;
        phi2.Set(px, py); 
        phi = phi2.Phi();
	SetCharge(charge); // TODO: For Lambda Only
	if(phi > PI) phi -= 2 * PI;
	if(phi < -PI) phi += 2 * PI;

        TVector2 vec2_dau1(chain->GetLeaf("fW0s.mDau1pxW0")->GetValue(idx),
	                   chain->GetLeaf("fW0s.mDau1pyW0")->GetValue(idx));
	phi_dau1 = vec2_dau1.Phi();
	if(phi_dau1 > PI) phi_dau1 -= 2 * PI;
	if(phi_dau1 < -PI) phi_dau1 += 2 * PI;

        TVector2 vec2_dau2(chain->GetLeaf("fW0s.mDau2pxW0")->GetValue(idx),
	                   chain->GetLeaf("fW0s.mDau2pyW0")->GetValue(idx));
	phi_dau2 = vec2_dau2.Phi();
	if(phi_dau2 > PI) phi_dau2 -= 2 * PI;
	if(phi_dau2 < -PI) phi_dau2 += 2 * PI;
    }

    void SetCharge(int charge_) { charge = charge_; }

    float pt;
    float px;
    float py; 
    float pz;
    float eta;
    float mass;
    float nsigma_dau1;
    float nsigma_dau2;
    int id_dau1;
    int id_dau2;
    float dcaGlobal;
    float declen;
    float charge_dau1;
    float charge_dau2;
    float pt_dau1;
    float pt_dau2;
    float eta_dau1;
    float eta_dau2;
    float phi_dau1;
    float phi_dau2;
    float dca_dau1;
    float dca_dau2;
    int prmatch_dau1;
    int prmatch_dau2;
    float dca1to2;
    float phi; 
    float charge;
};
#endif
