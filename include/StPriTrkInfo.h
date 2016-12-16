#ifndef STPRITRKINFO_H
#define STPRITRKINFO_H
//#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "TChain.h"
#include "TVector2.h"
#include "TLeaf.h"
#include "constants.h"
struct StPriTrkInfo{// nHits
    StPriTrkInfo(){}
    StPriTrkInfo(TChain* chain, int idx){
	id = chain->GetLeaf("fV0s.mTrackIdV0")->GetValue(idx);
	flag = chain->GetLeaf("fV0s.mTrackFlagV0")->GetValue(idx);
	eta = chain->GetLeaf("fV0s.mEtaV0")->GetValue(idx);
	pt = chain->GetLeaf("fV0s.mPtprimaryV0")->GetValue(idx);
	dcaGlobal = chain->GetLeaf("fV0s.mDCAglobalV0")->GetValue(idx);
	charge = chain->GetLeaf("fV0s.mChargeV0")->GetValue(idx);
	nsigma_e = chain->GetLeaf("fV0s.mnSigmaEV0")->GetValue(idx);
	nsigma_p = chain->GetLeaf("fV0s.mnSigmaPV0")->GetValue(idx);
	nsigma_pi = chain->GetLeaf("fV0s.mnSigmaPiV0")->GetValue(idx);
	nsigma_k = chain->GetLeaf("fV0s.mnSigmaKV0")->GetValue(idx);
	p = chain->GetLeaf("fV0s.mPV0")->GetValue(idx);
	double tof = chain->GetLeaf("fV0s.mTofV0")->GetValue(idx);
	double path = chain->GetLeaf("fV0s.mPathlenV0")->GetValue(idx);
        beta = path / tof / 29.9792458;// TODO: Need Check;
	//cout << "beta: " << beta << endl;
	m2 = p * p * (1.0 / beta / beta - 1.0);

	phi = chain->GetLeaf("fV0s.mPhiV0")->GetValue(idx);
        //TVector3 vec3_pritrk(chain->GetLeaf("fV0s.mPxprimaryV0")->GetValue(idx),
        //                    chain->GetLeaf("fV0s.mPyprimaryV0")->GetValue(idx),
        //                    chain->GetLeaf("fV0s.mPzprimaryV0")->GetValue(idx));
        //phi = vec3_pritrk.Phi(); 
    }

    short id;
    short flag;
    double eta;
    double pt;
    double dcaGlobal;
    short charge;
    double nsigma_e; 
    double nsigma_p;
    double nsigma_pi;
    double nsigma_k;
    double p;
    float beta;
    double m2;
    double phi;
};
#endif
