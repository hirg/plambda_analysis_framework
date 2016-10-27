#include "StPriKTrkCuts.h"
#include "../include/StPriTrkInfo.h"
#include <iostream>
ClassImp(StPriKTrkCuts)
bool StPriKTrkCuts::PassPIDCuts(const StPriTrkInfo& p) const{
    bool nsigma_flag = (p.nsigma_k < NSigmaKUpperBound()) && (p.nsigma_k > NSigmaKLowerBound());

    bool tof_flag = true;
    if(p.p < 0.65){
	if(p.beta > 0){
	    if(p.m2 < M2UpperBound() && p.m2 > M2LowerBound()) 
		tof_flag = true;
	    else
		tof_flag = false;
	}
	else{
	    tof_flag = true;
	}
    }
    else{
        if(p.m2 < M2UpperBound() && p.m2 > M2LowerBound())
            tof_flag = true;
        else
            tof_flag = false;
    }

    return nsigma_flag && tof_flag;
}

void StPriKTrkCuts::DumpPIDCuts() const{
    cout << "NSigmaK: " << NSigmaKLowerBound() << " < nsigma_k < " << NSigmaKUpperBound() << endl;
    cout << "M2: " << M2LowerBound() << " < m2_k < " << M2UpperBound() << endl;
}
