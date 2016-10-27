#include "StPriPTrkCuts.h"
#include "../include/StPriTrkInfo.h"
#include <iostream>
ClassImp(StPriPTrkCuts)
bool StPriPTrkCuts::PassPIDCuts(const StPriTrkInfo& p) const{
    bool nsigma_flag = (p.nsigma_p < NSigmaPUpperBound()) && (p.nsigma_p > NSigmaPLowerBound());

    bool tof_flag = true;
    if(p.beta > 0){
	if(p.m2 < M2UpperBound() && p.m2 > M2LowerBound()) 
	    tof_flag = true;
	else
	    tof_flag = false;
    }
    else{
	tof_flag = true;
    }


    return (nsigma_flag && tof_flag);
}

void StPriPTrkCuts::DumpPIDCuts() const{
    cout << "NSigmaP: " << NSigmaPLowerBound() << " < nsigma_p < " << NSigmaPUpperBound() << endl;
    cout << "M2(if momentum > 0.65): " << M2LowerBound() << " < m2_p < " << M2UpperBound() << endl;
}
