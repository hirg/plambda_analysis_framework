#include "StPriPiTrkCuts.h"
#include "../include/StPriTrkInfo.h"
#include <iostream>
ClassImp(StPriPiTrkCuts)
bool StPriPiTrkCuts::PassPIDCuts(const StPriTrkInfo& p) const{
    bool nsigma_flag = (p.nsigma_pi < NSigmaPiUpperBound()) && (p.nsigma_pi > NSigmaPiLowerBound());

    bool tof_flag = true;
    if(p.p < 0.65){
	if(p.beta > 0)
	    if(p.m2 < M2UpperBound() && p.m2 > M2LowerBound()) 
		tof_flag = true;
	    else
		tof_flag = false;
	else
	    tof_flag = true;
    }
    else{
        if(p.m2 < M2UpperBound() && p.m2 > M2LowerBound())
            tof_flag = true;
        else
            tof_flag = false;
    }

    return (nsigma_flag && tof_flag);
}

void StPriPiTrkCuts::DumpPIDCuts() const{
    cout << "NSigmaPi: " << NSigmaPiLowerBound() << " < nsigma_pi < " << NSigmaPiUpperBound() << endl;
    cout << "M2: " << M2LowerBound() << " < m2_p < " << M2UpperBound() << endl;
}
