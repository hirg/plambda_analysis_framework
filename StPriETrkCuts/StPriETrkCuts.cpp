#include "StPriETrkCuts.h"
#include "../include/StPriTrkInfo.h"
#include <iostream>
ClassImp(StPriETrkCuts)
bool StPriETrkCuts::PassPIDCuts(const StPriTrkInfo& p) const{// the cuts used here only apply to low pt electrons < 1GeV
    bool nsigma_flag = (p.nsigma_e < NSigmaEUpperBound()) && (p.nsigma_e > NSigmaELowerBound());

    bool tof_flag = true;
    if(p.p < 0.65){
	if(p.beta > 0){
	    if(fabs(1 - 1 / p.beta) < 0.03) 
		tof_flag = true;
	    else
		tof_flag = false;
	}
	else
	    tof_flag = true;
    }
    else{
        if(fabs(1 - 1 / p.beta) < 0.03)
            tof_flag = true;
        else
            tof_flag = false;
    }

    return nsigma_flag && tof_flag;
}

void StPriETrkCuts::DumpPIDCuts() const{
    cout << "NSigmaE: " << NSigmaELowerBound() << " < nsigma_e < " << NSigmaEUpperBound() << endl;
    cout << " |1 - 1 / beta| < 0.03" << endl;
    //cout << "Beta: " << M2LowerBound() << " < m2_e < " << M2UpperBound() << endl;
}
