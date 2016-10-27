#include <iostream>
#include <set>
#include "../include/StPriTrkInfo.h"
#include "StPriTrkCuts.h"
ClassImp(StPriTrkCuts)

bool StPriTrkCuts::PassFlagCuts(const StPriTrkInfo& p) const{
    return p.flag > FlagLowerBound(); 
}

bool StPriTrkCuts::PassEtaCuts(const StPriTrkInfo& p) const{
    return EtaLowerBound() < p.eta && p.eta < EtaUpperBound();
}

bool StPriTrkCuts::PassPtCuts(const StPriTrkInfo& p) const{
    return PtLowerBound() < p.pt && p.pt < PtUpperBound();
}

bool StPriTrkCuts::PassNHitsCuts(const StPriTrkInfo& p) const{
    //return (p.nhits > NHitsLowerBound()) && (p.nhitsratio > NHitsRatioLowerBound());
}

bool StPriTrkCuts::PassChargeCuts(const StPriTrkInfo& p) const{
    //for(set<int>::const_iterator iter = ChargeCuts().begin(); iter != ChargeCuts().end(); iter++){
     //   if(p.charge == *iter)
     //      return true;
    //}
    //return false;
}

bool StPriTrkCuts::PassDcaCuts(const StPriTrkInfo& p) const{
    return DcaLowerBound() < p.dcaGlobal && p.dcaGlobal < DcaUpperBound();
    return true; // This is for primary only
}

bool StPriTrkCuts::PassAllCuts(const StPriTrkInfo& p) const{
    return PassEtaCuts(p) 
        && PassPtCuts(p)
        && PassDcaCuts(p)
	//&& PassChargeCuts(p)
	&& PassFlagCuts(p)
	//&& PassNHitsCuts(p)
	&& PassPIDCuts(p); 
}

bool StPriTrkCuts::PassAllCuts(const StV0TrkInfo& v0, int no_dau) const{
    if(no_dau == 1){
	return v0.pt_dau1 > PtLowerBound() && v0.pt_dau1 < PtUpperBound()
	    && v0.dca_dau1 > DcaLowerBound() && v0.dca_dau1 < DcaUpperBound()
	    && v0.eta_dau1 > EtaLowerBound() && v0.eta_dau1 < EtaUpperBound();
    }
    else if(no_dau == 2){
	return v0.pt_dau2 > PtLowerBound() && v0.pt_dau2 < PtUpperBound()
	    && v0.dca_dau2 > DcaLowerBound() && v0.dca_dau2 < DcaUpperBound()
	    && v0.eta_dau2 > EtaLowerBound() && v0.eta_dau2 < EtaUpperBound();
    }
    else{
        return false;
    }
}

void StPriTrkCuts::Dump() const{
    cout << "Cuts are listed below: " << endl;
    cout << "Flag: " << m_FlagLowerBound << " < trk_flag" << endl;
    cout << "NHits: " << m_NHitsLowerBound << " < nhits" << endl;
    cout << "NHitsRatio: " << m_NHitsRatioLowerBound << " < nhitsratio" << endl;
    cout << "Eta: " << m_EtaLowerBound << " < eta < " << m_EtaUpperBound << endl;
    cout << "Pt: " << m_PtLowerBound << " < pt < " << m_PtUpperBound << endl;
    cout << "Dca: " << m_DcaLowerBound << " < dca < " << m_DcaUpperBound << endl;
    cout << "Charge: "; 
    for(set<int>::const_iterator iter = ChargeCuts().begin(); iter != ChargeCuts().end(); iter++)
        cout << *iter << " ";
    cout << endl;
    
    DumpPIDCuts();
}
