#include "StEvtCuts.h"
//#include "StMuDSTMaker/COMMON/StMuDst.h"
//#include "StMuDSTMaker/COMMON/StMuEvent.h"
//#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "../StEvtInfo/StEvtInfo.h"
#include <iostream>
ClassImp(StEvtCuts)
bool StEvtCuts::PassCentralityCuts(const StEvtInfo& evtInfo) const {
    for(set<int>::iterator iter = m_SetCentralities.begin(); iter != m_SetCentralities.end(); iter++){
        if(evtInfo.Centrality() == *iter)
            return true;
    }
    return false; 
}

bool StEvtCuts::PassNullVtxCuts(const StEvtInfo& evtInfo) const {
    if(fabs(evtInfo.Vx()) < 1e-5 
	    && fabs(evtInfo.Vy()) < 1e-5 
	    && fabs(evtInfo.Vz()) < 1e-5 )
	return false;
    else
	return true;
}

bool StEvtCuts::PassVzCuts(const StEvtInfo& evtInfo) const {
    float vz = evtInfo.Vz();
    return (m_VzLowerBound < vz && vz < m_VzUpperBound);
}

bool StEvtCuts::PassVrCuts(const StEvtInfo& evtInfo, double x, double y) const {
    float vr = sqrt(pow(evtInfo.Vx() - x, 2) + pow(evtInfo.Vy() - y, 2));
    return (m_VrLowerBound < vr && vr < m_VrUpperBound);
}

bool StEvtCuts::PassRunNoCuts(const StEvtInfo& evtInfo) const {//TODO
    for(set<int>::const_iterator it = m_BadRunsDay.begin(); it != m_BadRunsDay.end(); it++){
        if(evtInfo.Day() == *it){
            return false;
	}
    }

    for(set<int>::const_iterator it = m_BadRunsDay2.begin(); it != m_BadRunsDay2.end(); it++){
        if(evtInfo.Day2() == *it){
            return false;
	}
    }

    for(set<int>::const_iterator it = m_BadRunsDay3.begin(); it != m_BadRunsDay3.end(); it++){
        if(evtInfo.Day3() == *it){
            return false;
	}
    }

    if(evtInfo.RunId() < RunLowerBound() || evtInfo.RunId() > RunUpperBound()){
//        cout << "runnumber = " << evtInfo.RunId() << endl;
//	cout << "badrun out of range" << endl;
        return false;
    }

    if(evtInfo.IsBadRunRef())
        return false;

    return true;
}

bool StEvtCuts::PassAllCuts(const StEvtInfo& evtInfo, int testNTracks, double testVz, double x, double y) const{
    return (PassCentralityCuts(evtInfo)
	    && PassVzCuts(evtInfo)
	    && PassVrCuts(evtInfo, x, y)
            && PassRunNoCuts(evtInfo));
}

void StEvtCuts::Dump() const{
    cout << "Event cuts are listed below: " << endl;
    cout << "    NullVtxCuts: " << endl;
    //cout << "    Triggers:";
    //for(int i = 0; i < m_TriggerIds.size(); i++) cout << " " << m_TriggerIds[i];
    //cout << ";" << endl;;
    cout << "    Vz Cuts: " << m_VzLowerBound << " < vz < " << m_VzUpperBound << endl;  
    cout << "    Vr Cuts: " << m_VrLowerBound << " < vr < " << m_VrUpperBound << endl;  
}
