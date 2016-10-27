/*
for local analysis only
*/
#ifndef STEVTCUTS_H
#define STEVTCUTS_H
//#include "StMuDSTMaker/COMMON/StMuDst.h"
//#include "StMuDSTMaker/COMMON/StMuEvent.h"
//#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "../StEvtInfo/StEvtInfo.h"
#include <iostream> 
#include <set>
class StEvtCuts{
public:
    StEvtCuts(): m_VzUpperBound(40), m_VzLowerBound(-40), m_VrUpperBound(9999.0),
                 m_VrLowerBound(0), m_RunUpperBound(999999999), m_RunLowerBound(-1){}
    virtual ~StEvtCuts(){}
    void LoadBadRunsDay(const set<int>& badruns) {m_BadRunsDay = badruns; }
    void LoadBadRunsDay2(const set<int>& badruns){ m_BadRunsDay2 = badruns; }
    void LoadBadRunsDay3(const set<int>& badruns) { m_BadRunsDay3 = badruns; }
    void SetCentralities(const set<int>& cent){ m_SetCentralities = cent; }

    void SetVzUpperBound(double ub){ m_VzUpperBound = ub; }
    void SetVzLowerBound(double lb){ m_VzLowerBound = lb; }
    void SetVrUpperBound(double ub){ m_VrUpperBound = ub; }
    void SetVrLowerBound(double lb){ m_VrLowerBound = lb; }
    void SetRunLowerBound(int lb){ m_RunLowerBound = lb; }
    void SetRunUpperBound(int ub){ m_RunUpperBound = ub; }

    double VzLowerBound() const { return m_VzLowerBound; }
    double VzUpperBound() const { return m_VzUpperBound; }
    int RunUpperBound() const { return m_RunUpperBound; }
    int RunLowerBound() const { return m_RunLowerBound; }

    virtual bool PassCentralityCuts(const StEvtInfo& evt) const;
    virtual bool PassNullVtxCuts(const StEvtInfo& evt) const;
    virtual bool PassVzCuts(const StEvtInfo& evt) const;
    virtual bool PassVrCuts(const StEvtInfo& evt, double x = 0, double y = 0) const;
    virtual bool PassRunNoCuts(const StEvtInfo& evt) const;
    virtual bool PassAllCuts(const StEvtInfo& evtInfo, int testNTracks = 0, double testVz = 0, double x = 0, double y = 0) const; 
    virtual void Dump() const;

private:
    set<int> m_SetCentralities; 
    double m_VzUpperBound;
    double m_VzLowerBound;
    double m_VrUpperBound;
    double m_VrLowerBound;
    int m_RunUpperBound;
    int m_RunLowerBound;
    set<int> m_BadRunsDay;
    set<int> m_BadRunsDay2;
    set<int> m_BadRunsDay3;

    ClassDef(StEvtCuts,1)
};
#endif
