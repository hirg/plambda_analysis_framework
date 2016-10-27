#ifndef STPRITRKCUTS_H
#define STPRITRKCUTS_H
#include <iostream>
#include <set>
#include "../include/StPriTrkInfo.h" 
#include "../include/constants.h"
#include "../include/StV0TrkInfo.h"
/*
 * Could be used for phi and strong decay particle reconstruction to select primary tracks we want
 * Usage: users have to initialize every cuts manually! E.g., by doing like `cut->SetEtaUpperBound(1.0);`
 */
class StPriTrkCuts{
public:
    StPriTrkCuts():
	m_FlagUpperBound(999), m_FlagLowerBound(0),
	m_NHitsLowerBound(15), m_NHitsRatioLowerBound(0.52),
	m_EtaUpperBound(1.0), m_EtaLowerBound(-1.0), 
	m_PtUpperBound(10.0), m_PtLowerBound(0.15), 
	m_DcaUpperBound(2.0), m_DcaLowerBound(-1.0),//TODO
	m_NSigmaEUpperBound(3.0), m_NSigmaELowerBound(-3.0),
	m_NSigmaPUpperBound(3.0), m_NSigmaPLowerBound(-3.0),
	m_NSigmaPiUpperBound(3.0), m_NSigmaPiLowerBound(-3.0),
	m_NSigmaKUpperBound(3.0), m_NSigmaKLowerBound(-3.0){}

    virtual ~StPriTrkCuts(){}

    void SetFlagLowerBound(short lb){ m_FlagLowerBound = lb; }
    void SetFlagUpperBound(short ub){ m_FlagUpperBound = ub; }
    void SetNHitsLowerBound(int lb){ m_NHitsLowerBound = lb; }
    void SetNHitsRatioLowerBound(double lb){ m_NHitsRatioLowerBound = lb; }
    void SetEtaUpperBound(double ub){ m_EtaUpperBound = ub; }
    void SetEtaLowerBound(double lb){ m_EtaLowerBound = lb; }
    void SetPtUpperBound(double ub){ m_PtUpperBound = ub; }
    void SetPtLowerBound(double lb){ m_PtLowerBound = lb; }
    void SetDcaUpperBound(double ub){ m_DcaUpperBound = ub; }
    void SetDcaLowerBound(double lb){ m_DcaLowerBound = lb; }
    void SetNSigmaEUpperBound(double ub){ m_NSigmaEUpperBound = ub; }
    void SetNSigmaELowerBound(double lb){ m_NSigmaELowerBound = lb; }
    void SetNSigmaPUpperBound(double ub){ m_NSigmaPUpperBound = ub; }
    void SetNSigmaPLowerBound(double lb){ m_NSigmaPLowerBound = lb; }
    void SetNSigmaPiUpperBound(double ub){ m_NSigmaPiUpperBound = ub; }
    void SetNSigmaPiLowerBound(double lb){ m_NSigmaPiLowerBound = lb; }
    void SetNSigmaKUpperBound(double ub){ m_NSigmaKUpperBound = ub; }
    void SetNSigmaKLowerBound(double lb){ m_NSigmaKLowerBound = lb; }
    void SetM2UpperBound(double ub){ m_M2UpperBound = ub; }
    void SetM2LowerBound(double lb){ m_M2LowerBound = lb; }
    void SetChargeCuts(set<int> charge){ m_SetChargeCuts = charge; }

    short FlagLowerBound() const { return m_FlagLowerBound; }
    short FlagUpperBound() const { return m_FlagUpperBound; }
    double EtaLowerBound() const { return m_EtaLowerBound; }
    double EtaUpperBound() const { return m_EtaUpperBound; }
    double PtLowerBound() const { return m_PtLowerBound; }
    double PtUpperBound() const { return m_PtUpperBound; }
    int NHitsLowerBound() const { return m_NHitsLowerBound; }
    double NHitsRatioLowerBound() const { return m_NHitsRatioLowerBound; }
    double DcaLowerBound() const { return m_DcaLowerBound; }
    double DcaUpperBound() const { return m_DcaUpperBound; }
    double NSigmaEUpperBound() const { return m_NSigmaEUpperBound; }
    double NSigmaELowerBound() const { return m_NSigmaELowerBound; }
    double NSigmaPUpperBound() const { return m_NSigmaPUpperBound; }
    double NSigmaPLowerBound() const { return m_NSigmaPLowerBound; }
    double NSigmaPiUpperBound() const { return m_NSigmaPiUpperBound; }
    double NSigmaPiLowerBound() const { return m_NSigmaPiLowerBound; }
    double NSigmaKUpperBound() const { return m_NSigmaKUpperBound; }
    double NSigmaKLowerBound() const { return m_NSigmaKLowerBound; }
    double M2UpperBound() const { return m_M2UpperBound; }
    double M2LowerBound() const { return m_M2LowerBound; }
    const set<int>& ChargeCuts() const { return m_SetChargeCuts; }

    bool PassFlagCuts(const StPriTrkInfo& p) const;

    bool PassEtaCuts(const StPriTrkInfo& p) const;

    bool PassPtCuts(const StPriTrkInfo& p) const;

    bool PassNHitsCuts(const StPriTrkInfo& p) const;

    bool PassChargeCuts(const StPriTrkInfo& p) const;

    bool PassDcaCuts(const StPriTrkInfo& p) const;

    bool PassAllCuts(const StPriTrkInfo& p) const;

    bool PassAllCuts(const StV0TrkInfo& v0, int no_dau) const; // Special for this application

    void Dump() const;
private: 
    virtual bool PassPIDCuts(const StPriTrkInfo& p) const = 0;
    virtual void DumpPIDCuts() const = 0;

    short m_FlagUpperBound;
    short m_FlagLowerBound;
    int m_NHitsLowerBound;
    double m_NHitsRatioLowerBound;
    double m_EtaUpperBound;
    double m_EtaLowerBound;
    double m_PtUpperBound;
    double m_PtLowerBound;
    double m_DcaUpperBound;
    double m_DcaLowerBound;
    double m_NSigmaEUpperBound;
    double m_NSigmaELowerBound;
    double m_NSigmaPUpperBound;
    double m_NSigmaPLowerBound;
    double m_NSigmaPiUpperBound;
    double m_NSigmaPiLowerBound;
    double m_NSigmaKUpperBound;
    double m_NSigmaKLowerBound;
    double m_M2UpperBound;
    double m_M2LowerBound;
    set<int> m_SetChargeCuts;

    ClassDef(StPriTrkCuts,1)
};
#endif
