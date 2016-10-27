#ifndef STV0TrkCUTS_H
#define STV0TrkCUTS_H
#include "TLorentzVector.h"
#include "../include/StV0TrkInfo.h"
#include <iostream>
class StV0TrkCuts{
public:
    StV0TrkCuts():m_PtLowerBound(.5), m_PtUpperBound(5.0), 
                    m_EtaLowerBound(-1.0), m_EtaUpperBound(1.0),
                    m_DcaLowerBound(-1), m_DcaUpperBound(0.6), 
		    m_Dau1NSigmaLowerBound(-3.0), m_Dau1NSigmaUpperBound(3.0), 
		    m_Dau2NSigmaLowerBound(-3.0), m_Dau2NSigmaUpperBound(3.0), 
		    m_DecLenLowerBound(6.0), m_DecLenUpperBound(999999.0), 
		    m_Dau1DcaLowerBound(0.6), m_Dau1DcaUpperBound(999999.0), 
		    m_Dau2DcaLowerBound(1.8), m_Dau2DcaUpperBound(999999.0),
		    m_Dca1to2LowerBound(-1), m_Dca1to2UpperBound(0.7),
		    m_MassLowerBound(1.115683 - 0.004), m_MassUpperBound(1.115683 + 0.004){}

    virtual ~StV0TrkCuts(){}

    void SetPtLowerBound(double lb) { m_PtLowerBound = lb; }
    void SetPtUpperBound(double ub) { m_PtUpperBound = ub; }
    void SetEtaLowerBound(double lb) { m_EtaLowerBound = lb; }
    void SetEtaUpperBound(double ub) { m_EtaUpperBound = ub; }
    void SetDcaLowerBound(double lb) { m_DcaLowerBound = lb; }
    void SetDcaUpperBound(double ub) { m_DcaUpperBound = ub; }
    void SetDau1NSigmaLowerBound(double lb) { m_Dau1NSigmaLowerBound = lb; }
    void SetDau1NSigmaUpperBound(double ub) { m_Dau1NSigmaUpperBound = ub; }
    void SetDau2NSigmaLowerBound(double lb) { m_Dau2NSigmaLowerBound = lb; } 
    void SetDau2NSigmaUpperBound(double ub) { m_Dau2NSigmaUpperBound = ub; }
    void SetDecLenLowerBound(double lb) { m_DecLenLowerBound = lb; }
    void SetDecLenUpperBound(double ub) { m_DecLenUpperBound = ub; }
    void SetDau1DcaLowerBound(double lb) { m_Dau1DcaLowerBound = lb; }
    void SetDau1DcaUpperBound(double ub) { m_Dau1DcaUpperBound = ub ;}
    void SetDau2DcaLowerBound(double lb) { m_Dau2DcaLowerBound = lb; }
    void SetDau2DcaUpperBound(double ub) { m_Dau2DcaUpperBound = ub ;}
    void SetDca1to2LowerBound(double lb) { m_Dca1to2LowerBound = lb; }
    void SetDca1to2UpperBound(double ub) { m_Dca1to2UpperBound = ub; }
    void SetMassLowerBound(double lb){ m_MassLowerBound = lb; }
    void SetMassUpperBound(double ub){ m_MassUpperBound = ub; }
    void SetDipAngleLowerBound(double lb){ m_DipAngleLowerBound = lb; }
    void SetDipAngleUpperBound(double ub){ m_DipAngleUpperBound = ub; }
    void SetRapLowerBound(double lb){ m_RapLowerBound = lb; }
    void SetRapUpperBound(double ub){ m_RapUpperBound = ub; }

    double RapLowerBound() const { return m_RapLowerBound; }
    double RapUpperBound() const { return m_RapUpperBound; }
    double DipAngleLowerBound() const { return m_DipAngleLowerBound; }
    double DipAngleUpperBound() const { return m_DipAngleUpperBound; }
    double MassLowerBound() const { return m_MassLowerBound; }
    double MassUpperBound() const { return m_MassUpperBound; }

//    virtual bool PassDipAngleCuts(const StV0TrkInfo& v0);

//    virtual bool PassIdCuts(const StV0TrkInfo& v0);

//    virtual bool PassRapCuts(const StV0TrkInfo& v0);

//    virtual bool PassMassCuts(const StV0TrkInfo& v0);

    virtual bool PassAllCuts(const StV0TrkInfo& v0) const;
    virtual void Dump() const {
//        cout << m_PtLowerBound << " < pt < " << m_PtUpperBound << endl; 
    }
private:
    double m_PtLowerBound;
    double m_PtUpperBound;
    double m_EtaLowerBound;
    double m_EtaUpperBound;
    double m_RapLowerBound;
    double m_RapUpperBound;
    double m_DcaLowerBound;
    double m_DcaUpperBound;
    double m_Dau1NSigmaLowerBound;
    double m_Dau1NSigmaUpperBound;
    double m_Dau2NSigmaLowerBound;
    double m_Dau2NSigmaUpperBound;
    double m_DecLenLowerBound;
    double m_DecLenUpperBound;
    double m_Dau1DcaLowerBound;
    double m_Dau1DcaUpperBound;
    double m_Dau2DcaLowerBound;
    double m_Dau2DcaUpperBound;
    double m_Dca1to2LowerBound;
    double m_Dca1to2UpperBound;
    double m_MassLowerBound;
    double m_MassUpperBound;
    double m_DipAngleLowerBound;
    double m_DipAngleUpperBound;

    //virtual TLorentzVector FormPair(const StPriTrkInfo& trk1, const StPriTrkInfo& trk2);

    //virtual double DipAngle(const StPriTrkInfo& trk1, const StPriTrkInfo& trk2);
    ClassDef(StV0TrkCuts,1)
};
#endif
