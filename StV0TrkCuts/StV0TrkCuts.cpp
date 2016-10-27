#include "StV0TrkCuts.h"
#include "../include/constants.h"
#include "TLorentzVector.h"
#include "../include/StV0TrkInfo.h"
ClassImp(StV0TrkCuts)

bool StV0TrkCuts::PassAllCuts(const StV0TrkInfo& v0) const{
    bool flag_pt = (m_PtLowerBound < v0.pt && v0.pt < m_PtUpperBound);
    bool flag_eta = (m_EtaLowerBound < v0.eta && v0.eta < m_EtaUpperBound);
    bool flag_dca = (m_DcaLowerBound < v0.dcaGlobal && v0.dcaGlobal < m_DcaUpperBound);
    bool flag_declen = (m_DecLenLowerBound < v0.declen && v0.declen < m_DecLenUpperBound);
    bool flag_dca1to2 = (m_Dca1to2LowerBound < v0.dca1to2 && v0.dca1to2 < m_Dca1to2UpperBound);
    bool flag_dau1nsigma = (m_Dau1NSigmaLowerBound < v0.nsigma_dau1 && v0.nsigma_dau1 < m_Dau1NSigmaUpperBound);  
    bool flag_dau2nsigma = (m_Dau2NSigmaLowerBound < v0.nsigma_dau2 && v0.nsigma_dau2 < m_Dau2NSigmaUpperBound);
    bool flag_dau1dca = (m_Dau1DcaLowerBound < v0.dca_dau1 && v0.dca_dau2 < m_Dau1DcaUpperBound);
    bool flag_dau2dca = (m_Dau2DcaLowerBound < v0.dca_dau2 && v0.dca_dau2 < m_Dau2DcaUpperBound);
    bool flag_mass = (m_MassLowerBound < v0.mass && v0.mass < m_MassUpperBound);
    bool flag_lowptla = true;
    if((v0.pt < 0.6 && v0.dca_dau1 < 0.7)  || (v0.pt < 0.6 && v0.dca_dau2 < 2.5))
        flag_lowptla = false;

    return flag_pt && flag_eta && flag_dca & flag_declen && flag_dca1to2 && flag_dau1nsigma && flag_dau2nsigma
                   && flag_dau1dca && flag_dau2dca && flag_mass && flag_lowptla;
}
