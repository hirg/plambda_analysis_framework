#include "StEvtInfo.h"
#include "../include/StPriTrkInfo.h"
#include "../include/StV0TrkInfo.h"
ClassImp(StEvtInfo)
void StEvtInfo::AddPrimaryTrk(const StPriTrkInfo& trk){ m_VecPriTrks.push_back(trk); }
void StEvtInfo::AddAlphaTrk(const StPriTrkInfo& trk){ m_VecAlphaTrks.push_back(trk); }
void StEvtInfo::AddBetaTrk(const StV0TrkInfo& trk){ m_VecBetaTrks.push_back(trk); }
