#include "StPriTrkGeneralCuts.h"
#include "../include/StPriTrkInfo.h"
#include <iostream>
ClassImp(StPriTrkGeneralCuts)
bool StPriTrkGeneralCuts::PassPIDCuts(const StPriTrkInfo& p) const{ return true; }

void StPriTrkGeneralCuts::DumpPIDCuts() const{}
