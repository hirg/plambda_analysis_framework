#ifndef STPRITRKGENERALCUTS_H
#define STPRITRKGENERALCUTS_H
#include <iostream>
#include "../StPriTrkCuts/StPriTrkCuts.h"
#include "../include/StPriTrkInfo.h"
class StPriTrkGeneralCuts: public StPriTrkCuts{
public:
    StPriTrkGeneralCuts(){}
    virtual ~StPriTrkGeneralCuts(){}
private:
    virtual bool PassPIDCuts(const StPriTrkInfo& p) const;

    virtual void DumpPIDCuts() const;
    ClassDef(StPriTrkGeneralCuts,1)
};
#endif
