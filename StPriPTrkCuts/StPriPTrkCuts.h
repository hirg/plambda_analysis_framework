#ifndef STPRIPTRPCUTS_H
#define STPRIPTRPCUTS_H
#include <iostream>
#include "../StPriTrkCuts/StPriTrkCuts.h"
#include "../include/StPriTrkInfo.h"
#include "../include/constants.h"
class StPriPTrkCuts: public StPriTrkCuts{
public:
    StPriPTrkCuts(){}
    virtual ~StPriPTrkCuts(){}
private:
    virtual bool PassPIDCuts(const StPriTrkInfo& p) const;

    virtual void DumpPIDCuts() const;
    ClassDef(StPriPTrkCuts,1)
};
#endif
