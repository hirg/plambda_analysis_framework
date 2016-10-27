#ifndef STPRIPITRPCUTS_H
#define STPRIPITRPCUTS_H
#include <iostream>
#include "../StPriTrkCuts/StPriTrkCuts.h"
#include "../include/StPriTrkInfo.h"
class StPriPiTrkCuts: public StPriTrkCuts{
public:
    StPriPiTrkCuts(){}
    virtual ~StPriPiTrkCuts(){}
private:
    virtual bool PassPIDCuts(const StPriTrkInfo& p) const;

    virtual void DumpPIDCuts() const;
    ClassDef(StPriPiTrkCuts,1)
};
#endif
