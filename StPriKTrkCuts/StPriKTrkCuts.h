#ifndef STPRIKTRKCUTS_H
#define STPRIKTRKCUTS_H
#include <iostream>
#include "../StPriTrkCuts/StPriTrkCuts.h"
#include "../include/StPriTrkInfo.h"
class StPriKTrkCuts: public StPriTrkCuts{
public:
    StPriKTrkCuts(){}
    virtual ~StPriKTrkCuts(){}
private:
    virtual bool PassPIDCuts(const StPriTrkInfo& p) const;

    virtual void DumpPIDCuts() const;
    ClassDef(StPriKTrkCuts,1)
};
#endif
