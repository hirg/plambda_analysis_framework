#ifndef STPRIETRKCUTS_H
#define STPRIETRKCUTS_H
#include <iostream>
#include "../StPriTrkCuts/StPriTrkCuts.h"
#include "../include/StPriTrkInfo.h"
class StPriETrkCuts: public StPriTrkCuts{
public:
    StPriETrkCuts(){}
    virtual ~StPriETrkCuts(){}
private:
    virtual bool PassPIDCuts(const StPriTrkInfo& p) const;

    virtual void DumpPIDCuts() const;
    ClassDef(StPriETrkCuts,1)
};
#endif
