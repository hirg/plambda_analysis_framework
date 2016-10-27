#ifndef STCUTS_H
#define STCUTS_H
#include "StEventInfo.h"
class StCuts{
    public:
        StCuts(){}
        virtual ~StCuts(){}
        virtual Bool_t CheckCuts(const StEventInfo& evtInfo) = 0;
    private:
}
#endif
