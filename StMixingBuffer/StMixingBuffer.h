#ifndef STMIXINGBUFFER_H
#define STMIXINGBUFFER_H
#include <vector>
#include "TRandom3.h"
#include "../StEvtInfo/StEvtInfo.h"
class StMixingBuffer{
public:
    StMixingBuffer(int nCells, int cellSize): m_VecEvt(nCells),  m_CellSize(cellSize), m_RndSeed(0){ cout << "happy" << endl; }
    virtual ~StMixingBuffer(){}
    virtual void AddEntry(const StEvtInfo& evtInfo, int cellIdx);

    virtual const vector<StEvtInfo>& GetCell(int cellIdx) const;
    virtual void DumpInfo() const;
    virtual size_t size() const { return m_VecEvt.size(); }
private:
    vector< vector<StEvtInfo> >  m_VecEvt;
    int m_CellSize;
    int m_RndSeed;

    ClassDef(StMixingBuffer, 1);
};
#endif
