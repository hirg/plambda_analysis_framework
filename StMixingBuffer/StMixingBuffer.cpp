#include "StMixingBuffer.h"
#include "../StEvtInfo/StEvtInfo.h"
#include <iostream>
#include "TRandom3.h"
#include <vector>
ClassImp(StMixingBuffer);
void StMixingBuffer::AddEntry(const StEvtInfo& evtInfo, int cellIdx){
    vector<StEvtInfo>& cell = m_VecEvt[cellIdx];
    TRandom3* gRandom = new TRandom3(m_RndSeed++);

    int idx = -1;
    if(cell.size() == m_CellSize){
	do{
	    idx = (int)(m_CellSize * (1.0 - gRandom->Rndm()));
	}while(idx < 0 || idx >= m_CellSize);
	cell.erase(cell.begin() + idx);
    }

    cell.push_back(evtInfo);

    delete gRandom;
}

const vector<StEvtInfo>& StMixingBuffer::GetCell(int cellIdx) const{
    return m_VecEvt[cellIdx];
}

void StMixingBuffer::DumpInfo() const{
    cout << "Total number of cells: " << m_VecEvt.size() << endl;
    cout << "Total number of mixing evts: " << m_VecEvt[0].size() << endl;
}
