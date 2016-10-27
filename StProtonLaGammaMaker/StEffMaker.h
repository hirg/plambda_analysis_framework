#ifndef STEFFMAKER_H
#define STEFFMAKER_H
#include <string>
using namespace std;
class StEffMaker{
    public:
        StEffMaker():m_Initialized(false){}
        virtual Bool_t Init() = 0;
        virtual Double_t GetEfficiency(Int_t cent, Double_t pt) const = 0;

        Bool_t IsInitialized() const{ return m_Initialized; }
        unsigned int GetCentBins() const{ return m_CentBins; }
        Double_t GetPtThreshold() const{ return m_PtThreshold; }
        string GetEffFileName() const{ return m_EffFileName; } 

        void SetInitializationStat(Bool_t status){ m_Initialized = status; }
        void SetCentBins(unsigned int centBins){ m_CentBins = centBins; }
        void SetPtThreshold(Double_t pt){ m_PtThreshold = pt; }
        void SetEffFileName(const string& name){ m_EffFileName = name; }

    private:
        Bool_t m_Initialized;
        unsigned int m_CentBins;
        Double_t m_PtThreshold;
        string m_EffFileName;
};
#endif
