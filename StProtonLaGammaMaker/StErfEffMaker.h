#ifndef STERFEFFMAKER_H
#define STERFEFFMAKER_H
#include <string>
#include <fstream>
#include <vector>
#include "TF1.h"
#include "StEffMaker.h"
using namespace std;
class StErfEffMaker: public StEffMaker{
    public: 
	StErfEffMaker(){}

        virtual Bool_t Init(){
            string datfilename = GetEffFileName();
	    ifstream df(datfilename.c_str());
	    if(!df.good()){
		cout << "ERROR: BAD Erf Efficiency File!" << endl;
		SetInitializationStat(false);
		return true;
	    }

	    Int_t cent;
	    Double_t p0;
	    Double_t p1;
	    Double_t p2;

	    while(df >> cent){
		df >> p0 >> p1 >> p2;
		TF1 feff("feff", "[0] * TMath::Erf([1] * x - [2])", 0.5, 2.6);
		feff.SetParameter(0, p0); 
		feff.SetParameter(1, p1); 
		feff.SetParameter(2, p2); 
		m_VecEffFunc_Exp.push_back(feff);
		if(m_VecEffFunc_Exp.size() == GetCentBins())
		    break;
	    }

	    while(df >> cent){
		df >> p0 >> p1 >> p2;
		TF1 feff("feff", "[0] * TMath::Erf([1] * x - [2])", 0.5, 2.6);
		feff.SetParameter(0, p0); 
		feff.SetParameter(1, p1); 
		feff.SetParameter(2, p2); 
		m_VecEffFunc_Flat.push_back(feff);
	    }

	    df.close();
            SetInitializationStat(true);
	    return true;
	}

        virtual Double_t GetEfficiency(Int_t cent, Double_t pt) const{
             return IsInitialized()? ((pt >= m_PtThreshold)? m_VecEffFunc_Flat[cent - 1].Eval(pt):m_VecEffFunc_Exp[cent - 1].Eval(pt)) : 1; 
	}

    private: 
        vector<TF1> m_VecEffFunc_Flat;
        vector<TF1> m_VecEffFunc_Exp;
};
#endif
