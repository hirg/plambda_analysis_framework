#ifndef STEXPEFFMAKER_H
#define STEXPEFFMAKER_H
#include <string>
#include <fstream>
#include <vector>
#include "TF1.h"
#include "TMath.h"
#include "StEffMaker.h"
using namespace std;
class StTpcTofEffMaker: public StEffMaker{
    public: 
	StTpcTofEffMaker(){}

        virtual Bool_t Init(){
            string datfilename = GetEffFileName();
	    ifstream df(datfilename.c_str());
	    if(!df.good()){
		cout << "ERROR: BAD Efficiency File!" << endl;
		SetInitializationStat(false);
		return true;
	    }

	    Int_t cent;
	    Double_t p0;
	    Double_t p1;
	    Double_t p2;

	    while(df >> cent){
		df >> p0 >> p1 >> p2;
		TF1 feff("ftrking", "[0] * exp(-pow([1] / x, [2]))", 0, 3.0);
		feff.SetParameter(0, p0); 
		feff.SetParameter(1, p1); 
		feff.SetParameter(2, p2); 
		m_VecEffFunc_Tpc.push_back(feff);
                if(cent == 9)
                    break;
	    }

            while(df >> cent){
                df >> p0;//>> p1 >> p2;
                //TF1 feff("feff", "[0] * TMath::Erf([1] * x - [2])", 0, 3.0);
                //TF1 feff("ftof", "[0] + [1] * x + [2] * x * x", 0, 3.0);
                TF1 feff("ftof", "[0]", 0, 3.0);
                feff.SetParameter(0, p0);
                //feff.SetParameter(1, p1);
                //feff.SetParameter(2, p2);
                m_VecEffFunc_Tof.push_back(feff); 
	    }

	    df.close();
            SetInitializationStat(true);
	    return true;
	}

        virtual Double_t GetEfficiency(Int_t cent, Double_t pt) const{
             return IsInitialized()?  (m_VecEffFunc_Tpc[cent - 1].Eval(pt) * m_VecEffFunc_Tof[cent - 1].Eval(pt)):1.0;
	}

    private: 
        vector<TF1> m_VecEffFunc_Tpc;
        vector<TF1> m_VecEffFunc_Tof;
};
#endif
