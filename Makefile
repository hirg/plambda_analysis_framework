run:
	@for i in 1 2 3 4 5 6 7 8 9 ; do root -l -b -q compile.C\($$i\) ; done
dist:
	tar -cvzf plambda_auau200GeVrun11.tar.gz include StEvtCuts StEvtInfo StMixingBuffer StPriETrkCuts StPriKTrkCuts StPriPiTrkCuts StPriPTrkCuts StPriTrkCuts StPriTrkGeneralCuts StProtonLaGammaMaker StRefMultCorr StV0TrkCuts compile.C *dat Makefile

clean:
	rm -rf *tar.gz
