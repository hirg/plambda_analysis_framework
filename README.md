Author: Liwen Wen
Email: wenliwen64@ucla.edu / wenliwen64@gmail.com
**************

## Contents
README.md                   --this file
Makefile                    --Makefile
include                     --consist of constants.h StV0TrkInfo.h StPriTrkInfo.h 
StEvtCuts                   --event-wise cuts
StEvtInfo                   --class used to store event-scope information
StPriTrkCuts                --abstract base class for specific primary track cuts classes
StPriETrkCuts               --electron track cuts class
StPriKTrkCuts               --kaon track cuts class
StPriPTrkCuts               --proton track cuts class
StPriTrkGeneralCuts         --general primary cuts class used for event plane reconstruction and so on
StRefMultCorr               --STAR's centrality definition class
StV0TrkCuts                 --v0 particle cuts class
compile.C                   --engine
*dat                        --efficiency data

## Manual

### Introduction
This repository includes the data analysis framework used for proton-lambda gamma correlation for auau200GeV run11 dataset.

### Configuration:
In `compile.C`, you can tune almost every cut applied in this analysis. You should specify the dataset's location.
