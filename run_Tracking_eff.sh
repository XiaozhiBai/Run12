#!/bin/bash

starver SL12d;
root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StTrackingEfficiency.C'("Embedding/Single_electron_HT.list","Embedding/Tracking_eff_HT.root")' > Embedding/Log_tracking_efficiency_HT & 
root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StTrackingEfficiency.C'("Embedding/Single_electron_MB.list","Embedding/Tracking_eff_MB.root")' > Embedding/Log_tracking_efficiency_MB & 


