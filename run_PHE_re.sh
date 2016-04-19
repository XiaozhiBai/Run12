#!/bin/bash

#root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Gamma_HT.list","Embedding/Gamma_HT.root")' > Embedding/LogGamma_HT &
#root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Pi0_dalitz_HT.list","Embedding/Pi0_dalitz_HT.root")' > Embedding/LogPi0_HT &
#root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Eta_dalitz_HT.list","Embedding/Eta_dalitz_HT.root")' > Embedding/LogEta_HT &

root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Gamma_MB.list","Embedding/Gamma_MB.root")' > Embedding/LogGamma_MB &
