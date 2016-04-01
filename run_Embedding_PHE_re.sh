#!/bin/bash

root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Gamma_HT.list","Embedding/Gamma_HT.root")'
root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Pi0_dalitz.list","Embedding/Pi0_dalitz_HT.root")'
root4star -l -b -q StRoot/StEmbeddingAna/macors/run_StPhotonicReEff.C'("Embedding/Eta_dalitz.list","Embedding/Eta_dalitz_HT.root")'

