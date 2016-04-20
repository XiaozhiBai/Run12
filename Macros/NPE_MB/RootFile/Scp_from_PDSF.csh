#!/bin/tcsh

date

set BaseDir=/global/homes/x/xiao00/work/run12/submit
set Date=$1
set Dir=Root_File_${Date}  
mkdir ${Dir} 



scp -r xiao00@pdsf.nersc.gov:${BaseDir}/StRoot/StNpeMaker/StNpeMaker.cxx  ${Dir} 
scp -r xiao00@pdsf.nersc.gov:${BaseDir}/StRoot/StNpeMaker/StNpeMaker.h  ${Dir} 
scp -r xiao00@pdsf.nersc.gov:${BaseDir}/StRoot/StNpeMaker/StCuts.h  ${Dir} 

scp -r xiao00@pdsf.nersc.gov:${BaseDir}/production//hist_${Date}.root ${Dir}




#cd ../

