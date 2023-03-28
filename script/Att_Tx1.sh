#!/bin/bash
#source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-centos7-gcc62-opt/setup.sh

inputxsec=$1
inputmass=$2
inputmediator=$3
inputenergy=$4
outfile=$5

cd ../bin/
#echo $conddfig $pdf $toy $outfile
./ADM_MC_front $inputxsec $inputmass $inputmediator $inputenergy $outfile
