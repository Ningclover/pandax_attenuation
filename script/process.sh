#!/bin/bash
#source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-centos7-gcc62-opt/setup.sh

mcfile=$1
outfile=$2
root -l -b -q creat_ture_sum_f.C\(\"${mcfile}\",\"${outfile}\"\)

mcfile2=$3
outfile2=$4
root -l -b -q creat_ture_sum_b.C\(\"${mcfile2}\",\"${outfile2}\"\)

outfile_f=$5
flux_txt_f=$6
#root -l -b -q creat_ture_sum_4plot_f.C\(\"${mcfile}\",\"${outfile_f}\",\"${flux_txt_f}\"\)

outfile_b=$7
flux_txt_b=$8
#root -l -b -q creat_ture_sum_4plot_b.C\(\"${mcfile2}\",\"${outfile_b}\",\"${flux_txt_b}\"\)

sigma=$9
mx=${10}
ms=${11}
file_txt=${12}

root -l -b -q Trans_Tr_form.C\($sigma,$mx,$ms,\"${outfile}\",\"${outfile2}\",\"${file_txt}\"\)

file_spec_root=${13}
root -l -b -q spec.C\(\"${file_txt}\",\"${file_spec_root}\"\)



