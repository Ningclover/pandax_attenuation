#!/bin/bash 
countjobs(){
    condor_q -format "1\n" "" $1 | wc -l
}

username=ningxuyang
# listfile
listfile=../dat/test.lst
dirtot=../output/
logdir=../log/

cat $listfile | while read line
do
declare -a tmp=($line)
mass=`printf "%.3f" ${tmp[0]}`
xsec=`printf "%.5e" ${tmp[1]}`
ms=`printf "%.2f" ${tmp[2]}`
#echo ${mass} ${xsec} $ms
br=`printf "%.5e" ${tmp[3]}`
gu=`printf "%.5e" ${tmp[4]}`
mcfile1=${dirtot}/Att_outfile/Mass${mass}MeV/DM${mass}mev_${xsec}cm2_ms${ms}_1GeV_
outfile1=${dirtot}Att_spec/Tx_files/Tx_${mass}mev_${xsec}cm2_ms${ms}_1GeV.root
mcfile2=${dirtot}/Att_outfile_back/Mass${mass}MeV/DM${mass}mev_${xsec}cm2_ms${ms}_1GeV_
outfile2=${dirtot}Att_spec/Tx_files_back/Tx_${mass}mev_${xsec}cm2_ms${ms}_1GeV.root

outfile_f=${dirtot}Att_spec/Tx_files/Tx_${mass}mev_${xsec}cm2_ms${ms}_1GeV_4p.root
flux_txt_f=${dirtot}Att_spec/Tx_files/flux_${xsec}_${mass}mev_ms${ms}_1GeV.txt
outfile_b=${dirtot}Att_spec/Tx_files_back/Tx_${mass}mev_${xsec}cm2_ms${ms}_1GeV_4p.root
flux_txt_b=${dirtot}Att_spec/Tx_files_back/flux_${xsec}_${mass}mev_ms${ms}_1GeV.txt

file_txt=${dirtot}Att_spec/cs_${xsec}_${mass}mev_ms${ms}_1GeV.txt
spec_root=${dirtot}Att_spec/ADM_Rate_${mass}MeV_${xsec}_ms${ms}_p4_1GeV.root

if [ ! -d "${dirtot}Att_spec/Tx_files" ];then
mkdir -p ${dirtot}Att_spec/Tx_files
fi
if [ ! -d "${dirtot}Att_spec/Tx_files_back" ];then
mkdir -p ${dirtot}Att_spec/Tx_files_back
fi


#echo ${mcfile1} ${outfile1} ${mcfile2} ${outfile2} ${outfile_f} ${flux_txt_f} ${outfile_b} ${flux_txt_b} $xsec $mass $ms $file_txt $spec_root
      while [ $(countjobs ${username}) -ge 300 ]
      do
        sleep 5;
      done

    condor_submit<<EOF
    universe = vanilla
    executable = process.sh
    arguments =  ${mcfile1} ${outfile1} ${mcfile2} ${outfile2} ${outfile_f} ${flux_txt_f} ${outfile_b} ${flux_txt_b} $xsec $mass $ms $file_txt $spec_root
    output = ${logdir}${mass}${xsec}.output
    error  = ${logdir}${mass}.err
    log =    ${logdir}${mass}.log
    queue
EOF
done

