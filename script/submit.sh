#!/bin/sh

countjobs() {
  condor_q -format "1\n" "" $1 | wc -l
}

countjobs0() {
        condor_q |grep / | wc -l
}

username=ningxuyang

listfile=../dat/test.lst
cat $listfile | while read line
do
declare -a tmp=($line)
inputenergy=../dat/dm_flux/flux_mass_5e-03.txt
inputxsec=`printf "%.5e" ${tmp[1]}`
inputmass=`printf "%.3f" ${tmp[0]}`
inputmediator=`printf "%.2f" ${tmp[2]}`
outdir=../output/Att_outfile/Mass${inputmass}MeV/
logdir=../log/

echo $inputxsec $inputmass $inputmediator
if [ ! -d "${logdir}" ];then
mkdir -p $logdir
else
echo "${logdir} already exist" 
fi
if [ ! -d "${outdir}" ];then
mkdir -p $outdir
else
echo "${outdir} already exist" 
fi


for nn in {1..2}
do
        outfile=${outdir}DM${inputmass}mev_${inputxsec}cm2_ms${inputmediator}_1GeV_${nn}.root
        #echo $inputxsec $inputmass $inputmediator $inputenergy $outfile
        while [ $(countjobs ${username}) -ge 300 ]
        do
        sleep 10;
        done
        condor_submit<<EOF
        universe = vanilla
        executable = Att_Tx1.sh
        arguments = $inputxsec $inputmass $inputmediator $inputenergy $outfile
        output = ${logdir}/fmass${inputmass}mev_${inputxsec}cm2.output
        error = ${logdir}/fmass${inputmass}mev_${inputxsec}cm2.err
        log = ${logdir}/fmass${inputmass}mev_${inputxsec}cm2.log
	queue
EOF
done
done
