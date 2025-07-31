#!/bin/bash -x

## read in data
inputfile=paramsinput_heatbath.csv
betas=($(tail -n +2 $inputfile | cut -d ',' -f1) )
Nss=($(tail -n +2 $inputfile | cut -d ',' -f2) )
Nts=($(tail -n +2 $inputfile | cut -d ',' -f3) )
xis=($(tail -n +2 $inputfile | cut -d ',' -f4) )
napes=($(tail -n +2 $inputfile | cut -d ',' -f5) )
alphas=($(tail -n +2 $inputfile | cut -d ',' -f6) )
cpus=($(tail -n +2 $inputfile | cut -d ',' -f7) )
betaone=($(tail -n +2 $inputfile | cut -d ',' -f8) )
fraction=($(tail -n +2 $inputfile | cut -d ',' -f9) )
skip=($(tail -n +2 $inputfile | cut -d ',' -f10) )
nmeas=($(tail -n +2 $inputfile | cut -d ',' -f11) )
nsave=($(tail -n +2 $inputfile | cut -d ',' -f12) )
offset=($(tail -n +2 $inputfile | cut -d ',' -f13) )
every=($(tail -n +2 $inputfile | cut -d ',' -f14) )
noverrelax=($(tail -n +2 $inputfile | cut -d ',' -f15) )
nheatbath=($(tail -n +2 $inputfile | cut -d ',' -f16) )
startheat=($(tail -n +2 $inputfile | cut -d ',' -f17) )

len=${#Nts[@]}
echo "$len"


pathtoanalysisscripts=$(readlink -f ..)

echo " " >> commandsmeasheatbath.txt
echo " " >> commandsRheatbath.txt
for (( i=0; i<${len}; i++))
do

##set heat
heat="cold"
heatbool="false"
if (( $(echo "${startheat[$i]}" | bc -l) ))
then 
heat="hot"
heatbool="true"
fi

## measure potentials for small or large volume
potentialplanar="true"
potentialnonplanar="false"
if [ "${Nss[$i]}" -eq "3" ]; then
potentialplanar="false"
potentialnonplanar="true"
fi


## print input file for generating configs
xifile=$(echo "scale=2; xifile=${xis[$i]}; if (xifile<1) print 0 ; xifile" | bc -l | awk '{printf "%.2f", $0}')
betafile=$(echo "scale=2; xifile=${betas[$i]}; if (xifile<1) print 0 ; xifile" | bc -l | awk '{printf "%.4f", $0}')
betasix=$(echo "scale=6; xifile=${betas[$i]}; if (xifile<1) print 0 ; xifile" | bc -l | awk '{printf "%.6f", $0}')
jobname=b${betas[$i]}Ns${Nss[$i]}Nt${Nts[$i]}nover${noverrelax[$i]}nheat${nheatbath[$i]}_$heat
## where to store config
confdirfile=/hiskp4/gross/masterthesis/analyse/configs/heatbath/confL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile.nover${noverrelax[$i]}.nheat${nheatbath[$i]}_$heat/
#~ confdirfile=/lustre/scratch/data/cgross2_hpc-matching/confs/confL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile.nover${noverrelax[$i]}.nheat${nheatbath[$i]}_$heat/
## where to store results
resdirfile=$confdirfile
## name of input file
inputfile=inputs/heatbathinputL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile.nover${noverrelax[$i]}.nheat${nheatbath[$i]}_$heat
## location of input file
inputfileqbig=/hiskp4/gross/masterthesis/analyse/inputs/heatbath/heatbathinputL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile.nover${noverrelax[$i]}.nheat${nheatbath[$i]}_$heat
#~ inputfileqbig=/lustre/scratch/data/cgross2_hpc-matching/input/heatbathinputL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile.nover${noverrelax[$i]}.nheat${nheatbath[$i]}_$heat

echo $confdirfile

beta=${betas[$i]}
betas[$i]=$(echo "$beta"| awk '{printf "%f", $0}')
xi=${xis[$i]}
xis[$i]=$(echo "$xi"| awk '{printf "%f", $0}')
nmeas=$(echo "scale=0; nmeas=${nmeas[$i]}+1 ; nmeas" | bc -l | awk '{printf "%d", $0}')
skip=${skip[$i]} 

## input file
printf "geometry:\n  X: ${Nss[$i]}\n  Y: ${Nss[$i]}\n  Z: 1\n  T: ${Nts[$i]}\n  ndims: 3\n\n" > $inputfile 
printf "monomials:\n  gauge:\n    beta: ${betas[$i]}\n    anisotropic:\n      xi: $xi\n\n" >> $inputfile
printf "heatbath_overrelaxation:\n  do_mcmc: true\n  n_meas: $nmeas\n  N_save: ${nsave[$i]}\n  restart_condition: $heat\n" >> $inputfile
# use this line if you want to append to already existing data
#~ printf "heatbath_overrelaxation:\n  do_mcmc: true\n  n_meas: $nmeas\n  N_save: ${nsave[$i]}\n  restart_condition: read\n" >> $inputfile
printf "  conf_dir: $confdirfile\n  n_overrelax: ${noverrelax[$i]}\n  n_heatbath: ${nheatbath[$i]}\n" >> $inputfile
printf "  lenghty_conf_name: true\n\n" >> $inputfile
#~ printf "  seed: 9874321\n\n" >> $inputfile
printf "omeas:\n  offline:\n    conf_dir: $confdirfile\n    lenghty_conf_name: true\n\n" >> $inputfile
printf "  res_dir: $resdirfile\n  icounter: 0\n  n_meas: $nmeas\n  nstep: ${nsave[$i]}\n  " >> $inputfile
printf "potential:\n    potentialplanar: ${potentialplanar}\n    potentialnonplanar: ${potentialnonplanar}\n    sizeWloops: ${fraction[$i]}\n" >> $inputfile
printf "    n_apesmear: ${napes[$i]}\n    alpha: ${alphas[$i]}\n    append: true" >> $inputfile


echo "sbatch --job-name=$jobname --cpus-per-task=${cpus[$i]} scriptyamlheatbath.sh $inputfileqbig" >> commandsmeasheatbath.txt

## write commands for R analysis scripts


if [ "${Nss[$i]}" -eq "3" ]; then
echo "Rscript ${pathtoanalysisscripts}/L3singleplaquette.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --skip $skip --plotpath . --shapiro"  >> commandsRheatbath.txt
else
## compute effective masses, with and without drawing bootstrapsamples
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6 --drawbootstrap" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R    --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6 --drawbootstrap" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R    --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6" >> commandsRheatbath.txt

## analyse potential, determine xi and r0, ...
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1" >> commandsRheatbath.txt
           
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRheatbath.txt
          
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1" >> commandsRheatbath.txt
           
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRheatbath.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRheatbath.txt

echo "sleep 1" >> commandsRheatbath.txt
fi

done

