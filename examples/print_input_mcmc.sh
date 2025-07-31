#!/bin/bash -x


## read in data
inputfile=paramsinput_contlimitL16short.csv
#~ inputfile=paramsinput_contlimitL3.csv
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
startheat=($(tail -n +2 $inputfile | cut -d ',' -f15) )

len=${#Nts[@]}
echo "$len"

pathtoanalysisscripts=$(readlink -f ..)

echo " " >> commandsmeasmcmc.txt
echo " " >> commandsRmcmc.txt
for (( i=0; i<($len); i++))
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

xifile=$(echo "scale=2; xifile=${xis[$i]}; if (xifile<1) print 0 ; xifile" | bc -l | awk '{printf "%.2f", $0}')
betafile=$(echo "scale=2; xifile=${betas[$i]}; if (xifile<1) print 0 ; xifile" | bc -l | awk '{printf "%.4f", $0}')
betasix=$(echo "scale=6; xifile=${betas[$i]}; if (xifile<1) print 0 ; xifile" | bc -l | awk '{printf "%.6f", $0}')
jobname=b${betas[$i]}Ns${Nss[$i]}Nt${Nts[$i]}
resdirfile=/hiskp4/gross/masterthesis/analyse/results/
confdirfile=/hiskp4/gross/masterthesis/analyse/configs/confL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile
inputfile=inputs/mcmcinputL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile
inputfileqbig=/hiskp4/gross/masterthesis/analyse/inputs/mcmcinputL${Nss[$i]}.T${Nts[$i]}.b$betafile.x$xifile

# echo $confdirfile

beta=${betas[$i]}
betas[$i]=$(echo "$beta"| awk '{printf "%f", $0}')
xi=${xis[$i]}
xis[$i]=$(echo "$xi"| awk '{printf "%f", $0}')
#~ echo $xi
delta=$(echo "delta=1/${betas[$i]} ; if ($xi < 0.5) delta*=$xi ; if (delta<1) print 0 ; delta" | bc -l | awk '{printf "%f", $0}')
nmeasmc=$(echo "scale=0; nmeas=${nmeas[$i]}+1 ; nmeas" | bc -l | awk '{printf "%d", $0}')
nmeasme=$(echo "scale=0; nmeas=${nmeas[$i]}/${nsave[$i]} ; nmeas" | bc -l | awk '{printf "%d", $0}')
skip=${skip[$i]}  
 

## input file
printf "geometry:\n  X: ${Nss[$i]}\n  Y: ${Nss[$i]}\n  Z: 1\n  T: ${Nts[$i]}\n  ndims: 3\n\n" > $inputfile 
printf "monomials:\n  gauge:\n    beta: ${betas[$i]}\n    anisotropic:\n      xi: $xi\n\n" >> $inputfile
printf "metropolis:\n  do_mcmc: true\n  conf_dir: $confdirfile\n  conf_basename: config_u1\n  n_meas: $nmeasmc\n  N_save: ${nsave[$i]}\n  delta: $delta\n  restart_condition: $heat\n\n" >> $inputfile
printf "  lenghty_conf_name: true\n\n" >> $inputfile
#~ printf "  seed: 9874321\n\n" >> $inputfile
printf "omeas:\n  offline:\n    conf_dir: $confdirfile\n    lenghty_conf_name: true\n\n" >> $inputfile
printf "  res_dir: $resdirfile\n  icounter: 0\n  n_meas: $nmeas\n  nstep: ${nsave[$i]}\n  " >> $inputfile
printf "potential:\n    potentialplanar: ${potentialplanar}\n    potentialnonplanar: ${potentialnonplanar}\n    sizeWloops: ${fraction[$i]}\n" >> $inputfile
printf "    n_apesmear: ${napes[$i]}\n    alpha: ${alphas[$i]}\n    append: true" >> $inputfile

echo "sbatch --job-name=$jobname --cpus-per-task=${cpus[$i]} scriptyamlheatbath.sh $inputfileqbig $confdirfile" >> commandsmeasmcmc.txt


## write commands for R analysis scripts
if [ "${Nss[$i]}" -eq "3" ]; then
echo "Rscript ${pathtoanalysisscripts}/L3singleplaquette.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --skip $skip --plotpath . --shapiro" >> commandsRmcmc.txt
else

## compute effective masses, with and without drawing bootstrapsamples
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6 --drawbootstrap" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R    --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6 --drawbootstrap" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R    --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --analyse --plotuwerr --uwerrs 6" >> commandsRmcmc.txt

## analyse potential, determine xi and r0, ...
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1" >> commandsRmcmc.txt
           
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysissubtracted.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRmcmc.txt
          
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1" >> commandsRmcmc.txt
           
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 0 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 0 --errortotpot" >> commandsRmcmc.txt
echo "Rscript ${pathtoanalysisscripts}/analysisrotated.R --myfunctions ${pathtoanalysisscripts} --respath $resdirfile -b ${betas[$i]} -r ${Nss[$i]} -t ${Nts[$i]} --xidiff --xi $xi --betaone ${betaone[$i]} --skip $skip --nsave ${nsave[$i]} --every ${every[$i]} --zerooffset ${offset[$i]} --bootl 1 --aic --scaletauint --dofit --plotuwerr --uwerrs 6 --omit 1 --lowlimpot 0 --lowlim 1 --errortotpot" >> commandsRmcmc.txt

echo "sleep 1" >> commandsRmcmc.txt
fi

done
