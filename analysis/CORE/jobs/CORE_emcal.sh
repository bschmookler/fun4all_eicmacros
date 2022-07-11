#!/usr/bin/bash
#-----#-----#----#----#----#----#----#----#----#

#Go into scratch directory
chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}

#Make subdirectory and move there
INPUT=$(( 0 + $1 ))
echo $INPUT
DIR=`printf "%04d" $INPUT`
mkdir $DIR
cd $DIR

#Copy CORE directory. Find a better way to do this in future
cp -r /sphenix/user/baschmoo/myfork/fun4all_eicmacros/detectors/CORE .
cd CORE

#Run simulation
echo "start running in directory $PWD"
echo "Running Job Number $1"

momentum=(0.5 1.0 5.0 10.0 20.0)

for i in "${momentum[@]}"
do
	#Electrons
	root -l -b -q 'Fun4All_G4_CORE_EEMC.C(5000,'$i',"e-")'
	mv -v G4COREDetector_g4tracking_eval.root /sphenix/user/baschmoo/analysis/CORE/output/electron_single_emcal/G4COREDetector_g4tracking_eval_${i}_${INPUT}.root
	mv -v G4COREDetector_g4eemc_eval.root /sphenix/user/baschmoo/analysis/CORE/output/electron_single_emcal/G4COREDetector_g4eemc_eval_${i}_${INPUT}.root

	#Pions
 	root -l -b -q 'Fun4All_G4_CORE_EEMC.C(5000,'$i',"pi-")'
        mv -v G4COREDetector_g4tracking_eval.root /sphenix/user/baschmoo/analysis/CORE/output/pion_single_emcal/G4COREDetector_g4tracking_eval_${i}_${INPUT}.root
        mv -v G4COREDetector_g4eemc_eval.root /sphenix/user/baschmoo/analysis/CORE/output/pion_single_emcal/G4COREDetector_g4eemc_eval_${i}_${INPUT}.root

done

echo "DONE!!!"
