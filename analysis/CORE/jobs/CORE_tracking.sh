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
root -l -b -q 'Fun4All_G4_CORE.C(10000)'

#Move output files and cleanup
echo "Cleaning Up..."
mv -v G4COREDetector_g4tracking_eval.root /sphenix/user/baschmoo/analysis/CORE/output/muon_single/G4COREDetector_g4tracking_eval_${INPUT}.root

echo "DONE!!!"
