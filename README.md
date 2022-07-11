# fun4all_eicmacros &ndash; CORE detector
<br/>

Getting the CORE repository
---------------------------
Do the following:
```
git clone https://github.com/bschmookler/fun4all_eicmacros.git
cd fun4all_eicmacros
git checkout CORE
```
Then check that you are on the CORE branch:
```
git branch -a
```
or
```
git status
```

N.B. You will need to clone via SSH if you plan to push changes to the GitHub repository.

Setting up the environment
--------------------------
On the SDCC (BNL) machines, do the following (assuming C shell):
```
source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/eic_setup.csh -n
```

On JLab (ifarm), do the following (again, assuming C shell):
```
setenv LANG C
source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin/eic_setup.csh -n
```

You can also work locally if you have Singularity and CVMFS installed. See [here](https://github.com/ECCE-EIC/Singularity).

Running single-particle simulation
----------------------------------
To run a simulation for 100 events, do the following:
```
root -l -b -q 'Fun4All_G4_CORE.C(100)'
```

Analysis examples
-----------------
Several example analysis and job submission codes can be found in this [directory](analysis/CORE).
