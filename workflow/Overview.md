# Using DP-PIMD (DeePMD-kit + Path Integral Molecular Dynamics)
This document describes how we can use DeepMD-kit and LAMMPS to run path integral molecular dynamics(PIMD). 
## 1. Install DeePMD-kit and LAMMPS
The easiest way to do install DeePMD-kit is via off-line packages from
```https://github.com/deepmodeling/deepmd-kit/releases/download/v2.0.0.b4/```

However, the development of LAMMPS for PIMD has not been contributed into the LAMMPS official repository, so it has not been included into the off-line packages. Therefore, we need to build DeePMD-kit and LAMMPS from source.

The system we use is:
```
Linux version 3.10.0-1160.25.1.el7.x86_64
```
To make all things clear, we set up a new `anaconda` environment and then install the compilers and cudatoolkit in it. We refer to the installation procedure in this `gist` to set up the `conda` environment, install DeePMD-kit and LAMMPS:
```
https://gist.github.com/y1xiaoc/5cc980252db59e484c39060179959c67
```
Please note that our installation differs from that of the `gist` in the following ways:
- We install our own modified lammps so that the `fix dp_pimd` style is included. Please clone the LAMMPS from Yifan's github instead of the official site:
```
git clone https://github.com/Yi-FanLi/lammps.git
```
After that, check out to `beads` branch:
```
git checkout beads
```
- The `REPLICA` pakcage of LAMMPS should be installed.
- In addition to the packages specified in the gist, we also need to enable `nvcc` via `conda install -c conda-forge cudatoolkit-dev`.
- In this project we do not need PLUMED, so installing PLUMED is not necessary. 



## 2. Train the forcefield
In this project we use the forcefield in the paper *Phase Diagram of a Deep Potential Water Model* (https://doi.org/10.1103/PhysRevLett.126.236001). The training data and the trained model can be downloaded from DP Library: 
```
http://dplibrary.deepmd.net/#/project_details?project_id=202010.001
````
DP Library is free to use, but we need to sign up an account and then download the data and model. 

After downloading the data, we can train the model from scratch. The input script has been included in the tarball. Here, we just convert the model trained with DeePMD-kit v1.1 to a v2.0-compatible model. The command we use is `dp convert-from`.
```
dp convert-from 1.2 -i frozen_model.pb
```
It will give us a file `convert_out.pb` which can be used by DeePMD-kit v2.0.

## 3. Compress the forcefield
To speed up the MD simulations, we can compress the model so that the force calculation can be much faster. We need the input script used when training the forcefield to compress the model, but the `input.json` given in the tarball is for DeePMD-kit v1.1, so we need to convert it.

First unarchive the tarball:
```
tar -xjvf 2020.09.final.clean.tar.bz2
cd 2020.09.final.clean/data
```
Then start the training process for several steps:
```
mkdir train
cd train
cp ../../train_scripts/input.json .
vi input.json
40G
dd
452G
dd
:wq
dp train input.json
```
After several seconds it will generate an `out.json` file. Then we use `Ctrl+C` to shut down the training and compress the model:
```
cd ..
mkdir compress
cd compress
cp ../train/out.json .
cp ../../../convert_out.pb .
dp compress out.json -i convert_out.pb
```
After this run finishes, we will obtain a `frozen_model_compressed.pb` file. This is the model we will use.

## 4. Run the PIMD task
We use the newly-developed fix style `dp_pimd` in LAMMPS to run this task. The `lmp_run` folder contains 8 `data` files, which are the LAMMPS initial configuration files of the 8 beads. This folder also include an `in.lammps` file, which tells lammps which commands to execute. Also, it includes a `run.slurm` script which submits tasks to the Princeton tiger cluster. You can change the configurations in it to run this task on your own cluster. The last line of `run.slurm` runs the LAMMPS task:
```
mpirun -np 8 lmp -in in.lammps -p 8x1 -log log -screen screen
```
Note that the `-p MxN` option of LAMMPS specifies the partition configuration of the LAMMPS run. `M` specifies how many replicas LAMMPS uses, and is the number of beads in PIMD. `N` denotes how many processes are used for each replica. Moreover, the number of processes specified by `-np` in the mpi task must be equal to `MxN`.

The followings are the commands I use to submit the task:
```
cd ../
cp 2020.09.final.clean/data/compress/frozen_model_compressed.pb lmp_run/
cd lmp_run
sbatch run.slrum
```

