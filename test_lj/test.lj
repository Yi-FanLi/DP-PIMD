# Test LAMMPS PIMD on the Lennard-Jones Potential
## 1. Install LAMMPS
The configuration toold and compilers I use are
| Tool      | Version |
| :-----------: | :-----------: |
| cmake      | 3.17.5       |
| gcc   | 9.3.1 20200408 (Red Hat 9.3.1-2)        |
| openmpi | 3.1.5 |
The installation process of LAMMPS is
```
git clone https://github.com/Yi-FanLi/lammps.git
cd lammps
mkdir build
cd build
cmake ../cmake/ -DCMAKE_INSTALL_PREFIX=`pwd` -DPKG_REPLICA=on
make -j32 install
```
Then an executable `lmp` should be generated in this folder.

## 2. Run PIMD simulation of the Lennard-Jones 
The test input is in the `DP-PIMD/rho=1_nvt` folder of this repo. The `in.lj_nvt` specifies the LAMMPS task. Use the `run.slurm` script to submit the LAMMPS task. Edit your `run.slurm` to make sure to set the `LAMMPS_DIR` as your real path to the `lmp` executable and add it to the `PATH` environment variable. Submit the job using
```
cd DP-PIMD/rho=1_nvt
sbatch run.slurm
```