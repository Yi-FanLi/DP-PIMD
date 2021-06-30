# Using DP-PIMD (DeepMD-kit + Path Integral Molecular Dynamics)

## 1. Install DeepMD-kit
### 1.1 Install bazel
To install tensorflow C++ API, we need to install bazel first. We can download bazel from github:
`wget https://github.com/bazelbuild/bazel/releases/download/3.7.2/bazel-3.7.2-dist.zip`

Unzip it:
`unzip bazel-3.7.2-dist.zip -d bazel-3.7.2`

Compile it:
`
    cd bazel-3.7.2
    ./compile.sh
`

Add bazel to PATH environment variable:
`export PATH=/scratch/gpfs/yifanl/Softwares/bazel-3.7.2/output:$PATH`

## 2. Install LAMMPS