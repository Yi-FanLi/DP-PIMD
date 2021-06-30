# Using DP-PIMD (DeepMD-kit + Path Integral Molecular Dynamics)
This document describes how you can use 

## 1. Install DeepMD-kit
### 1.1 Install bazel
To install tensorflow C++ API, we need to install bazel first. We can download bazel from github:

`wget https://github.com/bazelbuild/bazel/releases/download/3.7.2/bazel-3.7.2-dist.zip`

Unzip it:

`unzip bazel-3.7.2-dist.zip -d bazel-3.7.2`

Compile it:

```
cd bazel-3.7.2
./compile.sh
```

Add bazel to PATH environment variable:

`export PATH=/scratch/gpfs/yifanl/Softwares/bazel-3.7.2/output:$PATH`

### 1.2 Install tensorflow
The tensorflow C++ API version we are using is 2.5.0. We first clone it from github and swith into the `v2.5.0` branch:

```
git clone https://github.com/tensorflow/tensorflow.git
cd tensorflow
git checkout -b v2.5.0 v2.5.0
```

Configure it:

`./configure`

The configurations I am using are:

```
Please specify the location of python. [Default is /usr/bin/python3]: /scratch/gpfs/yifanl/Packages/tensorflow/venv/tf2.5-dp-0621v2/bin/python3.8


Found possible Python library paths:
  /scratch/gpfs/yifanl/Packages/tensorflow/venv/tf2.5-dp-0621v2/lib/python3.8/site-packages
Please input the desired Python library path to use.  Default is [/scratch/gpfs/yifanl/Packages/tensorflow/venv/tf2.5-dp-0621v2/lib/python3.8/site-packages]

Do you wish to build TensorFlow with ROCm support? [y/N]: N
No ROCm support will be enabled for TensorFlow.

Do you wish to build TensorFlow with CUDA support? [y/N]: y
CUDA support will be enabled for TensorFlow.

Do you wish to build TensorFlow with TensorRT support? [y/N]: N
No TensorRT support will be enabled for TensorFlow.

Asking for detailed CUDA configuration...

Please specify the CUDA SDK version you want to use. [Leave empty to default to CUDA 10]: 11.3


Please specify the cuDNN version you want to use. [Leave empty to default to cuDNN 7]: 8.2.0


Please specify the locally installed NCCL version you want to use. [Leave empty to use http://github.com/nvidia/nccl]: 


Please specify the comma-separated list of base paths to look for CUDA libraries and headers. [Leave empty to use the default]: /usr/local/cuda-11.3,/usr/local/cudnn/cuda-11.3/8.2.0/


Found CUDA 11.3 in:
    /usr/local/cuda-11.3/targets/x86_64-linux/lib
    /usr/local/cuda-11.3/targets/x86_64-linux/include
Found cuDNN 8 in:
    /usr/local/cudnn/cuda-11.3/8.2.0/lib64
    /usr/local/cudnn/cuda-11.3/8.2.0/include


Please specify a list of comma-separated CUDA compute capabilities you want to build with.
You can find the compute capability of your device at: https://developer.nvidia.com/cuda-gpus. Each capability can be specified as "x.y" or "compute_xy" to include both virtual and binary GPU code, or as "sm_xy" to only include the binary code.
Please note that each additional compute capability significantly increases your build time and binary size, and that TensorFlow only supports compute capabilities >= 3.5 [Default is: 3.5,7.0]: 6.0


Do you want to use clang as CUDA compiler? [y/N]: n
nvcc will be used as CUDA compiler.

Please specify which gcc should be used by nvcc as the host compiler. [Default is /usr/bin/gcc]: /scratch/gpfs/yifanl/Softwares/gcc-9.3.0/build/bin/gcc


Please specify optimization flags to use during compilation when bazel option "--config=opt" is specified [Default is -Wno-sign-compare]: 


Would you like to interactively configure ./WORKSPACE for Android builds? [y/N]: N
Not configuring the WORKSPACE for Android builds.

Preconfigured Bazel build configs. You can use any of the below by adding "--config=<>" to your build command. See .bazelrc for more details.
	--config=mkl         	# Build with MKL support.
	--config=mkl_aarch64 	# Build with oneDNN and Compute Library for the Arm Architecture (ACL).
	--config=monolithic  	# Config for mostly static monolithic build.
	--config=numa        	# Build with NUMA support.
	--config=dynamic_kernels	# (Experimental) Build kernels into separate shared objects.
	--config=v2          	# Build TensorFlow 2.x instead of 1.x.
Preconfigured Bazel build configs to DISABLE default on features:
	--config=noaws       	# Disable AWS S3 filesystem support.
	--config=nogcp       	# Disable GCP support.
	--config=nohdfs      	# Disable HDFS support.
	--config=nonccl      	# Disable NVIDIA NCCL support.
Configuration finished

```

## 2. Install LAMMPS