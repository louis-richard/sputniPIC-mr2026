
# sputniPIC
sputniPIC is a space plasma simulation software that uses the Particle-In-Cell (PIC) method.

## Get the code
Git clone the repository

```bash
$ git clone https://github.com/louis-richard/sputniPIC-mr2026.git sputniPIC-mr2026
cd sputniPIC-mr2026
```

## Building
### On Linux

Compile with make

```bash
$ make
$ make -j 4 # build with 4 threads
$ make -j # build with maximum threads
```

Alternatively you can use CMake to build sputniPIC. Create the build directory

```bash
$ mkdir build && cd build
```

Generate the make files

```bash
$ cmake ..
```

Compile with make

```bash
$ make
$ make -j 4 # build with 4 threads
$ make -j # build with maximum threads
```

### On Mac OS

Install homebrew 

```bash
$ /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Install sputniPIC dependencies

```bash
$ brew install hdf5-mpi
```

Compile with make

```bash
$ make
$ make -j 4 # build with 4 threads
$ make -j # build with maximum threads
```


## Usage

Create a run directory with subdirectories for the output and the restart files
 
```bash
$ cd /homelocal/username/sputniPIC-mr2026/
$ mkdir ./run00
$ mkdir ./run00/data
$ mkdir ./run00/restart
```

Copy the input files and change the path to the output data

```bash
$ scp ./inputfiles/GEM_2D.inp ../run00/GEM_2D.inp
$ vi ./run00/GEM_2D.inp

SaveDirName = /homelocal/username/sputniPIC-mr2026/run00/data
RestartDirName = /homelocal/username/sputniPIC-mr2026/run00/restart

```

## Run

```bash
$ ./bin/sputniPIC_CPU.out ./run00/GEM_2D.inp 


-------------------------
sputniPIC Sim. Parameters
-------------------------
Number of species    = 4
Number of particles of species 0 = 4096000	 (MAX = 4096000)  QOM = -64
Number of particles of species 1 = 4096000	 (MAX = 4096000)  QOM = 1
Number of particles of species 2 = 4096000	 (MAX = 4096000)  QOM = -64
Number of particles of species 3 = 4096000	 (MAX = 4096000)  QOM = 1
x-Length                 = 40
y-Length                 = 20
z-Length                 = 1
Number of cells (x)      = 256
Number of cells (y)      = 128
Number of cells (z)      = 1
Time step                = 0.25
Number of cycles         = 20
Results saved in: /homelocal/louisr/simulations/sputniPIC-mr2026/run00/data
Output directory /homelocal/louisr/simulations/sputniPIC-mr2026/run00/restart exists.
*************************************************
**  Initialize GEM Challenge with Pertubation  **
*************************************************
** B0x = 0.0195
** B0y = 0
** B0z = 0
** Delta (current sheet thickness) = 0.5
** rho species 0 = 1 CURRENT SHEET
** rho species 1 = 1 CURRENT SHEET
** rho species 2 = 0.1 BACKGROUND
** rho species 3 = 0.1 BACKGROUND
*************************************************
Writing ic_data...
GPU has free memory: 50740461568x2
CUDA return check: Disabled
Total number of MPI ranks: 1
Number of cores per rank: 8
Number of GPUs per rank: 2
Threads Per Block of GPUs: 256
Total number of particles: 16384000
Number of particles per rank: 16384000; 375 MB of data
Allocating 23 MB of memory for particles on gpu
batch_size per species of 256000 (5 MB)

***********************
   cycle = 1
***********************
***  MOVER  ITERATIONS = 3 - Species 0 *** on gpu 0 - 16 batches
***  MOVER  ITERATIONS = 3 - Species 1 *** on gpu 1 - 16 batches
***  MOVER  ITERATIONS = 3 - Species 2 *** on gpu 0 - 16 batches
***  MOVER  ITERATIONS = 3 - Species 3 *** on gpu 1 - 16 batches
***********************
*** DIVERGENCE CLEANING ***
Initial error: 0.120194
CG converged at iteration # 32
*** MAXWELL SOLVER ***
Initial residual: 0.174301 norm b vector (source) = 0.274216
GMRES converged at restart # 0; iteration #15 with error: 0.00078713
*** B CALCULATION ***
Timing Cycle 1 : 4e-08 0.122447 0 0.25759 3.31e-07

***********************
   cycle = 2
***********************
***  MOVER  ITERATIONS = 3 - Species 0 *** on gpu 0 - 16 batches
***  MOVER  ITERATIONS = 3 - Species 1 *** on gpu 1 - 16 batches
***  MOVER  ITERATIONS = 3 - Species 2 *** on gpu 0 - 16 batches
***  MOVER  ITERATIONS = 3 - Species 3 *** on gpu 1 - 16 batches
***********************
*** DIVERGENCE CLEANING ***
Initial error: 0.108081



```

## Output format
- Every Iterations
```Timing Cycle <Iterations> : <Mover Time> <Interp Time> <Field Time> <IO Time>```
- Finally
```
Mover: <average> <standard deviation>
Field: <average> <standard deviation>
IO: <average> <standard deviation>
```

# Postprocessing
## On Odin

Create virtual environment and activate

```bash
$ python -m venv .venv
$ source .venv/bin/activate
```

Install the required dependencies

```bash
$ pip install numpy xarray matplotlib vtk tqdm
```

Run the example gem_2d.py

```bash
$ python3 postprocessing/gem_2d.py ./run00/data/ --time 10
```



## License
The software is released under BSD 2-Clause license. See LICENSE for details.

## Cite us
If you find our implementation and [paper](https://arxiv.org/pdf/2008.04397.pdf) useful to your work, we would apprecite if you cite us with:
```
@INPROCEEDINGS{9235052,
  author={S. W. D. {Chien} and J. {Nylund} and G. {Bengtsson} and I. B. {Peng} and A. {Podobas} and S. {Markidis}},
  booktitle={2020 IEEE 32nd International Symposium on Computer Architecture and High Performance Computing (SBAC-PAD)}, 
  title={sputniPIC: An Implicit Particle-in-Cell Code for Multi-GPU Systems}, 
  year={2020},
  volume={},
  number={},
  pages={149-156},
  doi={10.1109/SBAC-PAD49847.2020.00030}}
```
