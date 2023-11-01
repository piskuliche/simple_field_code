# Electric Field Calculation

This program is used to calculate the electric field around the OHs of water molecules for use with computational spectroscopy codes. This code is generalizable to a variety of systems, and it outputs the information needed for calculating IR, Raman, and other higher-order spectroscopies using the empirical mapping approach.

## Installation

### Prequisites:

Installation requires the following:

1) cmake > 3.13
2) HDF5

### How to Install

Installation of the code is easy using cmake, just use the following commands to install the command.

```
mkdir build
cd build
cmake ../
make
```

### Other Programming Notes

In the Support_Codes/Fortran directory, there is a demonstration for how a fortran program can read the HDF5 files that are written by the program. This is located in `read_field.f90`. A more thorough primer on HDF5 can be found https://docs.hdfgroup.org/hdf5/develop/_getting_started.html 

## How-To Guide

Calculations with this code are relatively simple. 

### Required Files

What you need:
1) An XYZ trajectory file, ```traj.xyz```
    - right now, I plan to expand this to XTC files written by GROMACS; however, for now a python code ```grab_frames.py``` is provided in the Support_Codes/Python which converts XTC files to XYZ.
2) An inputfile ```field_input.in```
    - An example of this input file is located in the Examples directory.
3) An input file for each molecule
    - Two examples of this, ```popc.in``` and ```water.in``` are included in the Examples directory.
4) ```L.dat```
    - This is a file with the xyz dimensions of the box, it is written by the same ```grab_frames.py``` file that can be found in the Support_Codes/Python directory (when converting from XTC) or must be created separately. One line per configuration.
5) OPTIONAL: ```samples.in``` a file that includes the information about what molecules should be calculated.
    - This is a simple file, 1 line is 1 molecule ID number, numbered from 1. 
    - These molecules can be contiguous (e.g., 1, 2, 3) or randomized (e.g. 5, 10, 1)


### Running the Codes

To run the code, make sure all of these files are in the same directory, and then run the following command:

```
./get_field_calc
```

The result will be saved as a file ```field.h5``` om the same directory. This file is in the HDF5 file format, which is a binary archive format to store the data from these calculations efficiently.

The hierarchy of this file is as follows:

field.h5

- dot_ioh

- eoh_ioh

- z0_ioh

Here, ioh stands for the number of OHs, so for 200 waters, there would be 400 OHs. For instance the 101st OH is referenced as ```dot_101```.

### Sampling

Now these calculations have a BIG problem, which I hope to eventually solve by adding MPI support, which is that as you get to big #'s of molecules, the calculation gets incredibly expensive. This also begins to pose memory management issues as well. The way the code solves this is by using sampling, where it does a full field calculations - only for a subset of water molecules. This option is enabled by setting samples to a number, and having a ```samples.in``` file that has at least that many lines in it. Likewise, this option can be turned off by setting nsamples in ```field_input.in``` to 0. 

## Copyright, Zeke Piskulich, 2023
