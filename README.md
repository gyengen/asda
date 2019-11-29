**This is an updated code of the ASDA code in https://github.com/PyDL/ASDA**

## Cite:
Liu, J., Nelson, C, Erdelyi, R, Automated Swirl Detection Algorithm (ASDA) and Its Application to Simulation and Observational Data, ApJ, 872, 22, 2019 (https://iopscience.iop.org/article/10.3847/1538-4357/aabd34/meta).

This version is a revision by Norbert Gyenge and Jiajia Liu of version 1.0 in https://github.com/PyDL/ASDA.
This repository is created in order to get ASDA ready for merging into sunkit-image.
This version is on average **3 times faster** than the non-MPI version of version 1.0.

# ASDA
Automatic Swirl Detection Algorithms.

## System Requirements
### OS Requirements
ASDA can be run on Windows, Mac OSX or Linux systems with Python 3 and the following dependencies installed.

### Dependencies:
**Python 3** with libraries including numpy, scipy, matplotlib, scikit-image, warnings, itertools.</br>
**Python management softwares including Anaconda or Virtualenv are recommended**

## Hardware Requirements:
ASDA requires a standard computer with enough CPU and computation power depending on the dataset used.

## Installation Guide:
ASDA is a stand-alone Python package, no installation is needed.

## Description of Files (More information can be found in each file):
**asda.py**: Main programm of the implication of the swirl detection algorithms</br>
**utils.py**: utility functions</br>
**demo.py**: demo</br>
**data**: this folder contains data need for running the demo</br>

## Instructions for Use:
Suppose you have vx and vy for the horizontal velocity field. These two arrays are in **(y, x)** order.</br>
1. import neccessary libraries, `import asda`
2. Initialize the ASDA class `lo = asda.ASDA_Calc(vx, vy)` </br>
3. calculate gamma1 and gamma2 values (see the reference) with `gamma = lo.gamma_values()`. The variable `gamma` will be a tuple, where `gamma[..., 0]` is gamma1 and `gamma[..., 1]` is gamma2.</br>
4. perform the detection of vortices using `center_edge = lo.center_edge()`. The variable `center_edge` will be a dictionary with keywords `"center", "edge", "points", and "peak"`. center is a list containing the pixel location of all swirls in the velocity field. edge is a list of the edges of all vortices. point is a list of all points within swirls. peak is a list of the peak gamma1 values of all swirls. radius is a list of the effective radii of all swirls. Note that center, points and edge are stored in the order of [x, y] and in units of pixels.</br>
5. use `ve, vr, vc, ia = lo.vortex_property(image=image)` to calculate the expanding, rotating, center speeds of above vortices. ia is the average intensity from the observational image (if given) of all points within each vortex. Note that vc is stored in the order of [vx, vy].</br>

## Demo
A demo **demo.py** is available:
1. To run the demo, run `python demo.py`
2. The demo data consists 2 files: **vxvy.npz** stores the velocity field (vx, vy) determined by FLCT using SOT Ca II observations at 2007-03-05T05:48:06.737 and at 2007-03-05T05:48:06.737, and the observational data (data) at 2007-03-05T05:48:06.737. **correct.npz** contains the swirl detection result from ASDA version 1.0 (https://github.com/PyDL/ASDA) on exactly the same dataset.

### Expected Running Time
On an Intel I7 4.20 GHz CPU: the first function in **demo.py** takes ~2 seconds, and the second function in **demo.py** takes ~ 27 seconds.
