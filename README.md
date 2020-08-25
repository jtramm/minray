# minray

## Overview

A minimal implementation of the random ray method for neutral particle transport. This purpose of this code is to serve as an open source performance benchmark for random ray, while still being able to perform the real neutronics and produce real solutions.

Key Features:
  - Much simpler than the full application (ARRC) -- allowing it to act as a reference implementation
  - Complex enough to converge real solutions to neutron transport problems (see picture below generated with minray)
  - Simplified geometry treatment limits complexity (can only simulate 2D heterogeneous Cartesian geometries)
  - Hard coded reactor simulation problem of interest: 2D C5G7
  - Kernelized for easy porting to accelerator languages

![minray](docs/img/2D_C5G7_thermal_flux_wm.jpg)

## Installing, Compiling, and Running

Minray is written in C with no dependencies.

Use the included makefile to install. Several options (optimization, OpenMP usage, debugging) are available as toggles at the top of the makefile.

Run the appliation as `./minray` to get the default problem. Other options are:

Usage: ./minray <options>
Options:
    -r <rays>                    Number of discrete rays
    -d <distance per ray>        Travel distance per ray (cm)
    -i <inactive iterations>     Set fixed number of inactive power iterations
    -a <active iterations>       Set fixed number of active power iterations
    -s <seed>                    Random number generator seed (for reproducibility)
    -m <problem size multiplier> Multiplioer to increase/decrease problem size/resolution
    -p                           Enables plotting

## Citing minray

Papers wishing to cite random ray, ARRC, or minray should in general refer to:

John R. Tramm, Kord S. Smith, Benoit Forget, and Andrew R. Siegel.  The Random
Ray Method for neutral particle transport. Journal of Computational Physics,
342:229 - 252, 2017. https://doi.org/10.1016/j.jcp.2017.04.038

For convenience, the bibtex entry is given below:

@article{Tramm2017,
title = {{The Random Ray Method} for neutral particle transport},
journal = "Journal of Computational Physics",
volume = "342",
pages = "229 - 252",
year = "2017",
issn = "0021-9991",
doi = "https://doi.org/10.1016/j.jcp.2017.04.038",
url = "http://www.sciencedirect.com/science/article/pii/S0021999117303170",
author = "John R. Tramm and Kord S. Smith and Benoit Forget and Andrew R. Siegel",
}
