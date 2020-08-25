# minray

## Overview

A minimal implementation of the random ray method for neutral particle transport. This purpose of this code is to serve as an open source performance benchmark for random ray, while still being able to perform the real neutronics and produce real solutions.

Key Features:
  - Much simpler than the full application (ARRC) -- allowing it to act as a reference implementation
  - Complex enough to converge real solutions to neutron transport problems (see picture below generated with minray)
  - Simplified geometry treatment limits complexity (can only simulate 2D heterogeneous Cartesian geometries)
  - Hard coded reactor simulation problem of interest: 2D C5G7
  - Kernelized for easy porting to accelerator languages

![minray](docs/img/2D_C5G7_thermal_flux_wm.png)

## Citing minray

Papers citing the ARRC program should in general refer to:

John R. Tramm, Kord S. Smith, Benoit Forget, and Andrew R. Siegel.  The Random
Ray Method for neutral particle transport. Journal of Computational Physics,
342:229 - 252, 2017. https://doi.org/10.1016/j.jcp.2017.04.038

John R. Tramm, Kord S. Smith, Benoit Forget, and Andrew R. Siegel. ARRC: A
random ray neutron transport code for nuclear reactor simulation. Annals of
Nuclear Energy, 112 (Supplement C): 693 – 714, 2018.
https://doi.org/10.1016/j.anucene.2017.10.015

For convenience, the bibtext entries are given below:

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

@article{Tramm2018,
title = {{ARRC}: A random ray neutron transport code for nuclear reactor simulation},
journal = "Annals of Nuclear Energy",
volume = "112",
number = "Supplement C",
pages = "693 - 714",
year = "2018",
issn = "0306-4549",
doi = "https://doi.org/10.1016/j.anucene.2017.10.015",
url = "http://www.sciencedirect.com/science/article/pii/S0306454917303444",
author = "John R. Tramm and Kord S. Smith and Benoit Forget and Andrew R. Siegel",
}
