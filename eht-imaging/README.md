# eht-imaging M87 Stokes I Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** April 10, 2019

**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)

**Data Product Code:** [2019-D01-02](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

The pipeline reconstructs an image from uvfits files simultaneously released in the EHT website (data release ID: 2019-D01-01) using a python-interfaced imaging package eht-imaging (Chael et al. 2016,2018).

To run the pipeline, specify the input uvfits data file. Multiple bands can be specified separately. Additional flags control the output, which is only the reconstructed image as a FITS image by default.

Example call:

    python eht-imaging_pipeline.py -i ../../data/uvfits/SR1_M87_2017_101_lo_hops_netcal_StokesI.uvfits -i2 ../../data/uvfits/SR1_M87_2017_101_hi_hops_netcal_StokesI.uvfits --savepdf

For additional details, please read the help document associated in the imaging script: "python eht-imaging_pipeline.py --help".

**Additional References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)](https://doi.org/10.3847/2041-8213/ab0ec7)
- [The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)](https://doi.org/10.3847/2041-8213/ab0c96)
- [The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)](https://doi.org/10.3847/2041-8213/ab0c57)
- [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)
- [The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)](https://doi.org/10.3847/2041-8213/ab0f43)
- [The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)](https://doi.org/10.3847/2041-8213/ab1141)
- [Chael, A., Johnson, M., Narayan, R. et al. 2016, ApJ, 829, 11C](https://ui.adsabs.harvard.edu/abs/2016ApJ...829...11C/abstract)
- [Chael, A., Johnson, M., Bouman, K. et al. 2018, ApJ, 857, 23C](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...23C/abstract)
- [Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.1.0)](https://zenodo.org/record/2614016)
- [eht-imaging: https://github.com/achael/eht-imaging](https://github.com/achael/eht-imaging)
