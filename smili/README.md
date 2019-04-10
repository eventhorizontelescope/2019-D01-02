# SMILI M87 Stokes I Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** April 10, 2019

**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)

**Data Product Code:** [2019-D01-02](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

The pipeline reconstructs an image from a specified uvfits file simultaneously
released in the EHT website (data release ID: [2019-D01-01](https://eventhorizontelescope.org/for-astronomers/data)) using a python-interfaced imaging package SMILI (Akiyama et al. 2017a,b) version 0.0.0 (Akiyama et al. 2019). This script assumes to load calibrated data sets at Stokes I (see M87 Paper III). As described in M87 Paper IV Section 6.2.3, the SMILI pipeline also uses the eht-imaging library (Chael et al. 2016, 2018) version 1.1.0 (Chael et al. 2019) solely for pre-imaging calibrations including time averaging, LMT calibration to use the input visibility data sets consistent with the eht-imaging imaging pipeline. The script has been tested in Python 2.7 installed with Anaconda on Ubuntu 18.04 LTS and MacOS 10.13 & 10.14 (with macport or homebrew).

The pipeline will output three files.

- image (specified with -o option; assumed to be xxxx.fits here)
- pre-calibrated uvfits file (xxxx.precal.uvfits)
- self-calibrated uvfits file (xxxx.selfcal.uvfits)

For usage and detail parameters of the pipeline, please read the help document associated in the imaging script "*python smili_imaging_pipeline.py --help*".

We also include an example driver *example_driver.py* to run *smili_imaging_pipeline.py* on all released data sets. You can run it by "*python example_driver.py --uvfitsdir < xxxx/uvfits > --nproc < number of the processors >*" For details, please have a look at the help document by "*python example_driver.py --help*".

**Notes**

We note that, as described in 2019-D01-01, released visibility data sets are
slightly different from data sets used in Paper IV for two reasons.

(1) They have only Stokes I, while Paper IV data sets have dual polarization at
Stokes RR and LL. This slight difference in released data sets will change
self-calibration procedures slightly in this pipeline since R and L gains are
calibrated separately for the latter dual polarization data sets. We find that
this does not produce any significant differences (within < ~5%) in resultant images.

(2) In Paper IV, this pipeline strictly uses Stoke I data; JCMT which has a
single polarization is flagged when Stokes I visibilities are computed from RR
and LL, while JCMT is included for self-calibration at the corresponding
polarization. On the other hand, in the released data sets, JCMT has pseudo
Stokes RR/LL data, such that Stokes I computed from them is identical to the
original single Stokes. This will make JCMT be used in imaging, which change
images slightly. For reproducibility of Paper IV results using this data set, in
default, we will not use JCMT for imaging, as we did in Paper IV. If you want
to use it to use all of data sets, you can specify --keepsinglepol to keep it.

**References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)](https://doi.org/10.3847/2041-8213/ab0ec7)
- [The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)](https://doi.org/10.3847/2041-8213/ab0c96)
- [The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)](https://doi.org/10.3847/2041-8213/ab0c57)
- [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)
- [The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)](https://doi.org/10.3847/2041-8213/ab0f43)
- [The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)](https://doi.org/10.3847/2041-8213/ab1141)
- [Akiyama, K., Ikeda, S., Pleau, M. et al. 2017a, ApJ, 838, 1](https://ui.adsabs.harvard.edu/#abs/2017ApJ...838....1A)
- [Akiyama, K., Kuramochi, K., Ikeda, S. et al. 2017b, AJ, 153, 159](https://ui.adsabs.harvard.edu/#abs/2017AJ....153..159A)
- [Akiyama, K., Tazaki, F., Moriyama, K. et al. 2019, Zenodo (SMILI version 0.0.0)](https://zenodo.org/record/2616725)
- [Chael, A., Johnson, M., Narayan, R. et al. 2016, ApJ, 829, 11C](https://ui.adsabs.harvard.edu/abs/2016ApJ...829...11C/abstract)
- [Chael, A., Johnson, M., Bouman, K. et al. 2018, ApJ, 857, 23C](https://ui.adsabs.harvard.edu/abs/2018ApJ...857...23C/abstract)
- [Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.1.0)](https://zenodo.org/record/2614016)
- [Anaconda: https://www.anaconda.com/](https://www.anaconda.com/)
- [eht-imaging: https://github.com/achael/eht-imaging](https://github.com/achael/eht-imaging)
- [SMILI: https://github.com/astrosmili/smili](https://github.com/astrosmili/smili)
