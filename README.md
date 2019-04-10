# First M87 EHT Results: Imaging Pipelines

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** April 10, 2019

**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)

**Data Product Code:** [2019-D01-02](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**

We release three imaging pipelines (DIFMAP, eht-imaging and SMILI)
used in the parameter survey of M87 Paper IV Section 6 and later. All
imaging pipelines create images from calibrated uvfits files (see M87
Paper III) simultaneously released ([data product code:
2019-D01-01](https://eventhorizontelescope.org/for-astronomers/data)). For
more detailed instructions, please see the README file in the
sub-directory for each pipeline.

We note that, as described in
[2019-D01-01](https://eventhorizontelescope.org/for-astronomers/data),
released visibility data sets have only Stokes *I*, which are slightly
different from data sets used in Paper IV that have dual polarization
at Stokes `RR` and `LL`. This slight difference in released data sets
will provide no net changes in DIFMAP and eht-imaging pipelines, while
it will change self-calibration procedures slightly for SMILI
calibrating `R` and `L` gains separately for the latter dual
polarization data sets. We confirm that reconstructed images are
consistent with images presented in Paper IV on fiducial parameters
(see Paper IV Section 6) for all three pipelines, and will not affect
our conclusions in the M87 publications (Paper I, II, III, IV, V and
VI).

**Notes:**

These data files only include Stokes *I* visibilities, while the
published results used data files from the full SR1 release which
include Stokes `RR` and `LL`. The slight difference in the underlying
data from the conversion to Stokes *I* in the single-precision
`*.uvfits` files, as well as differences in the python dependencies
used by `eht-imaging` and `SMILI` (e.g. `numpy`, `astropy`, `scipy`),
may slightly affect the final image.

**References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)](https://doi.org/10.3847/2041-8213/ab0ec7)
- [The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)](https://doi.org/10.3847/2041-8213/ab0c96)
- [The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)](https://doi.org/10.3847/2041-8213/ab0c57)
- [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)
- [The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)](https://doi.org/10.3847/2041-8213/ab0f43)
- [The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)](https://doi.org/10.3847/2041-8213/ab1141)
