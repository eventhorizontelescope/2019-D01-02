# DIFMAP M87 Stokes I Imaging Pipeline for EHT observations in April 2017

**Authors:** The Event Horizon Telescope Collaboration et al.

**Date:** April 10, 2019

**Primary Reference:** [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)

**Data Product Code:** [2019-D01-02](https://eventhorizontelescope.org/for-astronomers/data)

**Brief Description:**
DIFMAP script for imaging EHT data using a given mask (set of cleaning windows)

Requires Caltech's DIFMAP software for CLEAN imaging reconstruction, version v2.5b or above, that can be obtained from [ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html](ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html)

To reproduce fiducial EHT DIFMAP images of M87:

1. Place i) uv FITS file "your_file_name.uvfits", ii) Difmap script EHT_Difmap (this file), and iii) mask CircMask_r30_x-0.002_y0.022.win in the same directory

2. Run Difmap and use the calling sequence:

    @EHT_Difmap your_file_name,CircMask_r30_x-0.002_y0.022,-10,0.5,0.1,2,-1

The header of the EHT_Difmap file contains further information regarding the choice of script parameters

**References:**

- [EHT Collaboration Data Portal Website](https://eventhorizontelescope.org/for-astronomers/data)
- [The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)](https://doi.org/10.3847/2041-8213/ab0ec7)
- [The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)](https://doi.org/10.3847/2041-8213/ab0c96)
- [The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)](https://doi.org/10.3847/2041-8213/ab0c57)
- [The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)](https://doi.org/10.3847/2041-8213/ab0e85)
- [The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)](https://doi.org/10.3847/2041-8213/ab0f43)
- [The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)](https://doi.org/10.3847/2041-8213/ab1141)
- [Shepherd, M., Pearson, T., & Taylor, G. B., 1994, BAAS, 26, 987S](https://ui.adsabs.harvard.edu/abs/1994BAAS...26..987S/abstract)
- [Shepherd, M., Pearson, T., & Taylor, G. B., 1995, BAAS, 27, 903S](https://ui.adsabs.harvard.edu/abs/1995BAAS...27..903S/abstract)
- [DIFMAP: ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html](ftp://ftp.astro.caltech.edu/pub/difmap/difmap.html)
