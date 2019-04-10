#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""SMILI M87 Stokes I Imaging Pipeline for EHT observations in April 2017

Authors: The Event Horizon Telescope Collaboration et al.
Date: April 10, 2019
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)
Data Product Code: 2019-D01-02

Brief Description:
The pipeline reconstructs an image from a specified uvfits file simultaneously
released in the EHT website (data release ID: 2019-D01-01) using a
python-interfaced imaging package SMILI (Akiyama et al. 2017a,b) version 0.0.0
(Akiyama et al. 2019).This script assumes to load calibrated data sets at Stokes
I (see M87 Paper III). As described in M87 Paper IV Section 6.2.3, the SMILI
pipeline also uses the eht-imaging library (Chael et al. 2016, 2018)
version 1.1.0 (Chael et al. 2019) solely for pre-imaging calibrations including
time averaging, LMT calibration to use the input visibility data sets consistent
with the eht-imaging imaging pipeline. The script has been tested in Python 2.7
installed with Anaconda on Ubuntu 18.04 LTS and MacOS 10.13 & 10.14 (with
macport or homebrew).

Notes:
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

The pipeline will output three files.
 - image (specified with -o option; assumed to be xxxx.fits here)
 - pre-calibrated uvfits file (xxxx.precal.uvfits)
 - self-calibrated uvfits file (xxxx.selfcal.uvfits)
For usage and detail parameters of the pipeline, please
read the help document associated in the imaging script
"python smili_imaging_pipeline.py --help".

References:
 - EHT Collaboration Data Portal Website:
   https://eventhorizontelescope.org/for-astronomers/data
 - The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)
 - The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)
 - The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)
 - The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)
 - The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)
 - The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)
 - Akiyama, K., Ikeda, S., Pleau, M., et al. 2017a, ApJ, 838, 1
 - Akiyama, K., Kuramochi, K., Ikeda, S., et al. 2017b, AJ, 153, 159
 - Akiyama, K., Tazaki, F., Moriyama, K., et al. 2019, Zenodo (SMILI version 0.0.0)
 - Chael, A., Johnson, M., Narayan, R., et al. 2016, ApJ, 829, 11C
 - Chael, A., Johnson, M., Bouman, K., et al. 2018, ApJ, 857, 23C
 - Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.1.0)
 - Anaconda: https://www.anaconda.com/
 - eht-imaging: https://github.com/achael/eht-imaging
 - SMILI: https://github.com/astrosmili/smili
"""


#-------------------------------------------------------------------------------
# Information of Authors
#-------------------------------------------------------------------------------
__author__ = "The Event Horizon Telescope Collaboration et al."
__copyright__ = "Copyright 2019, the Event Horizon Telescope Collaboration et al."
__license__ = "GPL version 3"
__version__ = "1.0"
__date__  = "April 10 2019"


#-------------------------------------------------------------------------------
# Modules
#-------------------------------------------------------------------------------
from smili import uvdata,imdata,imaging,util
import ehtim as eh
import pandas as pd
import numpy as np
import copy
import os
import argparse


#-------------------------------------------------------------------------------
# Loading command line arguments
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i','--input',metavar='input uvfits file',type=str,required=True,
                    help='input uvfits file name ')
parser.add_argument('-o','--output',metavar='output fits file',type=str,default=None,
                    help='output fits file name')
parser.add_argument('--day',metavar='observing day of April 2017',type=int,required=True,
                    help='Observational Day of April: 5, 6, 10 or 11. This must be specified correctly to apply the correct total flux density used in the network calibration.')
parser.add_argument("--prior",metavar='prior FWHM size',type=float,default=40.,
                    help='Prior FWHM (uas)')
parser.add_argument("--tfd",metavar='total flux density',type=float,default=0.5,
                    help='Target total flux density (Jy)')
parser.add_argument("--sys",metavar='fractional systematic error',type=float,default=0.,
                    help='Fractional systematic error to be added to visibilities in quadrature (percent)')
parser.add_argument("--keepsinglepol", action='store_true', default=False,
                    help='This option uses single pol station(s) for imaging; this will cause slight difference in the resultant images from M87 Paper IV.')
parser.add_argument("--lambl1",metavar='lambda l1',type=float,default=1.,
                    help='lambda l1')
parser.add_argument("--lambtv",metavar='lambda tv',type=float,default=1e3,
                    help='lambda tv')
parser.add_argument("--lambtsv",metavar='lambda tsv',type=float,default=1e4,
                    help='lambda tsv')
parser.add_argument("--nproc",metavar="Number of parallel processors. Default=1.",type=int,default=1)
args = parser.parse_args()

# Input uvfits file
inputuvfits = args.input

# Output directory
outputfits = args.output
if outputfits is None:
    outputfits = os.path.splitext(inputuvfits)[0]+".fits"
outputhead = os.path.splitext(outputfits)[0]

# Observational date
day = args.day
if day not in [5,6,10,11]:
    raise ValueError("The observational day (--day integer) must be (April) 5, 6, 10 or 11.")

# Prior size
prior_fwhm = args.prior

# Target total flux density
totalflux = args.tfd

# Systematic error
syserr = args.sys

# Lambda L1, TV and TSV
lambl1  = args.lambl1
lambtv  = args.lambtv
lambtsv = args.lambtsv

# Number of processors
nproc = args.nproc
os.environ['OMP_NUM_THREADS'] = '%d'%(nproc)

#-------------------------------------------------------------------------------
# Other tuning parameters fixed in the paper.
#-------------------------------------------------------------------------------
# The station(s) with single pol data
singlepol = ['JC']

# Number of iterative imaging inside a single selfcal cycle
Nweig = 3

# Number of closure imaging attempted in this script (must be < Nself)
Nflcl = 2

# Number of selfcals
Nself = 4

# Imaging pixel size
dx_uas = 2

# Number of pixels
nx = 64


#-------------------------------------------------------------------------------
# Set default imaging parameters
#-------------------------------------------------------------------------------
# Definition of the prior image
def gen_prior_gauss(x0=0.,y0=0.,fwhm=40.,dx=dx_uas,nx=nx):
    '''
    This function provides a Gaussain image with a specified FWHM in uas.
    '''
    # create a blank image
    outimage = imdata.IMFITS(dx=dx,nx=nx,angunit="uas")
    outimage = outimage.add_gauss(majsize=fwhm, overwrite=True)
    return outimage

# Set a circular Gaussian prior for L1
l1_prior = gen_prior_gauss(fwhm=prior_fwhm,dx=dx_uas,nx=nx)

# Default imaging parameters
imprm_init={}

# Iteration numbers
imprm_init["niter"] = 1000

# Regularization parameters
#   Flat prior TSV
if lambtsv > 0:
    imprm_init["tsv_lambda"] = lambtsv
else:
    imprm_init["tsv_lambda"] = -1

#   Flat prior TV
if lambtv > 0:
    imprm_init["tv_lambda"] = lambtv
else:
    imprm_init["tv_lambda"] = -1

#   Weighted L1
if lambl1 > 0:
    imprm_init["l1_lambda"] = lambl1
    imprm_init["l1_prior"] = l1_prior
else:
    imprm_init["l1_lambda"] = -1

#   Maximum Entropy Methods: not used in SMILI pipeline
#     GS-MEM
imprm_init["gs_lambda"] = -1
#     KL-MEM
imprm_init["kl_lambda"] = -1

# Total flux density constraint
imprm_init["totalflux"] = totalflux     # Target total flux density
imprm_init["tfd_lambda"] = 1            # Lambda parameter
imprm_init["tfd_tgterror"] = 1e-2       # Target fractional accuracy


#-------------------------------------------------------------------------------
# Step 1: Pre-calibration
#   To make sure that we use visibility data sets pre-calibrated consistently with
#   two other pipelines, this pipeline also uses eht-imaging library.
#-------------------------------------------------------------------------------
# Load the uvfits file
obs = eh.obsdata.load_uvfits(inputuvfits).switch_polrep('stokes')

# Rescale short baselines to excite contributions from extended flux.
#   Total flux density adopted in the Network Calibration
#   see EHT Collaboration et al. 2019c, M87 Paper III
if day < 7:
    netcal_tfd = 1.1  # 1.1 Jy was adopted for Apr  5 and  6
else:
    netcal_tfd = 1.2  # 1.2 Jy was adopted for Apr 10 and 11
#   Scaling the total flux density:
#     setting totalflux < netcal_tfd assumes there is an extended constant flux
#     component of netcal_tfd-totalflux Jy
for j in range(len(obs.data)):
    if (obs.data['u'][j]**2 + obs.data['v'][j]**2)**0.5 < 0.1e9:
        for k in range(-8,0):
            obs.data[j][k] *= totalflux/netcal_tfd

# Do scan averaging
obs.add_scans() # this seperates the data into scans, if it isn't done so already with an NX table
obs = obs.avg_coherent(0.,scan_avg=True) # average each scan coherantly

# Order stations. this is to create a minimal set of closure quantities with the highest snr
obs.reorder_tarr_snr()

# Self-calibrate the LMT to a Gaussian model if LMT is included in the input uvfits file
if "LM" in obs.data["t1"] or "LM" in obs.data["t2"]:
    obs_selfcal  = obs.copy()
    #   Add systematic noise for leakage and other non-closing errors
    #   (reminder: this must be done *after* any averaging)
    obs_selfcal = obs_selfcal.add_fractional_noise(syserr)
    obs_selfcal = obs_selfcal.switch_polrep('stokes')
    #   Set a circular Gaussian with a FWHM of 60 uas as the model image
    gausspriormodel = eh.image.make_square(obs, int(nx), nx*dx_uas*eh.RADPERUAS)           # Create empty image
    gausspriormodel = gausspriormodel.add_gauss(0.6, (60*eh.RADPERUAS, 60*eh.RADPERUAS, 0, 0, 0)) # Add gaussian
    #   Self-calibrate amplitudes
    for repeat in range(3):
        caltab = eh.self_cal.self_cal(
            obs_selfcal.flag_uvdist(uv_max=2e9),
            gausspriormodel,
            sites=['LM','LM'],
            method='vis',
            ttype='nfft',
            processes=nproc,
            caltable=True,
            gain_tol=1.0)
        obs_selfcal = caltab.applycal(obs_selfcal, interp='nearest', extrapolate=True)
    # Save calibrated uvfits data sets
    obs_selfcal.save_uvfits(outputhead+".precal.uvfits")
    del obs, obs_selfcal, caltab, gausspriormodel, netcal_tfd
else:
    obs.save_uvfits(outputhead+".precal.uvfits")
    del obs


#---------------------------------------------------------------------------
# Step 2: Generating Data Tables for Imaging
#---------------------------------------------------------------------------
# Create the full complex visibility table
vtable = uvdata.UVFITS(outputhead+".precal.uvfits").select_stokes("I").make_vistable()
if not args.keepsinglepol:
    for stname in singlepol:
        vtable = vtable.query("st1name != '"+stname+"' and st2name != '"+stname+"'").reset_index(drop=True)

# Check the number of data sets
if len(vtable.amp) == 0:
    raise ValueError("No data points exist in the input visibility data set.")

# Create tables for non-trivial/non-redundant bi-spectrum and closure amplitudes
btable = vtable.make_bstable(redundant=[["AA","AP"],["JC","SM"]])
ctable = vtable.make_catable(redundant=[["AA","AP"],["JC","SM"]])

# Check the number of data sets
if len(btable.amp) == 0:
    btable = None
if len(ctable.amp) == 0:
    ctable = None

#---------------------------------------------------------------------------
# Step 3: Imaging
#---------------------------------------------------------------------------
def edit_image(image,imprm):
    '''
    Image editing (shifting, blurring)
    '''
    image_edited = copy.deepcopy(image)
    orgtotalflux = image_edited.totalflux()

    # Shifting images such that the ring or dominant structure comes to the center
    if "vistable" not in imprm.keys():
        image_edited = image_edited.comshift()

    # Blurring with 20 uas Gaussian beam
    image_edited = image_edited.convolve_gauss(majsize=20,angunit="uas")

    # Scaling images
    if "totalflux" in imprm.keys():
        image_edited.data *= imprm["totalflux"]/image_edited.totalflux()
    else:
        image_edited.data *= orgtotalflux/image_edited.totalflux()

    image_edited.update_fits()
    return image_edited

# Loop for self-calibration
for iself in xrange(Nself):
    # Set the initial image
    #   At each iteration of self-calibrations, imaging starts from a circular
    #   Gaussian.
    initimage = imdata.IMFITS(
        dx=dx_uas,
        nx=nx,
        angunit="uas",
        uvfits=outputhead+".precal.uvfits",
    )
    initimage = initimage.add_gauss(totalflux=totalflux,majsize=prior_fwhm,overwrite=True)

    # Iterative Imaging:
    #   At each iteration of self-calibration, imaging iterations are attempted
    #   for Nweig(=5) times.
    for iweig in xrange(Nweig):
        # Imaging parameters: copied from the default parameter sets
        imprm = copy.deepcopy(imprm_init)

        # Data to be used in default
        #   Closure phases/amplitudes are used always to avoid over-fitting to
        #   full complex visibilities
        imprm["bstable"]=btable
        imprm["catable"]=ctable

        # Use amplitudes or full complex visibilities:
        #   Parameters depending on self-cal stages
        if iself < Nflcl:
            # The first half for self-calibration: using closure quantities

            # Get Amplitude Data sets
            atable = copy.deepcopy(vtable)

            # Flag the intra-site to avoid conflict with with the total-flux-density regularization.
            atable = atable.query("uvdist > 50e6")

            # Adding a fractional error to the data
            #  For LMT (adding 30% error)
            idx = atable["st1name"] == "LM"
            idx|= atable["st2name"] == "LM"
            atable.loc[idx,"sigma"] = np.sqrt(atable.loc[idx,"sigma"]**2+(atable.loc[idx,"amp"]*0.3)**2)
            #  For other stations (adding an error specified with amp_err)
            idx = idx == False
            atable.loc[idx,"sigma"] = np.sqrt(atable.loc[idx,"sigma"]**2+(atable.loc[idx,"amp"]*0.05)**2)

            # Use amplitudes
            imprm["amptable"]=atable
        else:
            # The latter half for self-calibration: using full complex
            # visibilities. Note that we still add closure quantities
            # to prevent over-fitting to full complex visibilities
            vtable_in = vtable.query("uvdist > 50e6")
            imprm["vistable"] = vtable_in

        # Imaging
        outimage = imaging.lbfgs.imaging(initimage,**imprm)
        initimage = edit_image(outimage,imprm)

    # Save image
    outimage.to_fits(outputfits)

    # If iself==Nself, self-calibration won't be performed and exit the script.
    if iself == Nself-1:
        break

    # Self-calibration
    # 1. Set the file name
    if iself == 0:
        uvfitsold = outputhead+".precal.uvfits"
    else:
        uvfitsold = outputhead+".selfcal.uvfits"
    uvfitsnew = outputhead+".selfcal.uvfits"

    # 2. Load uvfits
    uvfits = uvdata.UVFITS(uvfitsold)
    #
    # 3. Edit images
    #   Before selfcal, we bring the center of the mass to the image center
    #   and also normalize the image with the target total flux density.
    outimage = outimage.comshift()
    #   Normalize the total flux density
    #   *** Thanks to the total-flux-density regularization, this will provide only a
    #   *** slight change in the total flux density by a few pecent or even < 1%
    outimage.data[0,0] *= totalflux/outimage.totalflux()

    # 4. Do selfcal
    #   This is a regularized selfcal function using Gaussian priors
    #   for gain amplitudes and phases. std_amp is the standard deviation of
    #   the Gaussian prior for amplitude gains. std_amp=10000 works as almost
    #   the flat prior, so it is essentially selfcal without any gain constraints.
    caltable = uvfits.selfcal(outimage,std_amp=10000)
    uvfits = uvfits.apply_cltable(caltable)

    # 5. Save the self-calibrated uvfits file
    uvfits.to_uvfits(uvfitsnew)

    # 6. Get a new visibility table
    vtable = uvfits.select_stokes("I").make_vistable()
    if not args.keepsinglepol:
        for stname in singlepol:
            vtable = vtable.query("st1name != '"+stname+"' and st2name != '"+stname+"'").reset_index(drop=True)
