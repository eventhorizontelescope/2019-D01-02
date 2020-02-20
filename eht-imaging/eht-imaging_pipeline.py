"""eht-imaging M87 Stokes I Imaging Pipeline for EHT observations in April 2017

Authors: The Event Horizon Telescope Collaboration et al.
Date: April 10, 2019
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)
Data Product Code: 2019-D01-02

Brief Description:
The pipeline reconstructs an image from uvfits files simultaneously
released in the EHT website (data release ID: 2019-D01-01) using a
python-interfaced imaging package eht-imaging (Chael et al. 2016,2018).

To run the pipeline, specify the input uvfits data file. Multiple bands
can be specified separately. Additional flags control the output, which
is only the reconstructed image as a FITS image by default.

Example call:
python eht-imaging_pipeline.py -i ../SR1_M87_2017_101_lo_hops_netcal_StokesI.uvfits -i2 ../SR1_M87_2017_101_hi_hops_netcal_StokesI.uvfits --savepdf

For additional details, please read the help document associated in the
imaging script: "python eht-imaging_pipeline.py --help".

Additional References:
 - EHT Collaboration Data Portal Website:
   https://eventhorizontelescope.org/for-astronomers/data
 - The Event Horizon Telescope Collaboration, et al. 2019a, ApJL, 875, L1 (M87 Paper I)
 - The Event Horizon Telescope Collaboration, et al. 2019b, ApJL, 875, L2 (M87 Paper II)
 - The Event Horizon Telescope Collaboration, et al. 2019c, ApJL, 875, L3 (M87 Paper III)
 - The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)
 - The Event Horizon Telescope Collaboration, et al. 2019e, ApJL, 875, L5 (M87 Paper V)
 - The Event Horizon Telescope Collaboration, et al. 2019f, ApJL, 875, L6 (M87 Paper VI)
 - Chael, A., Johnson, M., Narayan, R., et al. 2016, ApJ, 829, 11C
 - Chael, A., Johnson, M., Bouman, K., et al. 2018, ApJ, 857, 23C
 - Chael, A., Bouman, K., Johnson, M. et al., Zenodo (eht-imaging version 1.1.0)
 - eht-imaging: https://github.com/achael/eht-imaging
"""

#-------------------------------------------------------------------------------
# Author Information
#-------------------------------------------------------------------------------
__author__ = "The Event Horizon Telescope Collaboration et al."
__copyright__ = "Copyright 2019, the Event Horizon Telescope Collaboration et al."
__license__ = "GPL version 3"
__version__ = "1.0"
__date__  = "April 10 2019"

#-------------------------------------------------------------------------------
# Modules
#-------------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')

import os
import argparse
import ehtim as eh
import numpy as np

#-------------------------------------------------------------------------------
# Load command-line arguments
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Fiducial eht-imaging script for M87")
parser.add_argument('-i', '--infile',    default="obs.uvfits",help="input UVFITS file")
parser.add_argument('-i2', '--infile2',  default="",          help="optional 2nd input file (different band) for imaging")
parser.add_argument('-o', '--outfile',   default='out.fits',  help='output FITS image')
parser.add_argument('--savepdf',         default=False,       help='saves image pdf (True or False)?',action='store_true')
parser.add_argument('--imgsum',          default=False,       help='generate image summary pdf',action='store_true')
args = parser.parse_args()

#-------------------------------------------------------------------------------
# Fiducial imaging parameters obtained from the eht-imaging parameter survey
#-------------------------------------------------------------------------------
zbl        = 0.60               # Total compact flux density (Jy)
prior_fwhm = 40.0*eh.RADPERUAS  # Gaussian prior FWHM (radians)
sys_noise  = 0.02               # fractional systematic noise
                                # added to complex visibilities

# constant regularization weights
reg_term  = {'simple' : 100,    # Maximum-Entropy
             'tv'     : 1.0,    # Total Variation
             'tv2'    : 1.0,    # Total Squared Variation
             'l1'     : 0.0,    # L1 sparsity prior
             'flux'   : 1e4}    # compact flux constraint

# initial data weights - these are updated throughout the imaging pipeline
data_term = {'amp'    : 0.2,    # visibility amplitudes
             'cphase' : 1.0,    # closure phases
             'logcamp': 1.0}    # log closure amplitudes


#-------------------------------------------------------------------------------
# Fixed imaging parameters
#-------------------------------------------------------------------------------
obsfile   = args.infile         # Pre-processed observation file
ttype     = 'nfft'              # Type of Fourier transform ('direct', 'nfft', or 'fast')
npix      = 64                  # Number of pixels across the reconstructed image
fov       = 128*eh.RADPERUAS    # Field of view of the reconstructed image
maxit     = 100                 # Maximum number of convergence iterations for imaging
stop      = 1e-4                # Imager stopping criterion
gain_tol  = [0.02,0.2]          # Asymmetric gain tolerance for self-cal; we expect larger values
                                # for unaccounted sensitivity loss
                                # than for unaccounted sensitivity improvement
uv_zblcut = 0.1e9               # uv-distance that separates the inter-site "zero"-baselines
                                # from intra-site baselines
reverse_taper_uas = 5.0         # Finest resolution of reconstructed features

# Specify the SEFD error budget
# (reported in First M87 Event Horizon Telescope Results III: Data Processing and Calibration)
SEFD_error_budget = {'AA':0.10,
                     'AP':0.11,
                     'AZ':0.07,
                     'LM':0.22,
                     'PV':0.10,
                     'SM':0.15,
                     'JC':0.14,
                     'SP':0.07}

# Add systematic noise tolerance for amplitude a-priori calibration errors
# Start with the SEFD noise (but need sqrt)
# then rescale to ensure that final results respect the stated error budget
systematic_noise = SEFD_error_budget.copy()
for key in systematic_noise.keys():
    systematic_noise[key] = ((1.0+systematic_noise[key])**0.5 - 1.0) * 0.25

# Extra noise added for the LMT, which has much more variability than the a-priori error budget
systematic_noise['LM'] += 0.15

#-------------------------------------------------------------------------------
# Define helper functions
#-------------------------------------------------------------------------------

# Rescale short baselines to excise contributions from extended flux.
# setting zbl < zbl_tot assumes there is an extended constant flux component of zbl_tot-zbl Jy
def rescale_zerobaseline(obs, totflux, orig_totflux, uv_max):
    multiplier = zbl / zbl_tot
    for j in range(len(obs.data)):
        if (obs.data['u'][j]**2 + obs.data['v'][j]**2)**0.5 >= uv_max: continue
        for field in ['vis','qvis','uvis','vvis','sigma','qsigma','usigma','vsigma']:
            obs.data[field][j] *= multiplier

# repeat imaging with blurring to assure good convergence
def converge(major=3, blur_frac=1.0):
    for repeat in range(major):
        init = imgr.out_last().blur_circ(blur_frac*res)
        imgr.init_next = init
        imgr.make_image_I(show_updates=False)

#-------------------------------------------------------------------------------
# Prepare the data
#-------------------------------------------------------------------------------

# Load a single uvfits file
if args.infile2 == '':
    # load the uvfits file
    obs = eh.obsdata.load_uvfits(obsfile)

    # scan-average the data
    # identify the scans (times of continous observation) in the data
    obs.add_scans()

    # coherently average the scans, which can be averaged due to ad-hoc phasing
    obs = obs.avg_coherent(0.,scan_avg=True)

# If two uvfits files are passed as input (e.g., high and low band) then use both datasets,
# but do not form closure quantities between the two datasets
else:
    # load the two uvfits files
    obs1 = eh.obsdata.load_uvfits(obsfile)
    obs2 = eh.obsdata.load_uvfits(args.infile2)

    # Average data based on individual scan lengths
    obs1.add_scans()
    obs2.add_scans()
    obs1 = obs1.avg_coherent(0.,scan_avg=True)
    obs2 = obs2.avg_coherent(0.,scan_avg=True)

    # Add a slight offset to avoid mixed closure products
    obs2.data['time'] += 0.00001

    # concatenate the observations into a single observation object
    obs = obs1.copy()
    obs.data = np.concatenate([obs1.data,obs2.data])

# Estimate the total flux density from the ALMA(AA) -- APEX(AP) zero baseline
zbl_tot   = np.median(obs.unpack_bl('AA','AP','amp')['amp'])
if zbl > zbl_tot:
    print('Warning: Specified total compact flux density ' +
          'exceeds total flux density measured on AA-AP!')

# Flag out sites in the obs.tarr table with no measurements
allsites = set(obs.unpack(['t1'])['t1'])|set(obs.unpack(['t2'])['t2'])
obs.tarr = obs.tarr[[o in allsites for o in obs.tarr['site']]]
obs = eh.obsdata.Obsdata(obs.ra, obs.dec, obs.rf, obs.bw, obs.data, obs.tarr,
                         source=obs.source, mjd=obs.mjd,
                         ampcal=obs.ampcal, phasecal=obs.phasecal)

obs_orig = obs.copy() # save obs before any further modifications

# Rescale short baselines to excize contributions from extended flux
if zbl != zbl_tot:
    rescale_zerobaseline(obs, zbl, zbl_tot, uv_zblcut)

# Order the stations by SNR.
# This will create a minimal set of closure quantities
# with the highest snr and smallest covariance.
obs.reorder_tarr_snr()

#-------------------------------------------------------------------------------
# Pre-calibrate the data
#-------------------------------------------------------------------------------

obs_sc = obs.copy() # From here on out, don't change obs. Use obs_sc to track gain changes
res    = obs.res()  # The nominal array resolution: 1/(longest baseline)

# Make a Gaussian prior image for maximum entropy regularization
# This Gaussian is also the initial image
gaussprior = eh.image.make_square(obs_sc, npix, fov)
gaussprior = gaussprior.add_gauss(zbl, (prior_fwhm, prior_fwhm, 0, 0, 0))

# To avoid gradient singularities in the first step, add an additional small Gaussians
gaussprior = gaussprior.add_gauss(zbl*1e-3, (prior_fwhm, prior_fwhm, 0, prior_fwhm, prior_fwhm))

# Reverse taper the observation: this enforces a maximum resolution on reconstructed features
if reverse_taper_uas > 0:
    obs_sc = obs_sc.reverse_taper(reverse_taper_uas*eh.RADPERUAS)

# Add non-closing systematic noise to the observation
obs_sc = obs_sc.add_fractional_noise(sys_noise)

# Make a copy of the initial data (before any self-calibration but after the taper)
obs_sc_init = obs_sc.copy()

# Self-calibrate the LMT to a Gaussian model
# (Refer to Section 4's "Pre-Imaging Considerations")
print("Self-calibrating the LMT to a Gaussian model for LMT-SMT...")

obs_LMT = obs_sc_init.flag_uvdist(uv_max=2e9) # only consider the short baselines (LMT-SMT)
if reverse_taper_uas > 0:
    # start with original data that had no reverse taper applied.
    # Re-taper, if necessary
    obs_LMT = obs_LMT.taper(reverse_taper_uas*eh.RADPERUAS)

# Make a Gaussian image that would result in the LMT-SMT baseline visibility amplitude
# as estimated in Section 4's "Pre-Imaging Considerations".
# This is achieved with a Gaussian of size 60 microarcseconds and total flux of 0.6 Jy
gausspriorLMT = eh.image.make_square(obs, npix, fov)
gausspriorLMT = gausspriorLMT.add_gauss(0.6, (60.0*eh.RADPERUAS, 60.0*eh.RADPERUAS, 0, 0, 0))

# Self-calibrate the LMT visibilities to the gausspriorLMT image
# to enforce the estimated LMT-SMT visibility amplitude
caltab = eh.selfcal(obs_LMT, gausspriorLMT, sites=['LM'], gain_tol=1.0,
                    method='both', ttype=ttype, caltable=True)

# Spply the calibration solution to the full (and potentially tapered) dataset
obs_sc = caltab.applycal(obs_sc, interp='nearest', extrapolate=True)

#-------------------------------------------------------------------------------
# Reconstruct an image
#-------------------------------------------------------------------------------


# First  Round of Imaging
#-------------------------
print("Round 1: Imaging with visibility amplitudes and closure quantities...")

# Initialize imaging with a Gaussian image
imgr = eh.imager.Imager(obs_sc, gaussprior, prior_im=gaussprior,
                        flux=zbl, data_term=data_term, maxit=maxit,
                        norm_reg=True, systematic_noise=systematic_noise,
                        reg_term=reg_term, ttype=ttype, cp_uv_min=uv_zblcut, stop=stop)

# Imaging
imgr.make_image_I(show_updates=False)
converge()

# Self-calibrate to the previous model (phase-only);
# The solution_interval is 0 to align phases from high and low bands if needed
obs_sc = eh.selfcal(obs_sc, imgr.out_last(), method='phase', ttype=ttype, solution_interval=0.0)


# Second  Round of Imaging
#-------------------------
print("Round 2: Imaging with visibilities and closure quantities...")

# Blur the previous reconstruction to the intrinsic resolution of ~25 uas
init = imgr.out_last().blur_circ(res)

# Increase the weights on the data terms and reinitialize imaging
data_term_intermediate = {'vis':imgr.dat_terms_last()['amp']*10,
                          'cphase':imgr.dat_terms_last()['cphase']*10,
                          'logcamp':imgr.dat_terms_last()['logcamp']*10}

imgr = eh.imager.Imager(obs_sc, init, prior_im=gaussprior, flux=zbl,
                        data_term=data_term_intermediate, maxit=maxit, norm_reg=True,
                        systematic_noise=systematic_noise, reg_term = reg_term, ttype=ttype,
                        cp_uv_min=uv_zblcut, stop=stop)
# Imaging
imgr.make_image_I(show_updates=False)
converge()

# Self-calibrate to the previous model starting from scratch
# phase for all sites; amplitudes for LMT
obs_sc = eh.selfcal(obs_sc_init, imgr.out_last(), method='phase', ttype=ttype)
caltab = eh.selfcal(obs_sc, imgr.out_last(), sites=['LM'], method='both', gain_tol=gain_tol,
                    ttype=ttype, caltable=True)
obs_sc = caltab.applycal(obs_sc, interp='nearest',extrapolate=True)


# Third and Fourth Rounds of Imaging
#-----------------------------------
print("Rounds 3+4: Imaging with visibilities and closure quantities...")

# Increase the data weights before imaging again
data_term_final = {'vis':imgr.dat_terms_last()['vis']*5,
                   'cphase':imgr.dat_terms_last()['cphase']*2,
                   'logcamp':imgr.dat_terms_last()['logcamp']*2}

# Repeat imaging twice
for repeat_selfcal in range(2):
    # Blur the previous reconstruction to the intrinsic resolution of ~25 uas
    init = imgr.out_last().blur_circ(res)

    # Reinitialize imaging now using complex visibilities; common systematic noise
    imgr = eh.imager.Imager(obs_sc, init, prior_im=gaussprior, flux=zbl,
                            data_term=data_term_final, maxit=maxit, norm_reg=True,
                            systematic_noise=0.01, reg_term=reg_term, ttype=ttype,
                            cp_uv_min=uv_zblcut, stop=stop)
    # Imaging
    imgr.make_image_I(show_updates=False)
    converge()

    # Self-calibrate
    caltab = eh.selfcal(obs_sc_init, imgr.out_last(), method='both',
                        sites=['SM','JC','AA','AP','LM','SP','AZ'],
                        ttype=ttype, gain_tol=gain_tol, caltable=True)
    obs_sc = caltab.applycal(obs_sc_init, interp='nearest',extrapolate=True)

#-------------------------------------------------------------------------------
# Output the results
#-------------------------------------------------------------------------------

# Save the reconstructed image
im_out = imgr.out_last().copy()

# If an inverse taper was used, restore the final image
# to be consistent with the original data
if reverse_taper_uas > 0.0:
    im_out = im_out.blur_circ(reverse_taper_uas*eh.RADPERUAS)

# Save the final image
im_out.save_fits(args.outfile)

# Optionally save a pdf of the final image
if args.savepdf:
    pdfout = os.path.splitext(args.outfile)[0] + '.pdf'
    im_out.display(cbar_unit=['Tb'],label_type='scale',export_pdf=pdfout)

# Optionally generate a summary of the final image and associated data consistency metrics
if args.imgsum:

    # Add a large gaussian component to account for the missing flux
    # so the final image can be compared with the original data
    im_addcmp = im_out.add_zblterm(obs_orig, uv_zblcut, debias=True)
    obs_sc_addcmp = eh.selfcal(obs_orig, im_addcmp, method='both', ttype=ttype)

    # Save an image summary sheet
    # Note! these chi2 values will not be identical to what is shown in the paper
    # because they are combining high and low band
    matplotlib.pyplot.close('all')
    outimgsum = os.path.splitext(args.outfile)[0] + '_imgsum.pdf'
    eh.imgsum(im_addcmp, obs_sc_addcmp, obs_orig, outimgsum ,cp_uv_min=uv_zblcut)
