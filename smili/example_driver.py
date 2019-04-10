#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Example driver for the SMILI M87 Stokes I Imaging Pipeline for EHT observations in April 2017

Authors: The Event Horizon Telescope Collaboration et al.
Date: April 10, 2019
Primary Reference: The Event Horizon Telescope Collaboration, et al. 2019d, ApJL, 875, L4 (M87 Paper IV)
Data Product Code: 2019-D01-02

Brief Description:
The script is an example driver of the SMILI M87 Stokes I Imaging Pipeline
(smili_imaging_pipeline.py) for EHT observations in April 2017 attached in the
same directory. For more detail instructions for smili_imaging_pipeline.py
please read the help document associated in the imaging script
"python smili_imaging_pipeline.py --help".

This sript reconstructs images on fiducial parameters (see Section 6 of M87
Paper IV) from calibratd uvfits files released in the Data Product Release
(2019-D01-01; see also M87 Paper III). You can run it by
"python example_driver.py --uvfitsdir <xxxx/uvfits> --nproc <number of the processors>."
For details, please have a look at the help document by
"python example_driver.py --help".

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
import os
import argparse
import glob
import tqdm
from multiprocessing import Pool


#-------------------------------------------------------------------------------
# Loading command line arguments
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--uvfitsdir',metavar='input uvfits directory',type=str,default='../../data/uvfits/',
                    help='input uvfits directory of the data product release 2019-XX-XXXX')
parser.add_argument('--outputdir',metavar='output fits file',type=str,default="./smili_reconstructions",
                    help='output directory. The default is "./smili_reconstructions".')
parser.add_argument('--pipeline',metavar='pipeline script',type=str,default="smili_imaging_pipeline.py",
                    help='the filename of the pipeline scripts. The default is "smili_imaging_pipeline.py"')
parser.add_argument("--nproc",metavar="Number of parallel processors. Default=4.",type=int,default=4)
args = parser.parse_args()

nproc = args.nproc
uvfitsdir = args.uvfitsdir
outputdir = args.outputdir
pipeline = args.pipeline


#-------------------------------------------------------------------------------
# Preparation
#-------------------------------------------------------------------------------
# create the output directory
if not os.path.isdir(outputdir):
    command = "mkdir -p %s"%(outputdir)
    print(command)
    os.system(command)

# get uvfits files
uvfitsfiles = sorted(glob.glob(os.path.join(uvfitsdir, "*_hops_netcal_StokesI.uvfits")))
if len(uvfitsfiles) < 1:
    raise ValueError("No input uvfits files in %s"%(uvfitsdir))

# copy uvfits files
for uvfitsfile in uvfitsfiles:
    uvfitsfilebase = os.path.basename(uvfitsfile)
    command = "cp -p %s %s"%(uvfitsfile, os.path.join(outputdir, uvfitsfilebase))
    print(command)
    os.system(command)

# generating commands
commands = []
for uvfitsfile in uvfitsfiles:
    uvfitsfilebase = os.path.basename(uvfitsfile)

    # check observing errors
    obsdate = None
    for day in [95,96,100,101]:
        if "_%03d_"%(day) in uvfitsfilebase:
            obsdate = day - 90
            break
    if obsdate is None:
        raise ValueError("Apparently there is no obsday information in the filename of the uvfits: %s"%(uvfitsfilebase))

    command = []
    command.append("python")
    command.append(pipeline)
    command.append("-i %s"%(os.path.join(outputdir, uvfitsfilebase)))
    command.append("--day %d"%(obsdate))
    command.append("--nproc 1")
    commands.append(" ".join(command))

#-------------------------------------------------------------------------------
# Running the pipeline
#-------------------------------------------------------------------------------
# OMP NUM THREADS must be set to 1, so that each command will use only a single processor.
os.environ['OMP_NUM_THREADS'] = '1'

# running imaging in parallel
pool = Pool(processes=nproc)
for _ in tqdm.tqdm(pool.imap_unordered(os.system,commands),total=len(commands)):
    pass
exit()
