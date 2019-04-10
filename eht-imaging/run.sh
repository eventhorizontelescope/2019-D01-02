#!/usr/bin/env bash
#
# Copyright (C) 2019 The Event Horizon Telescope Collaboration
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

for d in 095 096 100 101; do
    python eht-imaging_pipeline.py \
        -i  ../data/uvfits/SR1_M87_2017_${d}_lo_hops_netcal_StokesI.uvfits \
        -i2 ../data/uvfits/SR1_M87_2017_${d}_hi_hops_netcal_StokesI.uvfits \
        -o                 SR1_M87_2017_${d}.fits
done
