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

REPO=https://github.com/eventhorizontelescope/2019-D01-01/raw/master
FILE=EHTC_FirstM87Results_Apr2019_uvfits.tgz

echo "Download uvfits files for the EHT first M87 results"
wget $REPO/$FILE

echo "Extract uvfits files to disk"
mkdir -p data
tar -vzxf $FILE -C data --strip-components=1

echo "Run imaging pipelines"
for d in difmap/ eht-imaging/ smili/; do
    pushd $d
    ./run.sh
    popd > /dev/null
done
