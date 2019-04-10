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

if [ $# -eq 0 ]; then
    echo "usage: $0 [input].uvfits"
    exit 0
fi

in_name=${1%.uvfits}

expect <<EOF
set timeout -1

spawn difmap

expect "*0>"
send -- "@EHT_Difmap ${in_name},CircMask_r30_x-0.002_y0.022,-10,0.5,0.1,2,-1\r"

expect "*0>"
send -- "exit\r"

expect "*quit without saving: "
send -- "\r"

expect eof
EOF
