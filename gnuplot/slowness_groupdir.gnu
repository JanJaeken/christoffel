#########################################################################
#                                                                       #
# Copyright (C) 2016 Jan Jaeken <jan.jaeken@gmail.com>                  #
#                                                                       #
# This file is part of Christoffel.                                     #
#                                                                       #
# Christoffel is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# Christoffel is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with Christoffel.  If not, see <http://www.gnu.org/licenses/>.  #
#                                                                       #
#########################################################################

reset
set terminal unknown
plot "slow_secondary.dat" u 0:6;
S_MIN = GPVAL_Y_MIN;
S_MAX = GPVAL_Y_MAX;
plot "fast_secondary.dat" u 0:6;
if (GPVAL_Y_MIN < S_MIN) S_MIN = GPVAL_Y_MIN;
if (GPVAL_Y_MAX > S_MAX) S_MAX = GPVAL_Y_MAX;

plot "primary.dat" u 0:6;
P_MIN = GPVAL_Y_MIN;
P_MAX = GPVAL_Y_MAX;

NUM_VALUES = GPVAL_DATA_X_MAX #This is a float!
NUM_THETA = 0.5*sqrt(NUM_VALUES)
SKIP = floor(NUM_THETA/12.0)

set terminal png truecolor enhanced font 'Arial, 22' size 2250,750
set output "slowness_groupdirection.png"

set palette defined (-1 "red", 0 "white", 1 "blue")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
set hidden3d front
unset border
unset xtics
unset ytics
unset ztics
set cblabel "Slowness (s/km)"

set view 60,120,1.3

set multiplot layout 1,3

set cbrange[1.0/S_MAX:1.0/S_MIN]

set title "Slow Secondary"
splot "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

set title "Fast Secondary"
splot "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

set cbrange[1.0/P_MAX:1.0/P_MIN]

set title "Primary"
splot "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(1.0/$6):(1.0/$6) w p palette ps 2 pt 3 notitle, "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u ((1.0/$6)*sin($1)*cos($2)):((1.0/$6)*sin($1)*sin($2)):((1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(0.03*($13/$11)):(0.03*($14/$11)):(0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (-(1.0/$6)*sin($1)*cos($2)):(-(1.0/$6)*sin($1)*sin($2)):(-(1.0/$6)*cos($1)):(-0.03*($13/$11)):(-0.03*($14/$11)):(-0.03*($15/$11)) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

unset multiplot
unset output

reset

