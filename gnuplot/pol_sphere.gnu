#########################################################################
#                                                                       #
# Copyright (C) 2016 Jan Jaeken <jan.jaeken@ugent.be>                   #
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

set terminal png truecolor enhanced font Arial 22 size 2250,750
set output "phase_polarization.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
set hidden3d front
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p (km/s)"

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front lw 2
set arrow from first 0,1,0 to first 0,1.4,0 front lw 2
set arrow from first 1,0,0 to first 1.4,0,0 front lw 2

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

set cbrange[S_MIN:S_MAX]

set title "Slow Secondary"
splot "slow_secondary.dat" u (0.98*sin($1)*cos($2)):(0.98*sin($1)*sin($2)):(0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, "slow_secondary.dat" u (-0.98*sin($1)*cos($2)):(-0.98*sin($1)*sin($2)):(-0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(0.03*$8):(0.03*$9):(0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(-0.03*$8):(-0.03*$9):(-0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(0.03*$8):(0.03*$9):(0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(-0.03*$8):(-0.03*$9):(-0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

set title "Fast Secondary"
splot "fast_secondary.dat" u (0.98*sin($1)*cos($2)):(0.98*sin($1)*sin($2)):(0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, "fast_secondary.dat" u (-0.98*sin($1)*cos($2)):(-0.98*sin($1)*sin($2)):(-0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(0.03*$8):(0.03*$9):(0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(-0.03*$8):(-0.03*$9):(-0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(0.03*$8):(0.03*$9):(0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(-0.03*$8):(-0.03*$9):(-0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

set cbrange[P_MIN:P_MAX]

set title "Primary"
splot "primary.dat" u (0.98*sin($1)*cos($2)):(0.98*sin($1)*sin($2)):(0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, "primary.dat" u (-0.98*sin($1)*cos($2)):(-0.98*sin($1)*sin($2)):(-0.98*cos($1)):6:6 w p palette ps 2 pt 3 notitle, "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(0.03*$8):(0.03*$9):(0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):(-0.03*$8):(-0.03*$9):(-0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(0.03*$8):(0.03*$9):(0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle, "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):(-0.03*$8):(-0.03*$9):(-0.03*$10) every SKIP:SKIP w vectors lw 2 lc 2 notitle;

unset multiplot
unset output

reset

