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
set terminal unknown;
plot "slow_secondary.dat" u 0:7;
S_min = floor(GPVAL_DATA_Y_MIN);
plot "fast_secondary.dat" u 0:7;
S_max = ceil(GPVAL_DATA_Y_MAX);
plot "primary.dat" u 0:7;
P_min = floor(GPVAL_DATA_Y_MIN);
P_max = ceil(GPVAL_DATA_Y_MAX);

if (-S_min > S_max) S_max = -S_min;
if (-P_min > P_max) P_max = -P_min;

set terminal png truecolor enhanced font 'Arial, 22' size 2250,750
set output "phase_velocity_relative_cube.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_p - v_{iso} (%)"

set view 60,120,0.95
set arrow from first -1,-1,+1 to first -1,-1,1.3 front
set arrow from first -1,+1,-1 to first -1,1.4,-1 front
set arrow from first +1,-1,-1 to first 1.5,-1,-1 front

set label 'X' at first 1.45,-0.8,-1 front 
set label 'Y' at first -0.65,1.3,-1 front
set label 'Z' at first -1.1,-0.9,1.2 front

set arrow from +1,+1,+1 to -1,+1,+1 nohead front
set arrow from +1,+1,+1 to +1,-1,+1 nohead front
set arrow from +1,+1,+1 to +1,+1,-1 nohead front

set arrow from -1,+1,+1 to -1,+1,-1 nohead front
set arrow from -1,+1,+1 to -1,-1,+1 nohead front

set arrow from +1,-1,+1 to -1,-1,+1 nohead front
set arrow from +1,-1,+1 to +1,-1,-1 nohead front

set arrow from +1,+1,-1 to -1,+1,-1 nohead front
set arrow from +1,+1,-1 to +1,-1,-1 nohead front

set multiplot layout 1,3
set cbrange[-S_max:S_max]
set title "Slow Secondary"
splot "slow_secondary.dat" u 3:4:5:7 w pm3d notitle, "slow_secondary.dat" u (-$3):(-$4):(-$5):7 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u 3:4:5:7 w pm3d notitle, "fast_secondary.dat" u (-$3):(-$4):(-$5):7 w pm3d notitle;
set cbrange[-P_max:P_max]
set title "Primary"
splot "primary.dat" u 3:4:5:7 w pm3d notitle, "primary.dat" u (-$3):(-$4):(-$5):7 w pm3d notitle;
unset multiplot
unset output

reset

