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
set terminal unknown;
plot "slow_secondary.dat" u 0:12;
S_min = floor(GPVAL_DATA_Y_MIN);
plot "fast_secondary.dat" u 0:12;
S_max = ceil(GPVAL_DATA_Y_MAX);
plot "primary.dat" u 0:12;
P_min = floor(GPVAL_DATA_Y_MIN);
P_max = ceil(GPVAL_DATA_Y_MAX);

if (-S_min > S_max) S_max = -S_min;
if (-P_min > P_max) P_max = -P_min;

set terminal png truecolor enhanced font Arial 22 size 2250,750
set output "group_velocity_relative_eqar.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g - v_{iso} (%)"

set arrow from first 1,0 to first 1.2,0 back
set arrow from first 0,1 to first 0,1.25 back

set label 'X' at first 1.1,0.1 front 
set label 'Y' at first 0.1,1.15 front

set view 0,0,1.1

set multiplot layout 1,3
set title "Slow Secondary"
set cbrange [-S_max:S_max]
splot "slow_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):12:12 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):12:12 w pm3d notitle;
set title "Primary"
set cbrange [-P_max:P_max]
splot "primary.dat" u (sqrt(2)*sin(0.5*$1)*cos($2)):(sqrt(2)*sin(0.5*$1)*sin($2)):12:12 w pm3d notitle;
unset multiplot
unset output

reset

