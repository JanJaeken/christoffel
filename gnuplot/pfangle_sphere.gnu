#########################################################################
#                                                                       #
# Copyright (C) 2016 Jan Jaeken <jan.jaeken@gmail.com>                   #
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
set terminal png truecolor enhanced font Arial 22 size 2250,750
set output "powerflow_sphere.png"

set palette defined (0 "black", 1 "green")
set pm3d depthorder explicit
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "PF angle (Deg)"

set view 60,120,1.3
set arrow from first 0,0,1 to first 0,0,1.35 front
set arrow from first 0,1,0 to first 0,1.4,0 front
set arrow from first 1,0,0 to first 1.4,0,0 front

set label 'X' at first 1.7,0,0 front 
set label 'Y' at first 0.2,1.35,0 front
set label 'Z' at first 0,0.07,1.3 front

set multiplot layout 1,3

set title "Slow Secondary"
splot "slow_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):16 w pm3d notitle, "slow_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):16 w pm3d notitle

set title "Fast Secondary"
splot "fast_secondary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):16 w pm3d notitle, "fast_secondary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):16 w pm3d notitle

set title "Primary"
splot "primary.dat" u (sin($1)*cos($2)):(sin($1)*sin($2)):(cos($1)):16 w pm3d notitle, "primary.dat" u (-sin($1)*cos($2)):(-sin($1)*sin($2)):(-cos($1)):16 w pm3d notitle

unset multiplot
unset output

reset

