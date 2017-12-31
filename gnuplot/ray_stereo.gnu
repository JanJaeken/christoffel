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
set terminal png truecolor enhanced font 'Arial, 22' size 2250,750
set output "ray_stereo.png"

set palette defined (-1 "blue", 0 "white", 1 "red")
set view equal xyz
set xyplane 0
unset border
unset xtics
unset ytics
unset ztics
set cblabel "v_g (km/s)"
set pm3d depthorder explicit

#set view map
set view 0,0,1.3

set multiplot layout 1,3
set title "Slow Secondary"
splot "slow_secondary.dat" u ($13/($11+$15)):($14/($11+$15)):11:11 w pm3d notitle;
set title "Fast Secondary"
splot "fast_secondary.dat" u ($13/($11+$15)):($14/($11+$15)):11:11 w pm3d notitle;
set title "Primary"
splot "primary.dat" u ($13/($11+$15)):($14/($11+$15)):11:11 w pm3d notitle;
unset multiplot
unset output

reset

