% IET (Inverse Equilibrium Tool) is a MATLAB computational tool constituted by three main modules,
% each with a dedicated GUI. It allows to compute the coil currents needed to obtain a predetermined 
% plasma shape with defined plasma global parameters (i.e. total plasma current and total poloidal 
% magnetic flux at the boundary) by solving a constrained minimization problem. 
% 
%     Copyright (C) 2019 Domenico Abate - domenico.abate@igi.cnr.it
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

function  [sp]=computeShapeParameters(R_boundary, Z_boundary)
sp.num_points = length(R_boundary);

[sp.RZOU(1), iout] = max(R_boundary);
sp.RZOU(2) = Z_boundary(iout);

[sp.RZIN(1), iin] = min(R_boundary);
sp.RZIN(2) = Z_boundary(iin);

[sp.RZUP(2), iup] = max(Z_boundary);
sp.RZUP(1) = R_boundary(iup);

[sp.RZLO(2), ilo] = min(Z_boundary);
sp.RZLO(1) = R_boundary(ilo);

sp.rgeom=.5*( sp.RZOU(1)+ sp.RZIN(1));
 sp.zgeom=.5*( sp.RZOU(2)+ sp.RZIN(2));
sp.major_radius= sp.rgeom;
 sp.minor_radius=sp.major_radius- sp.RZIN(1);
 sp.upper_elongation=( sp.RZUP(2)- sp.zgeom)/( sp.RZOU(1)- sp.rgeom);
 sp.lower_elongation=( sp.zgeom- sp.RZLO(2))/( sp.RZOU(1)- sp.rgeom);
 sp.elongation=( sp.RZUP(2)- sp.RZLO(2))/( sp.RZOU(1)- sp.RZIN(1));
sp.horizontal_elongation = 1/ sp.elongation;
 sp.upper_triangularity=( sp.rgeom- sp.RZUP(1))/( sp.rgeom- sp.RZIN(1));
 sp.lower_triangularity=( sp.rgeom- sp.RZLO(1))/( sp.rgeom- sp.RZIN(1));
 sp.triangularity=( sp.rgeom-.5* sp.RZUP(1)-.5* sp.RZLO(1))/( sp.rgeom- sp.RZIN(1));

 end