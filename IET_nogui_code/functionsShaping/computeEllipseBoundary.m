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

function [r_ellipse,z_ellipse,A,B,k] = computeEllipseBoundary(shape_parameters)
            num = shape_parameters.num_points;
            a = shape_parameters.minor_radius;      %minor radius
            Rmax = shape_parameters.RZOU(1);
            Rmin = shape_parameters.RZIN(1);
            Zmax = shape_parameters.RZUP(2);
            Zmin = shape_parameters.RZLO(2);
            Rzmax = shape_parameters.RZUP(1);
            Rzmin = shape_parameters.RZLO(1);
            Zrmax = shape_parameters.RZOU(2);
            Zrmin = shape_parameters.RZIN(2);
            Rgeo = shape_parameters.rgeom;
            delta_u = shape_parameters.upper_triangularity;    %upper triangularity
            delta_l = shape_parameters.lower_triangularity;   %lower triangularity
            
            %% defining new elongations for each quadrant
            Zoff14 = Zrmax;
            Zoff23 = Zrmin;
            
             k_first = (Zmax-Zoff14)/a;
             k_second = (Zmax-Zoff23)/a;
             k_third = (Zoff23-Zmin)/a;
             k_fourth = (Zoff14-Zmin)/a;
            
            r = zeros(4,num);
            z = zeros(4,num);
             A = zeros(4,1);
             B = zeros(4,1);
            
             A(1) = a*(1+delta_u) ;
             B(1) =  k_first*a;
            r(1,:) = linspace(Rmax, Rzmax, num);
            z(1,:) = Zoff14+ B(1)*(1-((r(1,:)-Rzmax)/ A(1)).^2).^0.5;
            
            r(4,:) = linspace(Rzmin, Rmax, num);
             A(4) = a*(1+delta_l);
             B(4) =  k_fourth*a;
            z(4,:) = Zoff14- B(4)*(1-((r(4,:)-Rzmin)/ A(4)).^2).^0.5;
            
             A(2) = a*(1-delta_u) ;
             B(2) =  k_second*a;
            r(2,:) = linspace(Rzmax, Rmin, num);
            z(2,:) = Zoff23+ B(2)*(1-((-r(2,:)+Rzmax)/ A(2)).^2).^0.5;
            
             A(3) = a*(1-delta_l);
             B(3) =  k_third*a;
            r(3,:)=linspace(Rmin, Rzmin, num);
            z(3,:) = Zoff23- B(3)*(1-((-r(3,:)+Rzmin)/ A(3)).^2).^0.5;
            
             r_ellipse=real(r);
             z_ellipse=real(z);
             k = [k_first;k_second;k_third;k_fourth];
        end