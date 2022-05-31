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

function [x,y,A,B,k,n_exponent] = computeSuperEllipses(shape_parameters,csi)
            if size(csi,1)==0
            csi = zeros(4,1);
            end
            [r_ellipse,z_ellipse,A,B,k] = computeEllipseBoundary(shape_parameters);

            n_exponent = zeros(4,1);
            
            Rgeo =  shape_parameters.rgeom;
            a =  shape_parameters.minor_radius;      %minor radius
            Rzmax =  shape_parameters.RZUP(1);
            Rzmin =  shape_parameters.RZLO(1);
            Zrmax =  shape_parameters.RZOU(2);
            Zrmin =  shape_parameters.RZIN(2);
            Zoff14 = Zrmax;
            Zoff23 = Zrmin;
            
            epsi = a/Rgeo;     %inverse aspect ratio
            delta_u = (Rgeo-Rzmax)/a;    %upper triangularity
            delta_l = (Rgeo-Rzmin)/a;   %lower triangularity
            
            x(1,:) =  r_ellipse(1,:)-a*(1/epsi-delta_u);                   %first quadrant upper outer
            n_exponent(1) = -log(2)/log(1/(2^0.5) + csi(1)*(1-1/(2^0.5))) ;
            y(1,:) =  B(1)*(1-(x(1,:)/ A(1)).^n_exponent(1)).^(1/n_exponent(1)); %explicit form of equation 14 in base y
            zs(1,:) = y(1,:)+Zoff14;
            
            x(2,:) = a*(1/epsi-delta_u)- r_ellipse(2,:);                   %second quadrant upper inner
            n_exponent(2) = -log(2)/log((0.5)^0.5 + csi(2)*(1-1/(2^0.5)));
            y(2,:) =  B(2)*(1-(x(2,:)/ A(2)).^n_exponent(2)).^(1/n_exponent(2));
            zs(2,:) = y(2,:)+Zoff23;
            
            x(3,:) = a*(1/epsi-delta_l)- r_ellipse(3,:);                   %third quadrant lower inner
            n_exponent(3) = -log(2)/log(1/(2^0.5) + csi(3)*(1-1/(2^0.5)));
            y(3,:) =  B(3)*(1-(x(3,:)/ A(3)).^n_exponent(3)).^(1/n_exponent(3));
            zs(3,:) = Zoff23-y(3,:);
            
            x(4,:) = -a*(1/epsi-delta_l)+ r_ellipse(4,:);                   %fourth quadrant lower outer
            n_exponent(4) = -log(2)/log(1/(2^0.5) + csi(4)*(1-1/(2^0.5)));
            y(4,:) =  B(4)*(1-(x(4,:)/ A(4)).^n_exponent(4)).^(1/n_exponent(4));
            zs(4,:) = Zoff14-y(4,:);
            
            zs_smooth = real(zs);
            
            x = r_ellipse; 
            y=zs_smooth;
            
            
        end