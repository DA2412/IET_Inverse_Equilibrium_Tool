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
function [zs_Xpoint] = computeSuperellipsesXpoint(shape_parameters,csi)
            if size(csi,1)==0
            csi = zeros(4,1);
            end
            [r_ellipse,z_ellipse,A,B,k,n_exponent] = computeSuperEllipses(shape_parameters,csi);
            k_first = k(1);
            k_second = k(2);
            k_third = k(3);
            k_fourth = k(4);
            
            Rgeo = shape_parameters.rgeom;
            a = shape_parameters.minor_radius;      %minor radius
            Rzmax = shape_parameters.RZUP(1);
            Rzmin = shape_parameters.RZLO(1);
            Zrmax = shape_parameters.RZOU(2);
            Zrmin = shape_parameters.RZIN(2);
            Zoff14 = Zrmax;
            Zoff23 = Zrmin;
            
            epsi = a/Rgeo;     %inverse aspect ratio
            delta_u = (Rgeo-Rzmax)/a;    %upper triangularity
            delta_l = (Rgeo-Rzmin)/a;   %lower triangularity
            
            x0 = [ A(:,1)  B(:,1) n_exponent(:,1)];
            first=real(fsolve(@(s)constraints1(a,delta_u, k_first,csi,s),x0(1,:)));%first quad
            second=real(fsolve(@(s)constraints2(a,delta_u, k_second,csi,s),x0(2,:)));%second quad
            third=real(fsolve(@(s)constraints3(a,delta_l, k_third,csi,s),x0(3,:)));%third quad
            fourth=real(fsolve(@(s)constraints4(a,delta_l, k_fourth,csi,s),x0(4,:)));%fourth quad
            
            alpha=[first(1);second(1);third(1);fourth(1)];
            beta=[first(2);second(2);third(2);fourth(2)];
            n_exponent_Xpoint=[first(3);second(3);third(3);fourth(3)];
            
            X(1,:)= + r_ellipse(1,:)-a*(1/epsi+1)+alpha(1);                %first quadrant upper outer
            y(1,:) =  beta(1)*(1-(X(1,:)/alpha(1)).^n_exponent_Xpoint(1)).^(1/n_exponent_Xpoint(1)) ; %explicit form of equation 14 in base y
            zs(1,:) = y(1,:)+Zoff14;
            
            
            X(2,:)=a*(1/epsi-1)+alpha(2)- r_ellipse(2,:);
            y(2,:) = beta(2)*(1-(X(2,:)/alpha(2)).^n_exponent_Xpoint(2)).^(1/n_exponent_Xpoint(2)) ;
            zs(2,:) = y(2,:)+Zoff23;
            %plot(r(2,:),zs(2,:),'-r','linewidth',2); %plot second quadrant ellipse
            
            %r(3,:)=linspace(Rmin,Rzmin+5/100*Rzmin,num);
            X(3,:)=a*(1/epsi-1)+alpha(3)- r_ellipse(3,:); %sostituito r*a con r
            y(3,:) =  beta(3)*(1-(X(3,:)/alpha(3)).^n_exponent_Xpoint(3)).^(1/n_exponent_Xpoint(3)) ;
            zs(3,:) = Zoff23-y(3,:);
            
            %r(4,:) = linspace(Rmax, Rzmin-5/100*Rzmin, num);
            X(4,:)= + r_ellipse(4,:) -a*(1/epsi+1) +alpha(4);
            %X(4,:)= +r(4,:) -a*(1/epsi-1) -alpha(4);
            y(4,:) =  beta(4)*(1-(X(4,:)/alpha(4)).^n_exponent_Xpoint(4)).^(1/n_exponent_Xpoint(4)) ;
            zs(4,:) = Zoff14 - y(4,:);
            
            zs_Xpoint = real(zs);
            
            
        end