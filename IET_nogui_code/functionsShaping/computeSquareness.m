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


function [csi] = computeSquareness(R_boundary_ref,Z_boundary_ref)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
            [shape_parameters]=computeShapeParameters(R_boundary_ref, Z_boundary_ref);
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
            Rgeo =  shape_parameters.rgeom;

            Zoff14 = Zrmax;
            Zoff23 = Zrmin;
            
            
            [r_ellipse,z_ellipse] = computeEllipseBoundary(shape_parameters);
            num = length(r_ellipse);
                E = zeros(4,2);   %extremal quadrant point matrix
                O = zeros(4,2);   %center quadrant point matrix
                C = zeros(4,2);   %intersection point between diagonal and ellipse matrix
                D = zeros(4,2);   %intersection point between diagonal and real shape matrix
                Rdiag = zeros(4,num); %R diagonal coordinates
                Zdiag = zeros(4,num); %Z diagonal coordinates
                
                E(1,:) = [Rmax Zmax];               %first quadrant
                O(1,:) = [Rzmax Zoff14];
                Rdiag(1,:) = linspace(O(1,1),Rmax,num);
                Zdiag(1,:) = O(1,2)+(E(1,2)-O(1,2))/(E(1,1)-O(1,1))*(Rdiag(1,:)-O(1,1));
                [RC,ZC] = intersections(Rdiag(1,:),Zdiag(1,:),r_ellipse(1,:),z_ellipse(1,:));  %calculation of intersections between ellipse and diagonal
                C(1,1) = RC;
                C(1,2) = ZC;
                [RD,ZD] = intersections(Rdiag(1,:),Zdiag(1,:),R_boundary_ref,Z_boundary_ref); %calculation on fintersections between real shape and diagonal
                D(1,1) = RD;
                D(1,2) = ZD;
                
                E(2,:) = [Rmin Zmax];   %second quadrant
                O(2,:) = [Rzmax Zoff23];
                Rdiag(2,:) = linspace(Rmin,O(2,1),num);
                Zdiag(2,:) = O(2,2)+(E(2,2)-O(2,2))/(E(2,1)-O(2,1))*(Rdiag(2,:)-O(2,1));
                [RC,ZC] = intersections(Rdiag(2,:),Zdiag(2,:),r_ellipse(2,:),z_ellipse(2,:));  %calculation of intersections between ellipse and diagonal
                C(2,1) = RC;
                C(2,2) = ZC;
                [RD,ZD] = intersections(Rdiag(2,:),Zdiag(2,:),R_boundary_ref,Z_boundary_ref); %calculation on fintersections between real shape and diagonal
                D(2,1) = RD;
                D(2,2) = ZD;
                
                E(3,:) = [Rmin Zmin];       %third quadrant
                O(3,:) = [Rzmin Zoff23];
                Rdiag(3,:) = linspace(Rmin,O(3,1),num);
                Zdiag(3,:) = O(3,2)+(E(3,2)-O(3,2))/(E(3,1)-O(3,1))*(Rdiag(3,:)-O(3,1)); %defining diagonal between point O and point E
                [RC,ZC] = intersections(Rdiag(3,:),Zdiag(3,:),r_ellipse(3,:),z_ellipse(3,:));  %calculation of intersections between ellipse and diagonal
                C(3,1) = RC;
                C(3,2) = ZC;
                [RD,ZD] = intersections(Rdiag(3,:),Zdiag(3,:),R_boundary_ref,Z_boundary_ref); %calculation on fintersections between real shape and diagonal
                D(3,1) = RD;
                D(3,2) = ZD;
                
                E(4,:) = [Rmax Zmin];   %fourth quadrant
                O(4,:) = [Rzmin Zoff14];
                Rdiag(4,:) = linspace(O(4,1),Rmax,num);
                Zdiag(4,:) = O(4,2)+(E(4,2)-O(4,2))/(E(4,1)-O(4,1))*(Rdiag(4,:)-O(4,1)); %defining diagonal between point O and point E
                [RC,ZC] = intersections(Rdiag(4,:),Zdiag(4,:),r_ellipse(4,:),z_ellipse(4,:));  %calculation of intersections between ellipse and diagonal
                C(4,1) = RC;
                C(4,2) = ZC;
                [RD,ZD] = intersections(Rdiag(4,:),Zdiag(4,:),R_boundary_ref,Z_boundary_ref); %calculation on fintersections between real shape and diagonal
                D(4,1) = RD;
                D(4,2) = ZD;
                
                %% distances calculation for each quadrant
                dOE(1) = pdist2(O(1,:),E(1,:),'euclidean');%first quadrant
                dOD(1) = pdist2(O(1,:),D(1,:),'euclidean');
                dOC(1) = pdist2(O(1,:),C(1,:),'euclidean');
                
                dOE(2) = pdist2(O(2,:),E(2,:),'euclidean');%second quadrant
                dOD(2) = pdist2(O(2,:),D(2,:),'euclidean');
                dOC(2) = pdist2(O(2,:),C(2,:),'euclidean');
                
                dOE(3) = pdist2(O(3,:),E(3,:),'euclidean');%third quadrant
                dOD(3) = pdist2(O(3,:),D(3,:),'euclidean');
                dOC(3) = pdist2(O(3,:),C(3,:),'euclidean');
                
                dOE(4) = pdist2(O(4,:),E(4,:),'euclidean');%fourth quadrant
                dOD(4) = pdist2(O(4,:),D(4,:),'euclidean');
                dOC(4) = pdist2(O(4,:),C(4,:),'euclidean');
                              
                
                csi(1,1) = (dOD(1)-dOC(1))/(dOE(1)-dOC(1));
                csi(2,1) = (dOD(2)-dOC(2))/(dOE(2)-dOC(2));
                csi(3,1) = (dOD(3)-dOC(3))/(dOE(3)-dOC(3));
                csi(4,1) = (dOD(4)-dOC(4))/(dOE(4)-dOC(4));
                
%                                 figure;
%                 hold on
%                 plot(R_boundary_ref,Z_boundary_ref,'r-','linewidth',2.5);   %plot real boundary
%                 hold on,grid on, axis equal
%                 plot(O(1,1),O(1,2),'o','markersize',10,'markerfacecolor','k');
%                 plot(O(2,1),O(2,2),'o','markersize',10,'markerfacecolor','b');
%                 plot(O(3,1),O(3,2),'o','markersize',10,'markerfacecolor','g');
%                 plot(O(4,1),O(4,2),'o','markersize',10,'markerfacecolor','r');
%                 text(O(1,1),O(1,2),'  O1','FontSize',11);
%                 text(O(2,1),O(2,2),'  O2','FontSize',11);
%                 text(O(3,1),O(3,2),'O3   ','HorizontalAlignment','right','FontSize',11);
%                 text(O(4,1),O(4,2),'O4   ','HorizontalAlignment','right','FontSize',11);
%                 
%                 plot(E(1,1),E(1,2),'.')
%                 text(E(1,1),E(1,2),' E1','FontSize',12)
%                 plot(Rdiag(1,:),Zdiag(1,:),'k')
%                 plot (C(1,1),C(1,2),'.')
%                 text(C(1,1),C(1,2),'C1  ','HorizontalAlignment','right','FontSize',12)
%                 plot (D(1,1),D(1,2),'.')
%                 text(D(1,1),D(1,2),'   D1','FontSize',12)
%                 
%                 plot(E(2,1),E(2,2),'.')
%                 text(E(2,1),E(2,2),' E2','FontSize',12)
%                 plot(Rdiag(2,:),Zdiag(2,:),'k')
%                 plot (C(2,1),C(2,2),'.')
%                 text(C(2,1),C(2,2),'C2  ','HorizontalAlignment','right','FontSize',12)
%                 plot (D(2,1),D(2,2),'.')
%                 text(D(2,1),D(2,2),'   D2','FontSize',12)
%                 
%                 plot(E(3,1),E(3,2),'.')
%                 text(E(3,1),E(3,2),'  E3','FontSize',12)
%                 plot(Rdiag(3,:),Zdiag(3,:),'k')
%                 plot (C(3,1),C(3,2),'.')
%                 text(C(3,1),C(3,2),'C3  ','HorizontalAlignment','right','FontSize',12)
%                 plot (D(3,1),D(3,2),'.')
%                 text(D(3,1),D(3,2),'   D3','FontSize',12)
%                 
%                 plot(E(4,1),E(4,2),'.')
%                 text(E(4,1),E(4,2),' E4','FontSize',12)
%                 plot(Rdiag(4,:),Zdiag(4,:),'k')
%                 plot (C(4,1),C(4,2),'.')
%                 text(C(4,1),C(4,2),'C4  ','HorizontalAlignment','right','FontSize',12)
%                 plot (D(4,1),D(4,2),'.')
%                 text(D(4,1),D(4,2),'    D4','FontSize',12)
%                 
%                 plot(r_ellipse,z_ellipse,'b.');

end

