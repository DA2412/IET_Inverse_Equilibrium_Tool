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

clear all
close all
clc
load('./36922BoundaryGUI.mat') %load file from eqdisk
addpath functionsShaping
addpath functionsOptim
addpath functionsEquil
addpath('distmeshModNoGui');

mu0 = 4*pi*1e-7;
R_boundary = geometry.boundaryCoordinates(:,1);
Z_boundary = geometry.boundaryCoordinates(:,2);

%% define plasma boundary
disp('------------------------------------DEFINING PLASMA BOUNDARY------------------------------------')
csi=[];
csi(1,1)=-0.5;
csi(2,1)=-0.3;
csi(3,1)=0;
csi(4,1)= 0;

[shape_parameters]=computeShapeParameters(R_boundary, Z_boundary);
[x,y,A,B,n_exponent] = computeSuperEllipses(shape_parameters,csi);
[zs_Xpoint] = computeSuperellipsesXpoint(shape_parameters,csi); %for geometrical X-point calculation

figure;
plot(R_boundary,Z_boundary,'r-','linewidth',2.5);   %plot real boundary
hold on,grid on, axis equal
plot(x(1,:),y(1,:),'b-','linewidth',2); %plot first quadrant ellipse
plot(x(2,:),y(2,:),'b-','linewidth',2); %plot second quadrant ellipse
plot(x(1,:),zs_Xpoint(1,:),'k-','linewidth',2); %for geometrical X-point
plot(x(2,:),zs_Xpoint(2,:),'k-','linewidth',2); %for geometrical X-point
plot(x(3,:),y(3,:),'b-','linewidth',2);%plot third quadrant ellipse
plot(x(4,:),y(4,:),'b-','linewidth',2);%plot fourth quadrant ellipse

R_boundary = [x(1,:),x(2,:),x(3,:),x(4,:)]'; 
Z_boundary = [y(1,:),y(2,:),y(3,:),y(4,:)]'; %change y to zs_Xpoint for geometrical X-point 

%% define mesh
disp('------------------------------------DEFINING PLASMA MESH------------------------------------')
b_box = [min(R_boundary) - 5,min(Z_boundary) - 1; max(max(R_boundary)) + 5,max(Z_boundary) + 1]; %Bounding box [xmin,ymin; xmax,ymax]
h0 = 0.02; %0.02
tic
[p,t]=distmesh2d(@dpoly,@huniform,h0,b_box,[],[R_boundary Z_boundary]);
toc
e = boundedges(p,t);
b = unique(e);
node = p;
triplot(t,p(:,1),p(:,2))
hold on
[Area_duale, kk_nu] = compute_dual_mesh_line(p, t, mu0);
plot(R_boundary,Z_boundary,'-r','linewidth',2);
axis equal

OutputShape.p = p;
OutputShape.t = t;
OutputShape.b = b;
OutputShape.R_boundary_finale = R_boundary;
OutputShape.Z_boundary_finale = Z_boundary;
OutputShape.kk_nu = kk_nu;
OutputShape.Area_duale = Area_duale;

%% fixed boundary equilibrium calculation via cell method
disp('------------------------------------COMPUTING PLASMA EQUILIBRIUM------------------------------------')
% numerical inputs
Ipla = 50e3;
psiBoundary = -0.68;
JpType = 'Blum';
R0=     1.995;
Beta0= 0.3;
AlphaM= 1.5;
AlphaN= 1.001;
tol = 1e-12;
ngrid = 400;

IET_code_fb

figure;
semilogy(Residual_Picard,'b-x','linewidth',2)
 hold('on')
  semilogy(scarto,'r-o','linewidth',2)
 grid('on')
 xlabel('Iterations')
 axis([0 iter 0 max(Residual_Picard)])
 ylabel('Picard convergence')
 legend('error','residual')

figure;
contour(RR,ZZ,EqSolution.PSI_final,100)
 hold('on')
 axis equal
plot(R_boundary_finale,Z_boundary_finale,'-r','linewidth',2);
title('\Psi final equilibrium')
colormap('jet')
colorbar('vert')

figure;
JP.faces=tri;
JP.vertices=node;
JP.facevertexcdata=EqSolution.Jphi_final./(max(EqSolution.Jphi_final));
patch(JP,'facecolor','interp','edgecolor','none');
hold('on')
plot(R_boundary_finale,Z_boundary_finale,'-r','linewidth',2);
grid('on')
axis('equal')
title('Jp final')
colormap('jet')
colorbar('vert')

%% Inverse Equilibrium: coil current calculation via least squares 
disp('------------------------------------COMPUTING COIL CURRENTS------------------------------------')
load ('RFX_active_coils.mat');

IET_ie

figure;
ax1 = gca;
hold on
contour(ax1,RR_g,ZZ_g,PSI_LS,100,'linewidth',2)
contour(ax1,RR_g,ZZ_g,PSI_LS,linspace((1-dim)*psi_b,(1+dim)*psi_b,10), 'r', 'LineWidth' ,1);
contour(gca,RR_g,ZZ_g,PSI_LS,[psi_b psi_b], 'r', 'LineWidth' ,2);
plot(R_fw,Z_fw,'k','linewidth',2)
axis( 'equal');
colormap( jet);
colorbar( 'vert');
xlabel( 'R [m]')
ylabel( 'Z [m]')
title( 'Least Square solution')
plot(R_boundary_finale,Z_boundary_finale,'-m','LineWidth',2);
xlim( [min(min(RR_g)) max(max(RR_g))])
ylim( [min(min(ZZ_g)) max(max(ZZ_g))])
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)
