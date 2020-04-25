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
addpath('funOptim');
load('../../Examples/RFX_example2/OptimizationTool_Inputs_Outputs/output/IET_MultiObjective_NC.mat');
GG_plasma_nodes = OptimDataIn.GG_plasma_nodes;
GG_coils_nodes = OptimDataIn.GG_coils_nodes;
KONNAX = OptimDataIn.kon;
plasma = OptimDataIn.plasma;
p = OptimDataIn.p;
t = OptimDataIn.t;
RR = OptimDataIn.RR;
ZZ = OptimDataIn.ZZ;
R_boundary_finale = OptimDataIn.R_boundary_finale;
Z_boundary_finale = OptimDataIn.Z_boundary_finale;
Psi_tot_direct = OptimDataIn.psi_tot_direct;
%% selecting solution
mu1_realbest = NC_solution.mu1_realbest;
mu2_realbest = NC_solution.mu2_realbest;
currents = NC_solution.currents;

[index]= find(mu1_realbest<1.05 & mu2_realbest<0.71);

figure;
plot(mu1_realbest,mu2_realbest,'bo','markersize',13,'linewidth',1)
hold('on');
grid('on');
axis('square');
xlabel('MAPE')
ylabel('RMSE')
plot(mu1_realbest(index),mu2_realbest(index),'cd','markersize',16,'linewidth',2)
legend('Pareto Points','Selected Solution')
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',12)

Psi_tot_NC=GG_plasma_nodes.psi*plasma.current+...
    GG_coils_nodes.psi*KONNAX'*currents(:,index);

figure
subplot(1,2,1)
PSI_1star = griddata(p(:,1),p(:,2),Psi_tot_NC,RR,ZZ);
contour(RR,ZZ,PSI_1star,10,'linewidth',2)
hold on
plot(R_boundary_finale,Z_boundary_finale,'--r','linewidth',2);
grid on
axis equal
colormap(jet)
colorbar('vert');
xlabel('R [m]')
ylabel('Z [m]')
title('NC')
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)

%error
err_patch_NC = (Psi_tot_NC - Psi_tot_direct)./(Psi_tot_direct)*100;

subplot(1,2,2)
JP.faces=t;
JP.vertices=p;
JP.facevertexcdata=err_patch_NC;
patch(JP,'facecolor','interp','edgecolor','none');
hold on
plot(R_boundary_finale,Z_boundary_finale,'--r','linewidth',2);
axis equal
grid on
colormap(jet)
colorbar('vert');xlabel('R [m]')
ylabel('Z [m]')
title('\epsilon NC %')
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)

