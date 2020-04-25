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
load('../../Examples/RFX_example2/OptimizationTool_Inputs_Outputs/output/IET_MultiObjective_AOF.mat')
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
MAPE_lambda = AOF_solution.MAPE_lambda;
RMSE_lambda = AOF_solution.RMSE_lambda;
currents = AOF_solution.currents;
PFront = AOF_solution.ParetoFrontData; 

ioptim = find(MAPE_lambda<1.2& RMSE_lambda<0.75);

figure;
plot(MAPE_lambda,RMSE_lambda,'bo','markersize',13,'linewidth',1.5)
hold('on');
plot(PFront(:,1),PFront(:,2),'m*','linewidth',1.5);
grid('on');
axis('equal');
xlabel('MAPE','interpreter','latex')
ylabel('RMSE','interpreter','latex')
plot(MAPE_lambda(ioptim),RMSE_lambda(ioptim),'cd','markersize',12,'linewidth',2)
legend('AOF solutions','Pareto Front','Selected solution','location','best')
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',12)
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)


Psi_tot_Lambda_Scelto=GG_plasma_nodes.psi*plasma.current+...
   GG_coils_nodes.psi*KONNAX'*currents(ioptim,:)';



PSI_LambdaScelto = griddata(p(:,1),p(:,2),Psi_tot_Lambda_Scelto,RR,ZZ);
err_patch_Lambda_Scelto = (Psi_tot_Lambda_Scelto- Psi_tot_direct)./(Psi_tot_direct)*100;


figure
subplot(1,2,1)
contour(RR,ZZ,PSI_LambdaScelto,10,'linewidth',2)
hold on
plot(R_boundary_finale,Z_boundary_finale,'--r','linewidth',2);
grid on
axis equal
colormap(jet)
c= colorbar('vert');
%c.Label.String = '\Psi';
xlabel('R [m]')
ylabel('Z [m]')
title(['Selected solution'])
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)

subplot(1,2,2)
JP.faces=t;
JP.vertices=p;
JP.facevertexcdata=err_patch_Lambda_Scelto;
patch(JP,'facecolor','interp','edgecolor','none');
hold on
plot(R_boundary_finale,Z_boundary_finale,'--r','linewidth',2);
axis equal
grid on
colormap(jet)
c = colorbar('vert');     %c.Label.String = '\epsilon_{\Psi}';
xlabel('R [m]')
ylabel('Z [m]')
title('Selected solution \epsilon_{\Psi} [%]') 
set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)
