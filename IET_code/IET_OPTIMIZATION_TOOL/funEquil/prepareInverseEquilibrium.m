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
function [OptimDataIn] = prepareInverseEquilibrium(app)
psi_tot_direct = app.OutputEquil.Equil.psi_tot_direct;
Iphi_plasma = app.OutputEquil.Equil.Iphi_plasma;
psiBoundary = app.OutputEquil.Equil.psiBoundary;
rr = app.OutputEquil.OutputShape.node(:,1);
zz = app.OutputEquil.OutputShape.node(:,2);
p = app.OutputEquil.OutputShape.node;
t= app.OutputEquil.OutputShape.tri;
Area_duale = app.OutputEquil.OutputShape.Area_duale;
nodi_BCs = app.OutputEquil.OutputShape.nodi_BCs;
coils = app.coils;

%JPHI=zeros(size(rr,1),1);
%JPHI=Iphi_plasma ./ Area_duale;
%Jphi_interp=pdeInterpolant(p',tnew,JPHI);

% Triangle point indices
it1=t(:,1);
it2=t(:,2);
it3=t(:,3);
% Find centroids of triangles
rbar=(rr(it1)+rr(it2)+rr(it3))/3;
zbar=(zz(it1)+zz(it2)+zz(it3))/3;

%Jp_b = Jphi_interp.evaluate(rbar,zbar)';

% Calculating plasma boundary flux contribution
Jp_n = Iphi_plasma./Area_duale;                      %plasma current density nodes

Jp_b = pdeintrp(p',t',Jp_n); %interpolating J over the triangles centroids

Area_primale=pdetrg(p',t');   %calculating area of triangles

Ip_b = Jp_b .* Area_primale;    % calculating total current per element

space.RR=rr(nodi_BCs);  %selecting only boundary points
space.ZZ=zz(nodi_BCs);

plasma.R=rbar;  %creating new struct variable 'plasma'
plasma.Z=zbar;
plasma.current=Ip_b';
GG_plasma= fun_Field_Loop(plasma, space);

%psi_plasma = GG_plasma.psi * plasma.II;  %calculating plasma contribution to psi boundary
psi_plasma = GG_plasma.psi * plasma.current;
psi_coils = psi_tot_direct(nodi_BCs) - psi_plasma;    %calculating coils contribution to psi boundary


%% Plot psi boundary values over boundary angle coordinate

r0 = sum(rr(nodi_BCs))/(length(nodi_BCs));
z0 = sum(zz(nodi_BCs))/(length(nodi_BCs));
dr = rr(nodi_BCs)-r0;
dz = zz(nodi_BCs)-z0;

theta = atan2(dz,dr);
[theta_ord itheta]=sort(theta);
%theta_ord=radtodeg(theta_ord);

% figure(20)
% hold on
% grid on
plot(app.psiBoundaryPlot,theta_ord,psi_plasma(itheta),'bo','markersize',15);
hold(app.psiBoundaryPlot,'on');
grid(app.psiBoundaryPlot,'on');
plot(app.psiBoundaryPlot,theta_ord,psi_coils(itheta),'ro','markersize',15);
plot(app.psiBoundaryPlot,theta_ord,psi_coils(itheta)+psi_plasma(itheta),'go','markersize',15);
plot(app.psiBoundaryPlot,theta_ord,psiBoundary*ones(length(theta_ord),1),'m*','markersize',15);
xlabel(app.psiBoundaryPlot,'Poloidal direction')
ylabel(app.psiBoundaryPlot,'\psi')
title(app.psiBoundaryPlot,'Check \Psi_b contributions');
legend(app.psiBoundaryPlot,'\psi_{plasma}','\psi_{coils}','\psi_{coils}+\psi_{plasma}','\psi_{boundary}','location','eastoutside')
%    set(gca,'FontWeight','Normal','Fontname','times','Fontsize',18)
%% calcolo matrice green attivi->boundary
% load sources
% coils.current=coils.current;
% coils.current=ones(56,1);
% coils.type='thick';

% disp('calculating GG_solenoid_space');
% tic
% GG_coils=solenoidCalc(coils, space); %Green matrix coils-boundary
% toc

%% CONNESSIONI tra coils
% figure;
% for ii=1:length(coils.R)
% plot(coils.R(ii),coils.Z(ii),'rX')
% hold on
% plot(coils.DR,coils.DZ,'o')
% grid on
% axis equal
% title(num2str(ii))
% pause
% end

% KONNAX = [
% 1 zeros(1,4) 1 zeros(1,8); 
% 0 1 zeros(1,2) 1 zeros(1,9);
% 0 0 1 1  zeros(1,10);
% zeros(1,6) 1 zeros(1,4) 1 0 0;
% zeros(1,7) 1 zeros(1,2) 1 0 0 0;
% zeros(1,8) 1 1 zeros(1,4);
% zeros(1,12) 1 0;
% zeros(1,13) 1];
% KONNAX = eye(14,14);
if (isfield(coils,'kon'))
           KONNAX = coils.kon;
else
KONNAX = eye(12,12);
end

%calcolo routine matteo (piu veloce)
tic
GG_coils = fun_Field_Coil( coils, space, 10 );
toc

GG_coils_k=GG_coils.psi*KONNAX';
%save  ([dir_save 'horizontal_elongated_inverse_data_equilibrium.mat'],'GG_coils_k', 'psi_coils','plasma')

%% SOLUZIONE INVERSO

Green_coils = GG_coils_k;

Psi_coils = psi_coils;


% %%
% k_G = cond(Green_coils);
% if k_G > 10^3
%     disp('problema malcondizionato K(Green):')
%     disp(num2str(k_G))
% end
% 
% disp('rango Green:')
% disp(num2str(rank(Green_coils)))
% if (rank(Green_coils) == min(size(Green_coils)))
% disp('matrice Green ha rango pieno')
% end
% 
% 
%% calcolo matrice green attivi->spazio occupato dal plasma e plasma->spaziooccupato dal plasma
%check mappa flusso su plasma intero
space2.RR=p(:,1);
space2.ZZ=p(:,2);
tic
GG_coils_nodes = fun_Field_Coil( coils, space2, 10 );
toc
tic
GG_plasma_nodes= fun_Field_Loop(plasma, space2);
toc


%global G_boundary PSI_coils_boundary
G_boundary = Green_coils;
PSI_coils_boundary = Psi_coils;
% save([dir_save  'OPTIM_DATA_INPUT'], 'G_boundary', 'PSI_coils_boundary', 'psi_tot_direct','p','t','R_boundary_finale','Z_boundary_finale',...
%     'GG_plasma_nodes','plasma','GG_coils_nodes','KONNAX','psiBoundary','psi_ak'); 

OptimDataIn.G_boundary = G_boundary;
OptimDataIn.PSI_coils_boundary=PSI_coils_boundary;
OptimDataIn.psi_tot_direct =psi_tot_direct; 
OptimDataIn.p = p;
OptimDataIn.t = t;
OptimDataIn.R_boundary_finale = app.OutputEquil.OutputShape.R_boundary_finale;
OptimDataIn.Z_boundary_finale = app.OutputEquil.OutputShape.R_boundary_finale;

OptimDataIn.GG_plasma_nodes = GG_plasma_nodes;
OptimDataIn.plasma = plasma; 
OptimDataIn.GG_coils_nodes = GG_coils_nodes;
OptimDataIn.kon = KONNAX;
OptimDataIn.psiBoundary = psiBoundary;
OptimDataIn.psiAxis = app.OutputEquil.Equil.Psi_axis;

OptimDataIn.RR = app.OutputEquil.Equil.RR;
OptimDataIn.ZZ = app.OutputEquil.Equil.ZZ;
end

