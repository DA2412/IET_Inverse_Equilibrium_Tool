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

% renaming vars
psi_tot_direct = OutputEquil.Equil.psi_tot_direct;
Iphi_plasma = OutputEquil.Equil.Iphi_plasma;
psiBoundary = OutputEquil.Equil.psiBoundary;
rr = OutputEquil.OutputShape.p(:,1);
zz = OutputEquil.OutputShape.p(:,2);
p =  OutputEquil.OutputShape.p;
t=  OutputEquil.OutputShape.t;
Area_duale =  OutputEquil.OutputShape.Area_duale;
nodi_BCs =  OutputEquil.OutputShape.b;

% Triangle point indices
it1=t(:,1);
it2=t(:,2);
it3=t(:,3);
% Find centroids of triangles
rbar=(rr(it1)+rr(it2)+rr(it3))/3;
zbar=(zz(it1)+zz(it2)+zz(it3))/3;

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

psi_plasma = GG_plasma.psi * plasma.current;
psi_coils = psi_tot_direct(nodi_BCs) - psi_plasma;    %calculating coils contribution to psi boundary

%% COILS - calcolo matrice green attivi->boundary
if (isfield(coils,'kon'))
           KONNAX = coils.kon;
else
           KONNAX = eye(12,12);
end
tic
GG_coils = fun_Field_Coil( coils, space, 10 );
toc

%% first wall
rfw = 0.4905;
theta = linspace(0,2*pi,100);
R_fw = 1.995 + rfw*cos(theta);
Z_fw = rfw*sin(theta);

%% Maximum Constraints
nMW = 4;
nFSW = 10;
ImaxMW = 50e3;
ImaxFSW = 6.25e3;
v_ub = [ImaxMW*ones(1,nMW),  ImaxFSW*ones(1,4), 0, ImaxFSW, 0, ImaxFSW*ones(1,3)];
v_lb = [-ImaxMW*ones(1,nMW), -ImaxFSW*ones(1,4), 0, -ImaxFSW, 0, -ImaxFSW*ones(1,3)];

% %% LS solution: flux on boundary
GG_coils_flux_b = GG_coils.psi*KONNAX';
% 
Aeq = sparse(size(GG_coils_flux_b,2),size(GG_coils_flux_b,2));
Aeq (9,9) = 1; %FS4D
Aeq (11,11) = 1;%FS5D
beq = zeros(size(GG_coils_flux_b,2),1);

C = GG_coils_flux_b;
d = psi_coils;
options = optimoptions('lsqlin','OptimalityTolerance', 1e-9,'ConstraintTolerance',1e-9,'MaxIterations',800);
[x_least,resnorm,residual,exitflag,output,lambda]= lsqlin(C,d,[],[],Aeq,beq,v_lb,v_ub,[],options);

% plot
npoint = 200;
rmax = 1.995+0.5;
rmin = 1.995-0.5;
zmin = 0-0.5;
zmax = 0+0.5;
rgrid=linspace(rmin,rmax,npoint);
zgrid=linspace(zmin,zmax,npoint);
[RR_g,ZZ_g] = meshgrid(rgrid,zgrid);

space_g.RR=RR_g(:);
space_g.ZZ=ZZ_g(:);
tic
GG_coils_grid = fun_Field_Coil( coils, space_g, 10 );
toc
tic
GG_plasma_grid= fun_Field_Loop(plasma, space_g);
toc

Psi_LS = GG_coils_grid.psi*KONNAX'*x_least + GG_plasma_grid.psi*plasma.current;
BR_LS  = GG_coils_grid.br*KONNAX'*x_least  + GG_plasma_grid.br*plasma.current;
BZ_LS  = GG_coils_grid.bz*KONNAX'*x_least  + GG_plasma_grid.bz*plasma.current;

PSI_LS = reshape(Psi_LS,size(RR_g));
BR_grid = reshape(BR_LS,size(RR_g));
BZ_grid = reshape(BZ_LS,size(RR_g));

Bpol = sqrt(BR_grid.^2+BZ_grid.^2);
[Bsort,ibsort] = sort(Bpol(:));
psi_b = PSI_LS(ibsort(2));

dim=1/100;