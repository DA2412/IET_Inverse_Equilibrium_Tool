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
function [SOF_solution] = LS_solution(app)
G_boundary = app.OptimDataIn.G_boundary;
PSI_coils_boundary = app.OptimDataIn.PSI_coils_boundary;
GG_plasma_nodes = app.OptimDataIn.GG_plasma_nodes;
plasma = app.OptimDataIn.plasma;
KONNAX= app.OptimDataIn.kon;
GG_coils_nodes = app.OptimDataIn.GG_coils_nodes;
v_ub = app.OptimDataIn.UpperBounds;
v_lb = app.OptimDataIn.LowerBounds;
Psi_tot_direct = app.OutputEquil.Equil.psi_tot_direct;
RR = app.OutputEquil.Equil.RR;
ZZ = app.OutputEquil.Equil.ZZ;

Aeq = sparse(size(G_boundary,2),size(G_boundary,2));
beq = zeros(size(G_boundary,2),1);

C = G_boundary;
d = PSI_coils_boundary;

options = optimoptions('lsqlin','OptimalityTolerance', 1e-9,'ConstraintTolerance',1e-12,'MaxIterations',800);
[x_least,resnorm,residual,exitflag,output,lambda]= lsqlin(C,d,[],[],Aeq,beq,v_lb,v_ub,[],options);

Psi_tot_LS=GG_plasma_nodes.psi*plasma.current+...
    GG_coils_nodes.psi*KONNAX'*x_least;

rr = app.OutputEquil.OutputShape.node(:,1);
zz = app.OutputEquil.OutputShape.node(:,2);
PSI_LS = griddata(rr,zz,Psi_tot_LS,RR,ZZ);

if app.OutputEquil.Equil.psiBoundary <0 
    Psi_const = max(Psi_tot_direct);
    err_patch_LS = (Psi_tot_LS- Psi_tot_direct)./(Psi_tot_direct+Psi_const)*100;
else
    err_patch_LS = (Psi_tot_LS- Psi_tot_direct)./(Psi_tot_direct)*100;
end

%plot contour
contour(app.PsiInvPlot,RR,ZZ,PSI_LS,10,'linewidth',2)
hold(app.PsiInvPlot,'on');
plot(app.PsiInvPlot,app.OutputEquil.OutputShape.R_boundary_finale,app.OutputEquil.OutputShape.Z_boundary_finale,'-r','linewidth',2);
grid(app.PsiInvPlot,'on');
axis(app.PsiInvPlot,'equal');
colormap(app.PsiInvPlot,jet);
colorbar(app.PsiInvPlot,'vert');
xlabel(app.PsiInvPlot,'R [m]')
ylabel(app.PsiInvPlot,'Z [m]')
title(app.PsiInvPlot,'Flux surfaces')
ylim(app.PsiInvPlot,[min(zz) max(zz)]) 
xlim(app.PsiInvPlot,[min(rr) max(rr)]) 
set(app.PsiInvPlot,'FontWeight','Normal','Fontname','times','Fontsize',18)

%plot error
p = app.OutputEquil.OutputShape.node;
t= app.OutputEquil.OutputShape.tri;
JP.faces=t;
JP.vertices=p;
JP.facevertexcdata=err_patch_LS;
patch(app.PsiErrorPlot,JP,'facecolor','interp','edgecolor','none');
hold(app.PsiErrorPlot,'on');
plot(app.PsiErrorPlot,app.OutputEquil.OutputShape.R_boundary_finale,app.OutputEquil.OutputShape.Z_boundary_finale,'-r','linewidth',2);
grid(app.PsiErrorPlot,'on');
axis(app.PsiErrorPlot,'equal');
colormap(app.PsiErrorPlot,jet);
hb=colorbar(app.PsiErrorPlot,'vert','FontWeight','Normal','Fontname','times','Fontsize',18);
title(hb,'\epsilon_{\psi} [%]');
xlabel(app.PsiErrorPlot,'R [m]')
ylabel(app.PsiErrorPlot,'Z [m]')
title(app.PsiErrorPlot,'Error map')
ylim(app.PsiErrorPlot,[min(zz) max(zz)]) 
xlim(app.PsiErrorPlot,[min(rr) max(rr)]) 
set(app.PsiErrorPlot,'FontWeight','Normal','Fontname','times','Fontsize',18)

%plot correnti
plot(app.CurrentsPlot,1:length(v_ub),x_least,'mo','linewidth',2,'markersize',15)
hold(app.CurrentsPlot,'on');
grid(app.CurrentsPlot,'on');
plot(app.CurrentsPlot,1:length(v_ub),v_ub,'r-','linewidth',2)
plot(app.CurrentsPlot,1:length(v_ub),v_lb,'r-','linewidth',2)
ylabel(app.CurrentsPlot,'I [A]')
xlabel(app.CurrentsPlot,'# coil')
xticks(app.CurrentsPlot,1:length(v_ub))
xlim(app.CurrentsPlot,[1 length(v_ub)])
title(app.CurrentsPlot,'LS currents solution');
set(app.CurrentsPlot,'FontWeight','Normal','Fontname','times','Fontsize',24)

%output
SOF_solution.currents = x_least;
SOF_solution.Psi_LS = Psi_tot_LS;
SOF_solution.err_patch_LS  = err_patch_LS;

end

