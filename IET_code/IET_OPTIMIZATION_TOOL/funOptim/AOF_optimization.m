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
function [AOF_solution] = AOF_optimization(app)
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

Psi_axis = app.OutputEquil.Equil.Psi_axis;
Psi_boundary = app.OutputEquil.Equil.Psi_boundary;

%% COMMENTO 2: lambda global
n_pareto = app.nParetoPoints;
lambda = linspace(0,1,n_pareto);

[SOF_solution] = LS_optimization(app);
x_least = SOF_solution.currents;
[MAPE_LS,RMSE_LS] = scalarRMSEMAPE_Results(x_least',G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);

parfor ii=1:n_pareto
 options=optimset('MaxFunEvals',100000,'MaxIter',10000,'algorithm','interior-point');
%options=optimset('MaxFunEvals',100000,'MaxIter',10000);
[x_lambda_global(ii,:),f_lambda_global(ii)] = fmincon(@(W) scalarRMSEMAPE(W,lambda(ii),G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary),x_least',[],[],[],[],v_lb,v_ub,[],options);

end


for ii=1:n_pareto
    [MAPE_lambda(ii),RMSE_lambda(ii)] = scalarRMSEMAPE_Results(x_lambda_global(ii,:),G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);

% Psi_tot_Lambda_Global(ii,:)=GG_plasma_nodes.psi*plasma.current+...
%     GG_coils_nodes.psi*KONNAX'*x_lambda_global(ii,:)';

end

 [~,ivmin] =min(f_lambda_global);
% x_lambda_best = x_lambda_global(ivmin,:)';
% Psi_tot_LambdaBest=GG_plasma_nodes.psi*plasma.current+...
%     GG_coils_nodes.psi*KONNAX'*x_lambda_best;
% MAPE_lambda(ivmin)
% RMSE_lambda(ivmin)

F_sol = [MAPE_lambda' RMSE_lambda'];
[PFront, ~]=DominanceFilter(F_sol,x_lambda_global);

% [ParetoFronte,isort]= sort(PFront(:,1),1);

% pareto front
cla(app.ParetoFront)
plot(app.ParetoFront,MAPE_lambda,RMSE_lambda,'bo','markersize',13,'linewidth',1.5)
hold(app.ParetoFront,'on');
plot(app.ParetoFront,PFront(:,1),PFront(:,2),'m*','linewidth',1.5);
grid(app.ParetoFront,'on');
axis(app.ParetoFront,'equal');

plot(app.ParetoFront,MAPE_lambda(ivmin),RMSE_lambda(ivmin),'rx','markersize',22,'linewidth',2)
plot(app.ParetoFront,MAPE_LS,RMSE_LS,'gp','markersize',14,'markerfacecolor','g')
xlabel(app.ParetoFront,'MAPE','interpreter','latex')
ylabel(app.ParetoFront,'RMSE','interpreter','latex')

legend(app.ParetoFront,'AOF solutions','Pareto Front','AOF min','LS','location','best')
set(app.NormalizedDesignSpace,'FontWeight','Normal','Fontname','times','Fontsize',12)

AOF_solution.ParetoFrontData = PFront;
AOF_solution.MAPE_lambda = MAPE_lambda;
AOF_solution.RMSE_lambda = RMSE_lambda;
AOF_solution.currents = x_lambda_global;

end