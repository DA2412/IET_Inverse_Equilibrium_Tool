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
function [NC_solution] = NC_optimization(app)
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

%% MAPE mu1 Anchor point no. 1
% Lancio le ottimizzazioni in un dominio normalizzato 0-1. La conversione
% nel dominio effettivo avviene nelle funzioni onlymu1_norm e onlymu2_norm
dof = size(G_boundary,2);
i_startingPoints = app.initPoints;
pop_starting = lhsdesign(i_startingPoints,dof);
% pop_start_points = bsxfun(@plus,v_lb,bsxfun(@times,pop_starting,v_ub-v_lb));

%Run the optimizer
W_bestmu1 = zeros(i_startingPoints,dof);
f_bestmu1 = zeros(i_startingPoints,1);
W_bestmu2 = zeros(i_startingPoints,dof);
f_bestmu2 = zeros(i_startingPoints,1);

options=optimset('MaxFunEvals',100000,'MaxIter',10000,'Algorithm','sqp');%matlabpool open 2
% options=optimset('LargeScale','on','TolFun',.00001,'MaxIter',100000,'MaxFunEvals',100000);%matlabpool open 2
%options =[];
tic
parfor countRun = 1:i_startingPoints
    
    xpop=pop_starting(countRun,:)';
    
    [W_bestSetmu1,f_bestSetmu1,exitflag,output,lambda,grad,hessian] = fmincon(@(W_norm)onlymu1(W_norm,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary),xpop,[],[],[],[],zeros(dof,1),ones(dof,1),[],options);
    
    W_bestmu1(countRun,:) = W_bestSetmu1;
    f_bestmu1(countRun,:) = f_bestSetmu1;
    
    [W_bestSetmu2,f_bestSetmu2,exitflag,output,lambda,grad,hessian] = fmincon(@(W_norm)onlymu2(W_norm,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary),xpop,[],[],[],[],zeros(dof,1),ones(dof,1),[],options);
    
    W_bestmu2(countRun,:) = W_bestSetmu2;
    f_bestmu2(countRun,:) = f_bestSetmu2;
    
end
toc

[mu1star,ivmin_mu1] =min(f_bestmu1);
W1star_norm= W_bestmu1(ivmin_mu1,:)';


[mu2star,ivmin_mu2] =min(f_bestmu2);
W2star_norm = W_bestmu2(ivmin_mu2,:)';


%% STEP2 mappings
Utopia = [mu1star mu2star];

[mu12star,~] = onlymu1(W2star_norm,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);
[mu21star,~] = onlymu2(W1star_norm,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);

L1 = mu12star - mu1star;
L2 = mu21star - mu2star;

L = [L1,L2];

%% step3  utopia line vector

mubar1star = muNorm(W1star_norm,Utopia,L,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);
mubar2star = muNorm(W2star_norm,Utopia,L,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);
N1 = mubar2star - mubar1star;

%% step 4 normalized increments
m1 = app.mSolUtopia ; %nmber of solutions

delta1 = 1/ (m1-1);

%% step 5 generate utopia line points
lambda = 0:delta1:1;
lambda = flip(lambda);
alfa1j = lambda;
alfa2j = 1-lambda;


% Xpj = bsxfun(@times,[alfa1j;alfa1j],mubar1star) + bsxfun(@times,[alfa2j;alfa2j],mubar2star);
Xpj = bsxfun(@times,[alfa1j;alfa1j],mubar1star) + bsxfun(@times,[alfa2j;alfa2j],mubar2star);


%% step 6
% options = [];
options=optimset('MaxFunEvals',100000,'MaxIter',10000,'Display','iter-detailed','Algorithm','sqp');%matlabpool open 2

clc
% diary('Fmincon_results.dat')
for jj=1:m1
    % parfor jj=m1:-1:1
    
    if jj==1
        x0 = W1star_norm; %x_least;%ones(12,1)*10e3;
    else
        x0 = x_lambda_global(:,jj-1);
    end
    [x_lambda_global(:,jj),f_lambda_global(jj),exitflag,output,lambda,grad,hessian] = fmincon(@(W_norm)mu2Norm(W_norm,Utopia,L,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary),...
        x0,[],[],[],[],...
        zeros(dof,1),ones(dof,1),...
        @(W_norm)SQPnolcon(W_norm,Utopia,L,N1,Xpj(:,jj),v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary),...
        options);
    %  x0,[],[],[],[],v_lb,v_ub,@(W)SQPnolcon(W,Utopia,L,N1,Xpj(:,jj)),options);
    
end
% diary off

for kk=1:m1
    [mu1_realbest(kk),x_best_NC_real(:,kk)] = onlymu1(x_lambda_global(:,kk),v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);
    [mu2_realbest(kk),] = onlymu2(x_lambda_global(:,kk),v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);
    
    
    mubarstar_best(:,kk) = muNorm(x_lambda_global(:,kk),Utopia,L,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);
end

plot(app.NormalizedDesignSpace,mubar1star(1),mubar1star(2),'mo','markersize',20,'linewidth',1.5);
hold(app.NormalizedDesignSpace,'on');
grid(app.NormalizedDesignSpace,'on');
plot(app.NormalizedDesignSpace,mubar2star(1),mubar2star(2),'gx','markersize',20,'linewidth',1.5);
plot(app.NormalizedDesignSpace,Xpj(1,2:end),Xpj(2,2:end),'r+','linewidth',1.5)
% plot(app.NormalizedDesignSpace,Xpj(1,1),Xpj(2,1),'db','markersize',12,'linewidth',1.5)
% plot(Utopia(1,1),Utopia(1,2),'ks','markersize',20,'linewidth',1.5);
plot(app.NormalizedDesignSpace,mubarstar_best(1,1),mubarstar_best(2,1),'rd','markersize',20,'linewidth',2)
plot(app.NormalizedDesignSpace,mubarstar_best(1,:),mubarstar_best(2,:),'rp','markersize',13,'linewidth',2)
legend(app.NormalizedDesignSpace,'\mu1*','\mu2*','Utopia line','Pareto Front 1st point','Pareto Front other points')
set(app.NormalizedDesignSpace,'FontWeight','Normal','Fontname','times','Fontsize',12)
axis(app.NormalizedDesignSpace,'equal');

plot(app.GeneralDesignSpace,mu1_realbest,mu2_realbest,'bo','markersize',13,'linewidth',1)
hold(app.GeneralDesignSpace,'on');
grid(app.GeneralDesignSpace,'on');
plot(app.GeneralDesignSpace,mu1star,mu21star,'md','markersize',16,'linewidth',2)
plot(app.GeneralDesignSpace,mu12star,mu2star,'rd','markersize',16,'linewidth',2)
xlabel(app.GeneralDesignSpace,'MAPE')
ylabel(app.GeneralDesignSpace,'RMSE')
legend(app.GeneralDesignSpace,'Pareto Points','MAPE Anchor Point','RMSE ANchor Point')
set(app.GeneralDesignSpace,'FontWeight','Normal','Fontname','times','Fontsize',12)
axis(app.GeneralDesignSpace,'square');

NC_solution.mu1_realbest = mu1_realbest;
NC_solution.mu2_realbest = mu2_realbest;
NC_solution.currents = x_best_NC_real;
end
 
